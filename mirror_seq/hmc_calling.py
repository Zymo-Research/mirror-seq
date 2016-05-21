BISMARK_METH_CODE_TYPE_MAP = {
    "H": "CHH",
    "h": "CHH",
    "X": "CHG",
    "x": "CHG",
    "Z": "CpG",
    "z": "CpG",
}

def meth_call_for_read(read, overlap=True, min_qual=20):
    ''' Do methyltaion calling per read.

    Parameters
    ----------
    read : pysam.AlignedSegment
        The read from fetch()
    overlap : bool
        If it is True, only count the sites in read 2 for overlapped region.
    min_qual : int
        The minium quality score to do methyltaion calling.

    yields
    -------
    int
        reference ID
    int
        position
    str
        strand
    str
        methylation code

    Notes
    -----
    * Because of technical problem, call read 2 sites for the overlapped sites
    instead of read 1 sites.
    * If the overlapped read 2 sites do not pass min_qual, it won't call read 1
    sites instead.
    '''

    if read.is_paired:
        if read.flag in (99, 147):
            strand = '+'
        else:
            strand = '-'
    else:
        if read.is_reverse:
            strand = '-'
        else:
            strand = '+'

    positions = read.get_reference_positions()
    is_left = read.reference_start < read.next_reference_start
    meth_codes = read.get_tag('XM')
    quals = read.query_qualities
    for pos, meth_code, qual in zip(positions, meth_codes, quals):
        if meth_code in BISMARK_METH_CODE_TYPE_MAP and qual>=min_qual:
            if read.is_paired and overlap and pos>=read.next_reference_start and is_left:
                continue
            yield read.reference_id, pos, strand, meth_code

def meth_call_by_region(bam_filename, chrom=None, start=None, end=None):
    ''' Methylation call for a given region.

    Parameters
    ----------
    bam_filename : str
        The alignment BAM filename.
    chrom : str, optional
        The chromsome name of the region.
    start : int, optional
        The start position of the region.
    end : int, optional
        The end position of the region.

    Returns
    -------
    pandas.DataFrame
        columns is ['chrom', 'pos', 'strand', 'meth_code', 'meth_count', 'total_count'].
    '''

    import pandas as pd
    import pysam
    import numpy as np

    print 'Working on {}:{}-{}'.format(chrom, start, end)
    with pysam.AlignmentFile(bam_filename) as samfile:
        tid_chrom_d = {i: d['SN'] for i, d in enumerate(samfile.header['SQ'])}
        # Values are tuples of meth_count and totoal counts.
        coor_meth_calls_d = {}
        for read in samfile.fetch(chrom, start, end):
            for reference_id, pos, strand, meth_code in meth_call_for_read(read):
                meth_calls = coor_meth_calls_d.setdefault(
                    (reference_id, pos, strand, meth_code.upper()),
                    [0, 0]
                )
                if meth_code.isupper():
                    meth_calls[0] += 1
                meth_calls[1] += 1

    result_df = pd.DataFrame()
    if coor_meth_calls_d:
        result_df = pd.DataFrame(
            coor_meth_calls_d.values(),
            columns=['meth_count', 'total_count'],
            index=pd.MultiIndex.from_tuples(
                coor_meth_calls_d.keys(),
                names=['chrom', 'pos', 'strand', 'meth_code']
            ),
            dtype=np.uint32,
        )
        result_df.reset_index(inplace=True)
        result_df['chrom'] = result_df['chrom'].astype(np.uint32)
        result_df.sort_values(['chrom', 'pos', 'strand'], inplace=True)
        if start!=None:
            idx_start = result_df['pos'].searchsorted(start, 'left')[0]
            result_df = result_df[idx_start:]
        if end!=None:
            idx_end = result_df['pos'].searchsorted(end, 'right')[0]
            result_df = result_df[:idx_end]
        result_df['chrom'] = result_df['chrom'].replace(tid_chrom_d)
    return result_df

def write_meth_data_by_regions(bam_filename, out_dir, regions, rand_str=''):
    ''' Write the region methylation calling DataFrame into a file.

    Parameters
    ----------
    bam_filename : str
        The alignment BAM filename.
    out_dir : str
        The ouput directory.
    regions : List of tuples
        A list of (chromosome, start, end).
    rand_str : str, optional
        Add the rand_str in the prefix to tempfiles.

    Notes
    -----
    * The output file is a temp file, and you have to delete it manually.
    '''

    import pandas as pd
    import tempfile
    import pysam

    result_df = pd.DataFrame()
    for chrom, start, end in regions:
        result_df = result_df.append(meth_call_by_region(bam_filename, chrom, start, end))

    if not result_df.empty:
        for meth_code in result_df['meth_code'].unique():
            meth_type = BISMARK_METH_CODE_TYPE_MAP[meth_code]
            tmp_df = result_df[result_df['meth_code']==meth_code]
            tmp_df = tmp_df[['chrom', 'pos', 'strand', 'meth_count', 'total_count']]
            tmp_df = tmp_df.reset_index(drop=True)
            # Mirror-seq can only detect CpGs so do not convert non-CpGs.
            if meth_code=='Z':
                mirror_seq_conversion(tmp_df)

            prefix = 'tmp_{0}_'.format(rand_str)
            suffix = '_{0}'.format(meth_type)
            f = tempfile.NamedTemporaryFile(dir=out_dir, prefix=prefix, suffix=suffix, delete=False)
            tmp_df.to_csv(f.name, compression='gzip', index=False)

def get_regions_chunks(bam_filename, nts_in_regions=100000000):
    ''' Iterate regions lists to roughly fit "nts_in_regions".

    Parameters
    ----------
    bam_filename : str
        The alignment BAM filename.
    nts_in_regions : int, optional
        Number of total nucleotides in an iter of regions. It is an rough number
        so it is possible to get more than the number.

    Yields
    ------
    List
        List of tuples of chromosome, start, and end.
    '''
    import pysam

    with pysam.AlignmentFile(bam_filename) as samfile:
        chrom_sizes = [(d['SN'], d['LN']) for d in samfile.header['SQ']]
        chrom_sizes = sorted(chrom_sizes, key=lambda x: x[1], reverse=True)

    regions = []
    nts = 0
    for chrom, size in chrom_sizes:
        chunk_num = size / nts_in_regions
        if size % nts_in_regions != 0:
            chunk_num += 1

        start = 0
        for i in range(chunk_num):
            end = min(start + nts_in_regions, size)
            region_size = end - start
            nts += region_size
            regions.append((chrom, start, end))
            if nts>nts_in_regions:
                yield regions
                nts = 0
                regions = []
            start = end + 1
    if regions:
        yield regions

def parse_to_bed(data_filename, bed_filename, chunksize=1000000):
    ''' Parse the standard output format to BED format.

    Parameters
    ----------
    data_filename : str
        The data filename.
    bed_filename : str
        The output BED filename.
    chunksize : int, optional
        The chunk size when reading files.
    '''
    import pandas as pd
    import numpy as np
    import subprocess
    import os

    with open(bed_filename, 'w') as fw:
        for df in pd.read_csv(data_filename, chunksize=chunksize):
            df['end'] = df['pos'] + 1
            df['thick_start'] = 0
            df['thick_end'] = 0
            df['meth_ratio'] = df['meth_count'] / df['total_count']
            df['score'] = (df['meth_ratio'] * 1000).round().astype(np.uint32)
            df['name'] = df.apply(
                lambda r: '{0}/{1}({2:.0%})'.format(
                    r['meth_count'],
                    r['total_count'],
                    r['meth_ratio']
                ),
                axis=1,
            )
            df['rgb'] = df.apply(
                lambda r: '255,{0:.0f},0'.format(255*r['meth_ratio']),
                axis=1,
            )

            colnames = [
                'chrom',
                'pos',
                'end',
                'name',
                'score',
                'strand',
                'thick_start',
                'thick_end',
                'rgb'
            ]
            df[colnames].to_csv(fw, sep='\t', index=False, header=False)
    temp_folder = os.path.dirname(bed_filename)
    subprocess.check_call((
        "sort",
        "-T", temp_folder,
        "-k1,1",
        "-k2,2n",
        bed_filename,
        "-o", bed_filename,
    ))
    subprocess.check_call(('gzip', '-f', bed_filename))

def mirror_seq_conversion(df):
    ''' Convert methylation ratios and strands.

    Parameters
    ----------
    df : pandas.DataFrame
        with columns - strand, pos, meth_count, and total_count.
    '''
    import pandas as pd

    pos_series = pd.concat((
        df[df['strand']=='+']['pos'] + 1,
        df[df['strand']=='-']['pos'] - 1,
    ))
    df['pos'] = pos_series
    df['strand'] = df['strand'].replace({'+': '--', '-': '++'}).replace({'--': '-', '++': '+'})
    df['meth_count'] = df['total_count'] - df['meth_count']

def merge_n_parse(out_prefix, meth_type, filenames, create_bed_file):
    ''' The is a shortcut function, which is easier to be used by multiprocessing.

    Parameters
    ----------
    out_prefix : str
        The output prefix.
    meth_type : str
        The methylation type. Eg: CpG, CHG, and CHH.
    filenames : str
        The csv filenames to be mreged.
    create_bed_file : bool
        Create a bed file or not.
    '''
    import pandas as pd
    import os
    import subprocess

    full_filename = '{0}_{1}.csv'.format(out_prefix, meth_type)
    header = True
    for filename in filenames:
        pd.read_csv(filename, compression='gzip').to_csv(full_filename, header=header, mode='a', index=False)
        header = False
    subprocess.check_call(('gzip', '-f', full_filename))
    full_filename += '.gz'

    if meth_type=='CpG' and create_bed_file:
        bed_filename = full_filename.replace('.csv.gz', '.bed')
        parse_to_bed(full_filename, bed_filename)

def get_bs_conv_rate(filenames):
    '''Calculate the bisulfite conversion rate using CHH and CHG methylation tracks.

    Parameters
    ----------
    filenames : List of str
        the filenames of non-CpGs.

    Returns
    -------
    float
        The estimated bisulfite conversion rate.

    NOTES
    -----
    1. Get conversion rate for each non CpGs and average them.
    2. The esitmated bisulfite conversion rate is rounded to 2 decimal.
    3. Return None if no sites in all files.

    '''
    import pandas as pd
    import os

    meth_ratio_sum = 0
    count = 0

    for filename in filenames:
        if os.path.exists(filename):
            for df in pd.read_csv(filename, compression='gzip', usecols=['meth_count', 'total_count'], chunksize=1000000):
                meth_ratio_sum += (df['meth_count'] / df['total_count']).sum()
                count += len(df)

    try:
        bs_conv_rate = 1 - round(meth_ratio_sum / count, 2)
    except ZeroDivisionError:
        bs_conv_rate = None

    return bs_conv_rate


def main(bam_filename, out_prefix, create_bed_file, nts_in_regions=100000000):
    ''' Run the entire methylation calling.

    Parameters
    ----------
    bam_filename : str
        The alignment bam filename. The index file (.bai) must exist in the same folder.
    out_prefix : str
        The output file prefix. The output file is <out_prefix>_<METH_TYPE>.h5.
    create_bed_file : bool
        Create a bed file or not.
    nts_in_regions : int, optional
        Number of total nucleotides in an iter of regions. It is an rough number
        so it is possible to get more than the number.

    '''
    from multiprocessing import Pool
    import multiprocessing
    import subprocess
    import pysam
    import os, string, random
    import pandas as pd

    print('Wokring on hydroxymethylation calling...')
    out_dir = os.path.dirname(out_prefix)
    rand_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))

    p = Pool()
    for regions in get_regions_chunks(bam_filename, nts_in_regions):
        p.apply_async(
            write_meth_data_by_regions,
            (bam_filename, out_dir, regions, rand_str),
        )
    p.close()
    p.join()

    prefix = 'tmp_{0}_'.format(rand_str)
    meth_type_filenames_dict = {}

    for filename in os.listdir(os.path.join('.', out_dir)):
        if not filename.startswith(prefix):
            continue
        meth_type = filename.rsplit('_', 1)[-1]
        filenames = meth_type_filenames_dict.setdefault(meth_type, [])
        filenames.append(os.path.join(out_dir, filename))

    print('Merge files...')

    p = Pool()
    for meth_type, filenames in meth_type_filenames_dict.iteritems():
        p.apply_async(
            merge_n_parse,
            (out_prefix, meth_type, filenames, True),
        )
    p.close()
    p.join()

    cpg_filename = '{}_CpG.csv.gz'.format(out_prefix)
    chg_filename = '{}_CHG.csv.gz'.format(out_prefix)
    chh_filename = '{}_CHH.csv.gz'.format(out_prefix)
    # Calculate bisulfite conversion rate.
    conversion_rate = get_bs_conv_rate([
        chg_filename,
        chh_filename,
    ])
    if conversion_rate is not None:
        print('Bisuflite conversion rate: {:.0%}'.format(conversion_rate))
    else:
        print('Cannot estimate bisuflite conversion rate.')
    # Remove tmp files after everthing is done.
    for filenames in meth_type_filenames_dict.itervalues():
        for filename in filenames:
            os.remove(filename)

    try:
        os.remove(chg_filename)
    except OSError:
        pass
    try:
        os.remove(chh_filename)
    except OSError:
        pass

    print('Done!')
