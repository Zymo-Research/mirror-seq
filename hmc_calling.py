BISMARK_METH_CODE_TYPE_MAP = {
    "H": "CHH",
    "h": "CHH",
    "X": "CHG",
    "x": "CHG",
    "Z": "CpG",
    "z": "CpG",
}
CALLING_TYPES = (
    'normal',
    'mirror',
)

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
        result_df.sort(['chrom', 'pos', 'strand'], inplace=True)
        if start!=None:
            idx_start = result_df['pos'].searchsorted(start, 'left')
            result_df = result_df[idx_start:]
        if end!=None:
            idx_end = result_df['pos'].searchsorted(end, 'right')
            result_df = result_df[:idx_end]
        result_df['chrom'] = result_df['chrom'].replace(tid_chrom_d)
    return result_df

def write_meth_data_by_regions(bam_filename, out_dir, regions, calling_type='normal', rand_str=''):
    ''' Write the region methylation calling DataFrame into a HDF5 file.

    Parameters
    ----------
    bam_filename : str
        The alignment BAM filename.
    out_dir : str
        The ouput directory of HDF5 file.
    regions : List of tuples
        A list of (chromosome, start, end).
    calling_type : str, optional
        If it is "mirror", the function call "mirror_seq_conversion" function for
        CpG DataFrame.
    rand_str : str, optional
        Add the rand_str in the prefix to tempfiles.

    Notes
    -----
    * The output file is a temp file, and you have to delete it manually.
    '''

    import pandas as pd
    import tempfile
    import pysam

    with pysam.AlignmentFile(bam_filename) as samfile:
        max_chrom_len = max([len(d['SN']) for d in samfile.header['SQ']])

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
            if calling_type=='mirror' and meth_code=='Z':
                mirror_seq_conversion(tmp_df)

            prefix = 'tmp_{0}_'.format(rand_str)
            suffix = '_{0}'.format(meth_type)
            f = tempfile.NamedTemporaryFile(dir=out_dir, prefix=prefix, suffix=suffix, delete=False)
            tmp_df.to_hdf(
                f.name,
                'sites',
                data_columns=True,
                format='t',
                complevel=5,
                complib='blosc',
                append=True,
                min_itemsize={'chrom': max_chrom_len},
            )

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

def parse_to_bed(hdf_filename, bed_filename, table_name='sites', chunksize=1000000):
    ''' Parse the standard output format (HDF5) to conventional sqlite3 format.

    Parameters
    ----------
    hdf_filename : str
        The HDF5 filename.
    bed_filename : str
        The output BED filename.
    table_name : str, optional
        The table name of the dataframe.
    chunksize : int, optional
        The chunk size per read_hdf.
    '''
    import pandas as pd
    import numpy as np
    import subprocess
    import os

    with open(bed_filename, 'w') as fw:
        for df in pd.read_hdf(hdf_filename, table_name, chunksize=chunksize):
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

def mirror_seq_conversion(df):
    import pandas as pd

    pos_series = pd.concat((
        df[df['strand']=='+']['pos'] + 1,
        df[df['strand']=='-']['pos'] - 1,
    ))
    df['pos'] = pos_series
    df['strand'] = df['strand'].replace({'+': '--', '-': '++'}).replace({'--': '-', '++': '+'})
    df['meth_count'] = df['total_count'] - df['meth_count']

def merge_n_parse(out_prefix, meth_type, hdf_filenames, max_chrom_len, create_bed_file=True):
    ''' The is a shortcut function, which is easier to be used by multiprocessing.

    Parameters
    ----------
    h5_filename : str
        The HDF5 filename.
    '''
    import pandas as pd
    import os

    full_hdf_filename = '{0}_{1}.h5'.format(out_prefix, meth_type)
    for hdf_filename in hdf_filenames:
        df = pd.read_hdf(hdf_filename, 'sites')
        df.to_hdf(
            full_hdf_filename,
            'sites',
            data_columns=True,
            format='t',
            complevel=5,
            complib='blosc',
            append=True,
            min_itemsize={'chrom': max_chrom_len},
        )
        del df

    if create_bed_file:
        bed_filename = os.path.splitext(full_hdf_filename)[0]+'.bed'
        parse_to_bed(full_hdf_filename, bed_filename)

def get_chrom_sizes_file(bam_filename, chrom_sizes_filename):
    ''' Generate the chrom.sizes file from bam file.
    bam_filename : str
        The alignment BAM file.
    '''
    import pysam

    with pysam.AlignmentFile(bam_filename) as samfile, open(chrom_sizes_filename, 'w') as fw:
        for d in samfile.header['SQ']:
            fw.write('{0}\t{1}\n'.format(d['SN'], d['LN']))

def get_bs_conv_rate(filenames):
    '''Calculate the bisulfite conversion rate using CHH and CHG methylation tracks.

    Parameters
    ----------
    filenames : List of str
        the HDF5 filenames of non-CpGs.

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

    meth_ratio_sum = 0
    count = 0

    for filename in filenames:
        for df in pd.read_hdf(filename, 'sites', columns=['meth_count', 'total_count'], chunksize=1000000):
            meth_ratio_sum += (df['meth_count'] / df['total_count']).sum()
            count += len(df)

    try:
        bs_conv_rate = 1 - round(meth_ratio_sum / count, 2)
    except ZeroDivisionError:
        bs_conv_rate = None

    return bs_conv_rate


def main(bam_filename, calling_type, out_prefix, nts_in_regions=100000000):
    ''' Run the entire methylation calling.

    Parameters
    ----------
    bam_filename : str
        The alignment bam filename. The index file (.bai) must exist in the same folder.
    out_prefix : str
        The output file prefix. The output file is <out_prefix>_<METH_TYPE>.h5.
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

    out_dir = os.path.dirname(out_prefix)
    rand_str = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))

    chrom_sizes_filename = 'chrom.sizes'
    get_chrom_sizes_file(bam_filename, chrom_sizes_filename)
    with pysam.AlignmentFile(bam_filename) as samfile:
        max_chrom_len = max([len(d['SN']) for d in samfile.header['SQ']])

    p = Pool()
    for regions in get_regions_chunks(bam_filename, nts_in_regions):
        p.apply_async(
            write_meth_data_by_regions,
            (bam_filename, out_dir, regions, calling_type, rand_str),
        )
    p.close()
    p.join()

    prefix = 'tmp_{0}_'.format(rand_str)
    meth_type_filenames_dict = {}
    for filename in os.listdir(out_dir):
        if not filename.startswith(prefix):
            continue
        meth_type = filename.rsplit('_', 1)[-1]
        filenames = meth_type_filenames_dict.setdefault(meth_type, [])
        filenames.append(os.path.join(out_dir, filename))

    p = Pool(min(len(meth_type_filenames_dict), multiprocessing.cpu_count()))
    for meth_type, filenames in meth_type_filenames_dict.iteritems():
        p.apply_async(
            merge_n_parse,
            (out_prefix, meth_type, filenames, max_chrom_len, True),
        )

    p.close()
    p.join()

    # bedToBigBed seems use a lot of memory. Do it one by one.
    for meth_type in meth_type_filenames_dict:
        bed_filename = '{0}_{1}.bed'.format(out_prefix, meth_type)
        if os.path.exists(bed_filename):
            with pd.HDFStore('{0}_{1}.h5'.format(out_prefix, meth_type)) as store:
                row_count = store.root.sites.table.nrows
            if row_count < 550000000:
                subprocess.check_call((
                    'bedToBigBed',
                    bed_filename,
                    chrom_sizes_filename,
                    '{0}_{1}.bb'.format(out_prefix, meth_type),
                ))
            else:
                subprocess.check_call(('pigz', bed_filename))

    # Remove tmp files after everthing is done.
    for hdf_filenames in meth_type_filenames_dict.itervalues():
        for hdf_filename in hdf_filenames:
            os.remove(hdf_filename)

if __name__=='__main__':
    import argparse
    import os

    parser = argparse.ArgumentParser(
        description='Methylation calling for Methyl-seq alignment files.'
    )
    parser.add_argument(
        '-b',
        dest='bam_filename',
        required=True,
        help='The BAM filename.'
    )
    parser.add_argument(
        '-t',
        dest='calling_type',
        default=CALLING_TYPES[0],
        choices=CALLING_TYPES,
        help='''Specify "mirror" for mirror-seq alignment files,
        otherwise use the default (normal).'''
    )
    parser.add_argument(
        '--nts-in-regions',
        dest='nts_in_regions',
        default=100000000,
        type=int,
        help='''Number of total nucleotides in an iter of regions.
        It is an rough number so it is possible to get more than the number. default is 100M.'''
    )
    parser.add_argument(
        '-o',
        dest='out_prefix',
        default=None,
        help='''The output preifx of all output files. With absolute path is recommended.
        If None (default), use the pathname of bam file without extension.'''
    )
    args = parser.parse_args()

    if args.out_prefix:
        out_prefix = args.out_prefix
    else:
        out_prefix = os.path.splitext(args.bam_filename)[0]
    main(
        args.bam_filename,
        args.calling_type,
        out_prefix,
        args.nts_in_regions,
    )
