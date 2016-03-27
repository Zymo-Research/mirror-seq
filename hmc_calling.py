def mirror_seq_conversion(df):
    ''' Convert methylation ratios and switch to the other strands.

    Parameters
    ----------
    df : Pandas.DataFrame
        The dataframe must has columns - pos, strand, meth_count, and total_count.

    Returns
    -------
    Pandas.DataFrame
        The converted dataframe.
    '''
    import pandas as pd

    pos_series = pd.concat((
        df[df['strand']=='+']['pos'] + 1,
        df[df['strand']=='-']['pos'] - 1,
    ))
    del df['pos']
    df['pos'] = pos_series
    df['strand'] = df['strand'].replace({'+': '--', '-': '++'}).replace({'--': '-', '++': '+'})
    df['meth_count'] = df['total_count'] - df['meth_count']

    return df
