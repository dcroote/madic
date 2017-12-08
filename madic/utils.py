from __future__ import division
import numpy as np


def norm_col_sum_gb(df, column, grouping_cols=None):
    """ Normalizes column to sum of one, optionally grouping by columns

    Necessary for transition ratio calculations for each peptide. Also, handles
    interference by masking column values.

    Args:
        df (pd.DataFrame): containing column and grouping columns
        column (str): column to normalize
        grouping_cols (list): columns to group by when normalizing

    Returns:
        (pd.Series) normalized column
    """

    if 'interference' in df.columns:
        df = df.copy()
        # mask areas with nan where interference is True
        df.loc[df[df.interference].index, column] = np.nan

    if grouping_cols is None:
        series_norm = df[column] / df[column].sum()
    else:
        series_norm = df[column] / df.groupby(
            grouping_cols)[column].transform('sum')

    series_norm.name = '{}_frac'.format(column)

    return series_norm


def cols_to_category_dtype(df, columns=None):
    """ Converts list of pandas DataFrame columns to category dtype """

    if columns is not None:
        for col in columns:
            if not hasattr(df[col], 'cat'):
                df[col] = df[col].astype('category')

    return df
