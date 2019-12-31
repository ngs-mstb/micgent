from __future__ import unicode_literals
from __future__ import print_function
from __future__ import division
from __future__ import absolute_import

from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import *

import pandas as pd
import numpy as np


def series_equal(a,b):
    return a.sort_values().reset_index(drop=True).equals(
        b.sort_values().reset_index(drop=True))

def keys_equal(a,b,keys=None):
    if not keys:
        if not list(a) == list(b):
            return False
        keys = list(a)
    keys_a = keys
    keys_b = keys
    return a[keys_a].sort_values(keys_a).reset_index(drop=True).equals(
        b[keys_b].sort_values(keys_b).reset_index(drop=True))

def create_empty_df(columns):
    return pd.DataFrame(columns=columns)

def create_columns_in_df(df,columns):
    """Add columns to an existing dataframe, using specified dtypes and constant scalar values.
    This is similar to SQL `alter table` and works on empty data frames too.

    :param df: DataFrame, can be empty or not
    :param columns: dictionary { column name => scalar value }. The dtype of the new column
    will be deduced by pandas from teh supplied scalar value. To specify a particular dtype,
    one can create the value like `numpy.dtype('int32').type(0)`. 

    If adding columns to the empty data frame, the scalar values themselves will be ignored
    (no rows are ever appended here), and only be used to define the dtype of the new columns.
    """
    nrows = df.shape[0]
    for column in columns:
        df[column] = columns[column]
    assert df.shape[0] == nrows, \
        "Detected implementation change in Pandas - adding columns changed number of rows from {} to {}".\
        format(nrows,df.shape[0])
    return df

def add_column_name_prefix(x,pref,columns=None):
    if columns is None:
        columns = list(x.columns)
    return x.rename(columns=dict(((col,"{}{}".format(pref,col)) for col in columns)))
