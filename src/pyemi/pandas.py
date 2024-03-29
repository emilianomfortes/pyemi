import pandas as pd

def drop_duplicate_columns(df):
	"""
	Drop duplicate columns in a pandas dataframe. 
	"""
	# https://stackoverflow.com/questions/14984119/python-pandas-remove-duplicate-columns

	df_clean = df.loc[:,~df.columns.duplicated()].copy()
	return df_clean

def check_for_nas_in_df(df):
	answer = df.isnull().values.any()
	return answer


def monotone_df_cols(df, cols, increasing=True):
	"""Returns a boolean mask for a dataframe where the columns are monotonic."""
	if increasing:
		mon_inc = (df[cols].cummax().diff().fillna(.1) > 0).all(axis=1)
	else:
		mon_inc = (df[cols].cummin().diff().fillna(.1) < 0).all(axis=1)
	df2 = df[mon_inc]
	return df2


#
# DateTime
#
def find_last_date_of_month_changes(df):

    # Check that the index is a datetime index
    if not isinstance(df.index, pd.DatetimeIndex):
        raise TypeError("Index must be a datetime index")

    # Get the last date of each month
    last_dates = df.resample("M").last().index

    return last_dates

def find_first_date_of_month_changes(df):

    # Check that the index is a datetime index
    if not isinstance(df.index, pd.DatetimeIndex):
        raise TypeError("Index must be a datetime index")

    # Get the last date of each month
    last_dates = df.resample("M").first().index
    
    first_dates = []
    first_dates.append(df.index[0])
    for last_date in last_dates[:-1]:
        new_df = df.loc[last_date:]
        first_date = new_df.index[1]
        first_dates.append(first_date)
    
    try:
        new_df = df.loc[last_dates[-1]:]
        first_date = new_df.index[1]
        first_dates.append(first_date)
    except IndexError:
        pass

    first_dates = pd.DatetimeIndex(first_dates)

    return first_dates 
