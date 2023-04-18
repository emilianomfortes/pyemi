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
