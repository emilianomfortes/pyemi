import numpy as np
import pandas as pd


def cumprod_returns(data: (pd.DataFrame or pd.Series or np.ndarray), col=None):
    cumprod_r = None
    if isinstance(data, pd.DataFrame):
        cumprod_r = (data[col] + 1.0).cumprod() - 1.0
    elif isinstance(data, pd.Series):
        cumprod_r = (data + 1.0).cumprod() - 1.0
    elif isinstance(data, np.ndarray):
        cumprod_r = (data + 1.0).cumprod() - 1.0
    return cumprod_r


# Volatility
def ATR(df: pd.DataFrame, period: int = 14):
    """Average True Range"""
    _df = df.copy()

    true_range = np.zeros(len(df))
    true_range[0] = df["High"].iloc[0] - df["Low"].iloc[0]

    for i in range(1, len(df)):
        true_range[i] = max(
            df["High"].iloc[i] - df["Low"].iloc[i],
            abs(df["High"].iloc[i] - df["Close"].iloc[i - 1]),
            abs(df["Low"].iloc[i] - df["Close"].iloc[i - 1]),
        )

    _df["TR"] = true_range
    _df["ATR"] = _df["TR"].rolling(window=period).mean()
    return _df["ATR"]
