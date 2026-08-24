from scipy.stats import norm
import numpy as np
from scipy import optimize
import pandas as pd
import yfinance as yf
import matplotlib.pyplot as plt
from pathlib import Path 
import jump_diffusion as jd
from datetime import datetime, timezone

string = input("What ticker do you want? \n")
ticker = yf.Ticker(string)
expiry_dates = np.array(ticker.options)

sizes = np.zeros(len(expiry_dates))
for i in range(0, len(expiry_dates)):
    call_data = ticker.option_chain(expiry_dates[i]).calls
    sizes[i] = len(call_data["strike"])

maturities = np.zeros(20)
now = datetime.now(timezone.utc)

for i in range(0,len(maturities)):
    maturities[i] = (datetime.strptime(expiry_dates[i], "%Y-%m-%d").replace(tzinfo=timezone.utc) - now).total_seconds() / 31536000

strikes_grid = np.full((len(expiry_dates), 20), np.nan)
prices_grid = np.full((len(expiry_dates), 20), np.nan)

for i in range(0, len(expiry_dates)):
    data = ticker.option_chain(expiry_dates[i]).calls
    data = data[(data["ask"] > 0) & (data["bid"] > 0)] 
    prices = ((data["ask"] + data["bid"]) / 2).to_numpy()
    n = min(len(prices), 20)
    prices_grid[i, 0:n] = prices[0:n]
    strikes = (data["strike"]).to_numpy()
    strikes_grid[i, 0:n] = strikes[0:n]
    
def objective_call(params, S, K, r, T, N, prices):
    ji, mlj, jsv = params
    error = 0
    for i in range(0,len(T)):
        if T[i] < 0:
            continue
        for j in range(0, len(K[i])):
            if not(np.isnan(K[i][j]) or np.isnan(prices[i][j])):
                model_price = (jd.call(S, K[i][j], r, T[i], ji, mlj, jsv, N))
                error += (model_price - prices[i][j]) ** 2
    return error

spot = float(ticker.history(period="1d")["Close"].iloc[-1])
interest_rate = 0.0363
N = 5000

bounds = [
    (1e-4, 5.0),   
    (-0.5, 0.5),  
    (1e-4, 0.8)    
]

result = optimize.minimize(
    objective_call,
    x0=[0.5, -0.1, 0.2],
    args=(spot,strikes_grid,interest_rate,maturities,N, prices_grid),
    method="L-BFGS-B",
    bounds=bounds
)

print(result)