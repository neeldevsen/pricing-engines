from scipy.stats import norm
import numpy as np
from scipy import optimize
import pandas as pd
import yfinance as yf
import matplotlib.pyplot as plt
from pathlib import Path 
import heston_qe as hqe
from datetime import datetime, timezone

string = input("What ticker do you want? \n")
ticker = yf.Ticker(string)
expiry_dates = np.array(ticker.options)[:20]

sizes = np.zeros(len(expiry_dates))
for i in range(0, len(expiry_dates)):
    call_data = ticker.option_chain(expiry_dates[i]).calls
    sizes[i] = len(call_data["strike"])

maturities = np.zeros(len(expiry_dates))
now = datetime.now(timezone.utc)

for i in range(0,len(maturities)):
    maturities[i] = (datetime.strptime(expiry_dates[i], "%Y-%m-%d").replace(tzinfo=timezone.utc) - now).total_seconds() / 31536000

strikes_grid = np.full((len(expiry_dates), 20), np.nan)
prices_grid = np.full((len(expiry_dates), 20), np.nan)

for i in range(0, len(expiry_dates)):
    data = ticker.option_chain(expiry_dates[i]).calls
    data = data[(data["ask"] > 0) & (data["bid"] > 0) & (data["ask"] >= data["bid"])] 
    prices = ((data["ask"] + data["bid"]) / 2).to_numpy()
    n = min(20, len(prices))
    prices_grid[i, 0:n] = prices[0:n]
    strikes = (data["strike"]).to_numpy()
    strikes_grid[i, 0:n] = strikes[0:n]
    
def objective_call(params, S, K, r, T, M, N, threshold, gamma1, gamma2, prices):
    kappa, theta, xi, rho, v = params
    error = 0
    for i in range(0,len(T)):
        if T[i] <= 0:
            continue
        for j in range(0, len(K[i])):
            if not(np.isnan(K[i][j]) or np.isnan(prices[i][j])):
                model_price = (hqe.call(S, K[i][j], r, T[i], v, kappa, theta, xi, rho, M, N, threshold, gamma1, gamma2))
                error += (model_price - prices[i][j]) ** 2
    return error

spot = float(ticker.history(period="1d")["Close"].iloc[-1])
interest_rate = 0.0363
threshold = 1.5
gamma1 = 0.5
gamma2 = 0.5
M = 5000
N = 100

bounds = [
    (0.01, 10.0),   
    (0.001, 1.0),  
    (0.01, 5.0),    
    (-0.999, 0.999),
    (0.001, 1.0)    
]

result = optimize.minimize(
    objective_call,
    x0=[2.0, 0.04, 0.5, -0.7, 0.04],
    args=(spot,strikes_grid,interest_rate,maturities,M,N,threshold,gamma1,gamma2, prices_grid),
    method="L-BFGS-B",
    bounds=bounds
)

print(result)