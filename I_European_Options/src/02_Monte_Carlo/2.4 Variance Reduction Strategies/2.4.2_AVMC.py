import numpy as np
import AVMC
import pandas as pd
import os

print(os.getcwd())

def trial(nIterations, nRepetitions):
    simulations = np.logspace(2, 6, nIterations, dtype=int)
    simulations += simulations % 2
    prices_calls = np.zeros(nIterations)
    standardErrors_calls = np.zeros(nIterations)
    prices_puts = np.zeros(nIterations)
    standardErrors_puts = np.zeros(nIterations)
    for i in range(0,len(simulations)):
        for j in range(0, nRepetitions):
            calls = np.array(AVMC.monteCarloBSCallAV(100.0, 100.0, 0.05, 0.2, 1.0, simulations[i], j))
            puts = np.array(AVMC.monteCarloBSPutAV(100.0, 100.0, 0.05, 0.2, 1.0, simulations[i], j))
            prices_calls[i] += np.average(calls)
            prices_puts[i] += np.average(puts)
            avg_call_pairs = np.empty(simulations[i] // 2)
            avg_put_pairs = np.empty(simulations[i] // 2)
            for k in range(0, simulations[i], 2):
                avg_put_pairs[k // 2] = 0.5 * (puts[k] + puts[k+1])
                avg_call_pairs[k //2] = 0.5 * (calls[k] + calls[k+1])
            standardErrors_calls[i] += np.std(avg_call_pairs, ddof=1) / np.sqrt(avg_call_pairs.size)
            standardErrors_puts[i] += np.std(avg_put_pairs, ddof=1) / np.sqrt(avg_put_pairs.size)
        prices_calls[i] /= nRepetitions
        prices_puts[i] /= nRepetitions
        standardErrors_calls[i] /= nRepetitions
        standardErrors_puts[i] /= nRepetitions
    return simulations, prices_calls, prices_puts, standardErrors_calls,  standardErrors_puts

simulations, callPrices, putPrices, STDcallErrors, STDputErrors = trial(200, 25)

df = pd.DataFrame(
    {
        "# Simulations": simulations,
        "Call Prices": callPrices,
        "Put Prices": putPrices,
        "Standard Error of Calls": STDcallErrors,
        "Standard Error of Puts": STDputErrors
    }
)

df.to_csv("2.4.2_results.csv", index=False)
print("hi")
