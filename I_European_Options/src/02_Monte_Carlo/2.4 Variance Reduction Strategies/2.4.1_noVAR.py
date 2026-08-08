import numpy as np
import noRedMC 
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
            calls = np.array(noRedMC.monteCarloBSCall(100.0, 100.0, 0.05, 0.2, 1.0, simulations[i], j))
            puts = np.array(noRedMC.monteCarloBSPut(100.0, 100.0, 0.05, 0.2, 1.0, simulations[i], j))
            prices_calls[i] += np.average(calls)
            prices_puts[i] += np.average(puts)
            standardErrors_calls[i] += np.std(calls, ddof=1) / np.sqrt(simulations[i])
            standardErrors_puts[i] += np.std(puts, ddof=1) / np.sqrt(simulations[i])
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

df.to_csv("2.4.1_results.csv", index=False)
print("hi")