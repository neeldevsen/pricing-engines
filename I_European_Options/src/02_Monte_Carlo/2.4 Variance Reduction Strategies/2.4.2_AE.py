import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

dfNR = pd.read_csv("2.4.1_results.csv")
dfAV = pd.read_csv("2.4.2_results.csv")

simulationsNR = dfNR["# Simulations"]
call_pricesNR = dfNR["Call Prices"]
put_pricesNR = dfNR["Put Prices"]
call_errorsNR = dfNR["Standard Error of Calls"]
put_errorsNR = dfNR["Standard Error of Puts"]
abs_errors_callsNR = dfNR["ABS Error Calls"]
abs_errors_putsNR = dfNR["ABS Error Puts"]

simulationsAV = dfAV["# Simulations"]
call_pricesAV = dfAV["Call Prices"]
put_pricesAV = dfAV["Put Prices"]
call_errorsAV = dfAV["Standard Error of Calls"]
put_errorsAV = dfAV["Standard Error of Puts"]
abs_errors_callsAV = dfAV["ABS Error Calls"]
abs_errors_putsAV = dfAV["ABS Error Puts"]

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(
    simulationsAV,
    abs_errors_callsAV,
    color="red",
    linestyle="-.",
    marker="x",
    markersize=1.5,
    linewidth=0.5,
    label="Monte Carlo Call with Antithetic Variates"
)

ax.plot(
    simulationsAV,
    abs_errors_callsNR,
    color="blue",
    linestyle="-.",
    marker="x",
    markersize=1.5,
    linewidth=0.7,
    label="Vanilla Monte Carlo Call"
)

ax.set_xscale("log")
ax.set_xlabel("Number of simulations")
ax.set_ylabel("Absolute Error($)")
ax.set_yscale("log")
ax.set_title("Absolute Error vs Simulation Size")

ax.grid(
    True,
    which="both",
    color="0.82",
    linestyle=":",
    linewidth=0.6
)

ax.grid(True, which="major", linestyle=":", linewidth=0.6)
ax.legend(frameon=False)

fig.tight_layout()
fig.savefig("2.4.1_price_convergence.pdf", bbox_inches="tight")
fig.savefig("2.4.1_price_convergence.png", dpi=300, bbox_inches="tight")

plt.show()