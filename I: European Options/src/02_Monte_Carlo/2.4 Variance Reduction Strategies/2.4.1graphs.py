import matplotlib.pyplot as plt
import pandas as pd 
import numpy as np

df = pd.read_csv("2.4.1_results.csv")

simulations = df["# Simulations"]
call_prices = df["Call Prices"]
put_prices = df["Put Prices"]
call_errors = df["Standard Error of Calls"]
put_errors = df["Standard Error of Puts"]
abs_errors_calls = df["ABS Error Calls"]
abs_errors_puts = df["ABS Error Puts"]

fig, ax = plt.subplots(figsize=(10, 6))

ax.plot(
    simulations,
    abs_errors_calls,
    color="red",
    linestyle="-.",
    marker="x",
    markersize=1.5,
    linewidth=0.5,
    label="Monte Carlo Call"
)

ax.plot(
    simulations,
    abs_errors_puts,
    color="blue",
    linestyle="-.",
    marker="x",
    markersize=1.5,
    linewidth=0.7,
    label="Monte Carlo Put"
)

ax.set_xscale("log")
ax.set_xlabel("Number of simulations")
ax.set_ylabel("Absolute Error ($)")
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