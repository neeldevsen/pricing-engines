Please refer to the PDF for all of the derivations and code explanations
So far implemented:
- Black Scholes Closed Forms
- Monte Carlo Simulations
- Binomial Trees
- Trinomial Trees
- Finite Difference Methods

```text
The directory structure is shown below:
pricing-engines/
├── README.md
├── LICENSE
├── .gitignore
├── PricingEngines.pdf
├── SupplementEngines.pdf
│
├── I_Classical_Models_and_Methods/
│   ├── 01_Closed_Forms/
│   ├── 02_Monte_Carlo/
│   ├── 03_Binomial_Trees/
│   ├── 04_Trinomial_Trees/
│   ├── 05_Finite_Differences/
│   ├── Headers/
```
## Build

This project uses C++20.

Compile a source file with:

```bash
g++ -std=c++20 file.cpp -o exec
