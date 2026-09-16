<<<<<<< HEAD
For parts 1-5 please refer to this [PDF1](./I_Classical_Models_and_Methods/PDF/PricingEngines.pdf) for the .1 subparts and this [PDF2](./I_Classical_Models_and_Methods/PDF/SupplementEngines.pdf) for the .2 and .3 subparts.
=======
For parts 1-5 please refer to this [PDF1](./PricingEngines.pdf) for the x.1 subparts and this [PDF2](./SupplementEngines.pdf) for the x.2 and x.3 subparts.
>>>>>>> 2d784cef6d24102b5a801bfd82b7904bd8588539

So far implemented:
1. Black Scholes Closed Forms
2. Monte Carlo Simulations
3. Binomial Trees
4. Trinomial Trees
5. Finite Difference Methods
6. [Heston's Model](./II_Advanced_Models/06_Heston/Heston.md)
7. Merton's Jump Diffusion Model
- [Convergence Study](./ConvergenceStudy.pdf)
- Calibration 

The directory structure is shown below:
```text

pricing-engines/
├── README.md
├── LICENSE
├── .gitignore
├── PricingEngines.pdf
├── SupplementEngines.pdf
├── ConvegenceStudy.pdf
│
├── Headers/
│
├── I_Classical_Models_and_Methods/
│   ├── 01_Closed_Forms/
│   ├── 02_Monte_Carlo/
│   ├── 03_Binomial_Trees/
│   ├── 04_Trinomial_Trees/
│   ├── 05_Finite_Differences/
│
├── II_Advanced_Models/
│   ├── 06_Heston/
│   ├── 07_Jump_Diffusion/
│
├── Calibration/
│
├── I_Convergence_Study/
│   ├── C++/
│   ├── Data-Analysis/
│   ├── Excels/
│   ├── Experiments/
│   ├── Images/




```
## Build

This project uses C++20.

Compile a source file with:

```bash
g++ -std=c++20 file.cpp -o exec
```


