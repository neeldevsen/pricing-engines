For parts 1-5 please refer to this [PDF1](./I_Classical_Models_and_Methods/PDF/PricingEngines.pdf) for the .1 subparts and this [PDF2](./I_Classical_Models_and_Methods/PDF/SupplementEngines.pdf) for the .2 and .3 subparts.
=======
For parts 1-5 please refer to this [PDF1](./I_Classical_Models_and_Methods/PDF/PricingEngines.pdf) for the x.1 subparts and this [PDF2](./I_Classical_Models_and_Methods/PDF/SupplementEngines.pdf) for the x.2 and x.3 subparts.

P.S: there is a reason there is missing numbers, I don't create these in order but I have an overall order in mind and I stick to it

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

17. American Binomial Tree
18. American Trinomial Tree

- And then some other pricing engines such as exotic options:
  
26. Barrier Options
27. Asian Options
28. Digital Options
    
The directory structure is shown below:
```text

pricing-engines/
├── README.md
├── LICENSE
├── .gitignore
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
│   ├── PDF/
│   |  ├── PricingEngines.pdf
│   |  ├── SupplementEngines.pdf
|
├── II_Advanced_Models/
│   ├── 06_Heston/
│   ├── 07_Jump_Diffusion/
|
├── IV_American_Options/
│   ├── 17_American_Binomial_Tree/
│   ├── 18_American_Trinomial_Tree/
|
├── VI_Exotic_Options/
│   ├── 28_Barrier_Options/
│   ├── 29_Asian_Options/
│   ├── 30_Digital_Options/
|
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


