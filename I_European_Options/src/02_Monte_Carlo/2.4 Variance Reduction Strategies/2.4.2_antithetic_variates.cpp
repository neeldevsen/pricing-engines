#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "../../Headers/functions.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

std::vector<double> monteCarloBSCallAV(double spot, double strike, double rate, double volatility, double maturity, int simulations, unsigned long seed) 
{
    std::mt19937 mt{seed};
    std::normal_distribution W {static_cast<double>(0), std::sqrt(maturity)};
    std::vector<double> prices (simulations);

    for (int i {}; i < simulations / 2; ++i)
    {
        double Z = W(mt);
        double S_final_1 {spot * std::exp((rate - 0.5 * volatility * volatility) * maturity + (volatility * Z))};
        prices[2*i] = std::exp(-rate * maturity) * max(S_final_1 - strike, static_cast<double>(0));
        double S_final_2 {spot * std::exp((rate - 0.5 * volatility * volatility) * maturity - (volatility * Z))};
        prices[2*i+1] = std::exp(-rate * maturity) * max(S_final_2 - strike, static_cast<double>(0));
    }
    
    return prices;
}

std::vector<double> monteCarloBSPutAV(double spot, double strike, double rate, double volatility, double maturity, int simulations, unsigned long seed) 
{
    std::mt19937 mt{seed};
    std::normal_distribution W {static_cast<double>(0), std::sqrt(maturity)};
    std::vector<double> prices (simulations);

    for (int i {}; i < simulations / 2; ++i)
    {
        double Z = W(mt);
        double S_final_1 {spot * std::exp((rate - 0.5 * volatility * volatility) * maturity + (volatility * Z))};
        prices[2*i] = std::exp(-rate * maturity) * max(strike - S_final_1, static_cast<double>(0));
        double S_final_2 {spot * std::exp((rate - 0.5 * volatility * volatility) * maturity - (volatility * Z))};
        prices[2*i+1] = std::exp(-rate * maturity) * max(strike - S_final_2, static_cast<double>(0));
    }

    return prices;
}

PYBIND11_MODULE(AVMC, module)
{
    module.def("monteCarloBSCallAV", &monteCarloBSCallAV);
    module.def("monteCarloBSPutAV", &monteCarloBSPutAV);
}

