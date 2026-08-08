#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include "../../Headers/functions.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

std::vector<double> monteCarloBSCall(double spot, double strike, double rate, double volatility, double maturity, int simulations, unsigned long seed) 
{
    std::mt19937 mt{seed};
    std::normal_distribution W {static_cast<double>(0), std::sqrt(maturity)};
    std::vector<double> prices (simulations);

    for (int i {}; i < simulations; ++i)
    {
        double S_final {spot * std::exp((rate - 0.5 * volatility * volatility) * maturity + (volatility * W(mt)))};
        prices[i] = std::exp(-rate * maturity) * max(S_final - strike, static_cast<double>(0));
    }
    
    return prices;
}

std::vector<double> monteCarloBSPut(double spot, double strike, double rate, double volatility, double maturity, int simulations, unsigned long seed)
{
    std::mt19937 mt{seed};
    std::normal_distribution W {static_cast<double>(0), std::sqrt(maturity)};
    std::vector<double> prices (simulations);

    for (int i {}; i < simulations; ++i)
    {
        double S_final {spot * std::exp((rate - 0.5 * volatility * volatility) * maturity + (volatility * W(mt)))};
        prices[i] = std::exp(-rate * maturity) * max(strike - S_final, static_cast<double>(0));
    }
    
    return prices;
}

PYBIND11_MODULE(noRedMC, module)
{
    module.def("monteCarloBSCall", &monteCarloBSCall);
    module.def("monteCarloBSPut", &monteCarloBSPut);
}

