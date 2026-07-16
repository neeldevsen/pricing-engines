#include <iostream>
#include <cmath>
#include <type_traits>
#include "../Headers/functions.hpp"
#include "../Headers/black_scholes_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto blackScholesCall(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data) -> decltype(data.spot)
{
    using CommonType = decltype(data.spot);
    auto sqrt_maturity {std::sqrt(data.maturity)};

    auto d1 {(std::log(data.spot / data.strike) + (data.rate + (static_cast<CommonType>(0.5) * data.volatility * data.volatility)) * (data.maturity))/ (data.volatility * sqrt_maturity)};
    auto d2 {d1 - data.volatility * sqrt_maturity};

    auto C {data.spot * normalCDF(d1) - data.strike * std::exp(-data.rate * data.maturity) * normalCDF(d2)};
    return C;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto blackScholesPut(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data) -> decltype(data.spot)
{
    using CommonType = decltype(data.spot);
    auto sqrt_maturity {std::sqrt(data.maturity)};
    
    auto d1 {(std::log(data.spot / data.strike) + (data.rate + (static_cast<CommonType>(0.5) * data.volatility * data.volatility)) * (data.maturity)) / (data.volatility * sqrt_maturity)};
    auto d2 {d1 - data.volatility * sqrt_maturity};

    auto P  {data.strike * std::exp(-data.rate * data.maturity) * normalCDF(-d2) - data.spot * normalCDF(-d1)};
    return P;
}

int main()
{
    OptionDataBS myOption {100, 100, 0.01, 0.3, 30};
    std::cout << "Call price: $" << blackScholesCall(myOption) << '\n';
    std::cout << "Put price: $" << blackScholesPut(myOption) << '\n';

    return 0;  
}







