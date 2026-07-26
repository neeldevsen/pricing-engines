#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include "../Headers/black_scholes_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto monteCarloBSCall(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, int simulations) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), std::sqrt(data.maturity)};
    commonType values {};

    std::generate(values.begin(), values.end(), 
    [&] () 
    {
        double S_final {data.spot * std::exp((data.rate - 0.5 * data.volatility * data.volatility) * data.maturity + (data.volatility * W(mt)))};
        return std::exp(-data.rate * data.maturity) * std::max(S_final - data.strike, static_cast<double>(0));
    } )
    


    for (int i {}; i < simulations; ++i)
    {
        commonType S_final {data.spot * std::exp((data.rate - 0.5 * data.volatility * data.volatility) * data.maturity + (data.volatility * W(mt)))};
        payoffs += std::max(S_final - data.strike, static_cast<commonType>(0));
    }
    payoffs /= simulations;
    
    return std::exp(-data.rate * data.maturity) * payoffs;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto monteCarloBSPut(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, int simulations) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), std::sqrt(data.maturity)};
    commonType payoffs {};

    for (int i {}; i < simulations; ++i)
    {
        commonType S_final {data.spot * std::exp((data.rate - 0.5 * data.volatility * data.volatility) * data.maturity + (data.volatility * W(mt)))};
        payoffs += std::max(data.strike - S_final , static_cast<commonType>(0));
    }
    payoffs /= simulations;

    return std::exp(-data.rate * data.maturity) * payoffs;
}


int main()
{
    int simulations {10000000};
    OptionDataBS myOption {100, 100, 0.01, 0.3, 30};
    std::cout << "Call price: $" << monteCarloBSCall(myOption, simulations) << '\n';
    std::cout << "Call price: $" << monteCarloBSPut(myOption, simulations) << '\n';

    return 0;
}
