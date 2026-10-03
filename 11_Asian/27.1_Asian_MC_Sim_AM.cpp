#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include "../../Headers/black_scholes_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto monteCarloAsianCall(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, int simulations, int steps) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    commonType payoffs {};
    commonType deltat {data.maturity / static_cast<commonType>(paths)};
    commonType sqrt_deltat {std::sqrt(deltat)};

    for (int i {}; i < simulations; ++i)
    {
        commonType S_avg {};
        commonType S {data.spot};
        for (int j {}; j < steps; ++j)
        {
            S *= std::exp((data.rate - 0.5 * data.volatility * data.volatility) * deltat + data.volatility * W(mt) * sqrt_deltat);
            S_avg += S;
        }
        S_avg /= static_cast<commonType>(steps);
        payoffs += std::max(S_avg - data.strike, static_cast<commonType>(0));
    }
    payoffs /= static_cast<commonType>(simulations);
    
    return std::exp(-data.rate * data.maturity) * payoffs;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto monteCarloAsianPut(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, int simulations, int steps) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    commonType payoffs {};
    commonType deltat {data.maturity / static_cast<commonType>(paths)};
    commonType sqrt_deltat {std::sqrt(deltat)};

    for (int i {}; i < simulations; ++i)
    {
        commonType S_avg {};
        commonType S {data.spot};
        for (int j {}; j < steps; ++j)
        {
            S *= std::exp((data.rate - 0.5 * data.volatility * data.volatility) * deltat + data.volatility * W(mt) * sqrt_deltat);
            S_avg += S;
        }
        S_avg /= static_cast<commonType>(steps);
        payoffs += std::max(data.strike - S_avg, static_cast<commonType>(0));
    }
    payoffs /= static_cast<commonType>(simulations);
    
    return std::exp(-data.rate * data.maturity) * payoffs;
}

int main()
{
    int simulations {100000};
    int steps {1000};
    OptionDataBS myOption {100, 100, 0.01, 0.3, 1};
    std::cout << "Call price: $" << monteCarloAsianCall(myOption, simulations, steps) << '\n';
    std::cout << "Put price: $" << monteCarloAsianPut(myOption, simulations, steps) << '\n';

    return 0;
}