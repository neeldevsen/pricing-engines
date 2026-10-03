#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include "../../Headers/black76_structs.hpp"

template <typename F_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto monteCarloB76Call(const OptionDataB76<F_t, K_t, r_t, sigma_t, T_t>& data, int simulations) -> decltype(data.future)
{
    using commonType = decltype(data.future);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), std::sqrt(data.maturity)};
    commonType payoffs {};

    for (int i {}; i < simulations; ++i)
    {
        commonType F_final {data.future * std::exp(((data.volatility * W(mt)) - 0.5 * data.volatility * data.volatility * data.maturity))};
        payoffs += std::max(F_final - data.strike, static_cast<commonType>(0));
    }
    payoffs /= simulations;
    
    return std::exp(-data.rate * data.maturity) * payoffs;
}

template <typename F_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto monteCarloB76Put(const OptionDataB76<F_t, K_t, r_t, sigma_t, T_t>& data, int simulations) -> decltype(data.future)
{
    using commonType = decltype(data.future);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), std::sqrt(data.maturity)};
    commonType payoffs {};

    for (int i {}; i < simulations; ++i)
    {
        commonType F_final {data.future * std::exp(((data.volatility * W(mt)) - 0.5 * data.volatility * data.volatility * data.maturity))};
        payoffs += std::max(data.strike - F_final, static_cast<commonType>(0));
    }
    payoffs /= simulations;
    
    return std::exp(-data.rate * data.maturity) * payoffs;
}


int main()
{
    int simulations {1000000};
    OptionDataB76 myOption {80, 100, 0.01, 0.3, 30};
    std::cout << "Call price: $" << monteCarloB76Call(myOption, simulations) << '\n';
    std::cout << "Put price: $" << monteCarloB76Put(myOption, simulations) << '\n';

    return 0;
}
