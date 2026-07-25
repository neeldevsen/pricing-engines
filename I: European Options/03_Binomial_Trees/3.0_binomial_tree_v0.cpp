#include <iostream>
#include <cmath>
#include <type_traits>
#include <algorithm>
#include "../Headers/black_scholes_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto binomialTreeBSCall(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    auto sqrt_deltat {std::sqrt(data.maturity / static_cast<commonType>(N))};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(std::exp(data.rate * data.maturity / N) - d) / (u -d)};
    
    commonType sum {};

    for (int j {}; j <= N; ++j)
    {
        sum += binomial(N,j) * std::pow(p, j) * std::pow(1-p, N-j) * std::max(data.spot * std::pow(u, j) * std::pow(d, N-j) - data.strike , static_cast<commonType>(0));
    }
    return sum * std::exp(-data.rate * data.maturity);
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto binomialTreeBSPut(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    auto sqrt_deltat {std::sqrt(data.maturity / static_cast<commonType>(N))};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(std::exp(data.rate * data.maturity / N) - d) / (u -d)};
    
    commonType sum {};

    for (int j {}; j <= N; ++j)
    {
        sum += binomial(N,j) * std::pow(p, j) * std::pow(1-p, N-j) * std::max(data.strike - data.spot * std::pow(u, j) * std::pow(d, N-j), static_cast<commonType>(0));
    }
    return sum * std::exp(-data.rate * data.maturity);
}


int main()
{
    OptionDataBS myOption {100, 100, 0.01, 0.3, 30};
    int iterations {1030};
    std::cout << "Call price: $" << binomialTreeBSCall(myOption, iterations) << '\n';
    std::cout << "Put price: $" << binomialTreeBSPut(myOption, iterations) << '\n';


    return 0;  
}

// after n = 1030 it breaks