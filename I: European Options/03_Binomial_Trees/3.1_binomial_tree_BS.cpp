#include <iostream>
#include <cmath>
#include <vector>
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

    std::vector<commonType> V_current (N + 1);
    for (int j {}; j < N + 1; ++j)
        {
            V_current[j] = std::max(data.spot * std::pow(static_cast<double>(u) , static_cast<double>(j)) * std::pow(static_cast<double>(d), static_cast<double>(N-j)) - data.strike , static_cast<commonType>(0));
        }

    for (int i {}; i < N ; ++i)
    {
        std::vector<commonType> V_new (N-i);
        for (int k {}; k < N-i; ++k)
        {
            V_new[k] = p * V_current[k+1] + (1-p) * V_current[k];
        }
        std::swap(V_current, V_new);
        
    }   
   
    return V_current[0] * std::exp(-data.rate * data.maturity);
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto binomialTreeBSPut(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    
    auto sqrt_deltat {std::sqrt(data.maturity / static_cast<commonType>(N))};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(std::exp(data.rate * data.maturity / N) - d) / (u -d)};

    std::vector<commonType> V_current (N + 1);
    for (int j {}; j < N + 1; ++j)
        {
            V_current[j] = std::max(data.strike - data.spot * std::pow(static_cast<double>(u) , static_cast<double>(j)) * std::pow(static_cast<double>(d), static_cast<double>(N-j)) , static_cast<commonType>(0));
        }

    for (int i {}; i < N ; ++i)
    {
        std::vector<commonType> V_new (N-i);
        for (int k {}; k < N-i; ++k)
        {
            V_new[k] = p * V_current[k+1] + (1-p) * V_current[k];
        }
        std::swap(V_current, V_new);
        
    }   
    return V_current[0] * std::exp(-data.rate * data.maturity);
}

int main()
{
    OptionDataBS myOption {100, 100, 0.01, 0.3, 30};
    int iterations {10000};
    std::cout << "Call price: $" << binomialTreeBSCall(myOption, 10000) << '\n';
    std::cout << "Put price: $" << binomialTreeBSPut(myOption, 10000) << '\n';


    return 0;  
}
