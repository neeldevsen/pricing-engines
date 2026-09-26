#include <iostream>
#include <cmath>
#include <vector>
#include <type_traits>
#include <algorithm>
#include "../../Headers/black_scholes_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto americanBinomialTreeBSCall(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    
    auto deltat {data.maturity / static_cast<commonType>(N)};
    auto sqrt_deltat {std::sqrt(deltat)};
    auto exp_deltat {std::exp(-data.rate * deltat)};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(std::exp(data.rate * deltat) - d) / (u -d)};
    commonType intrinsic {};
    commonType continuation {};

    std::vector<commonType> V_current (N + 1);
    for (int j {}; j < N + 1; ++j)
        {
            V_current[j] = std::max(data.spot * std::pow(static_cast<commonType>(u) , static_cast<commonType>(j)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-j)) - data.strike , static_cast<commonType>(0));
        }

    for (int i {}; i < N ; ++i)
    {
        std::vector<commonType> V_new (N-i);
        for (int k {}; k < N-i; ++k)
        {
            continuation = exp_deltat * (p * V_current[k+1] + (1-p) * V_current[k]);
            intrinsic = std::max(data.spot * std::pow(static_cast<commonType>(u) , static_cast<commonType>(k)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-i-k-1)) - data.strike , static_cast<commonType>(0));
            V_new[k] = std::max(continuation, intrinsic);
        }
        std::swap(V_current, V_new);
        
    }   
   
    return V_current[0];
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto americanBinomialTreeBSPut(const OptionDataBS<S_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.spot)
{
    using commonType = decltype(data.spot);
    
    auto deltat {data.maturity / static_cast<commonType>(N)};
    auto sqrt_deltat {std::sqrt(deltat)};
    auto exp_deltat {std::exp(-data.rate * deltat)};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(std::exp(data.rate * deltat) - d) / (u -d)};
    
    commonType intrinsic {};
    commonType continuation {};

    std::vector<commonType> V_current (N + 1);
    for (int j {}; j < N + 1; ++j)
        {
            V_current[j] = std::max(data.strike - data.spot * std::pow(static_cast<commonType>(u) , static_cast<commonType>(j)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-j)) , static_cast<commonType>(0));
        }

    for (int i {}; i < N ; ++i)
    {
        std::vector<commonType> V_new (N-i);
        for (int k {}; k < N-i; ++k)
        {
            continuation = exp_deltat * (p * V_current[k+1] + (1-p) * V_current[k]);
            intrinsic = std::max(data.strike - data.spot * std::pow(static_cast<commonType>(u) , static_cast<commonType>(k)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-i-k-1)), static_cast<commonType>(0));
            V_new[k] = std::max(continuation, intrinsic);
        }
        std::swap(V_current, V_new);
        
    }   
    return V_current[0];
}

int main()
{
    OptionDataBS myOption {100, 100, 0.01, 0.3, 30};
    int iterations {10000};
    std::cout << "Call price: $" << americanBinomialTreeBSCall(myOption, 10000) << '\n';
    std::cout << "Put price: $" << americanBinomialTreeBSPut(myOption, 10000) << '\n';

    return 0;  
}
