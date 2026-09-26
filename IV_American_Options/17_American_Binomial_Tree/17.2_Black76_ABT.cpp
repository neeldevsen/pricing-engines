#include <iostream>
#include <cmath>
#include <vector>
#include <type_traits>
#include <algorithm>
#include "../../Headers/black76_structs.hpp"

template <typename F_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto americanBinomialTreeB76Call(const OptionDataB76<F_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.future)
{
    using commonType = decltype(data.future);
    
    auto deltat {data.maturity / static_cast<commonType>(N)};
    auto sqrt_deltat {std::sqrt(deltat)};
    auto exp_deltat {std::exp(-data.rate * deltat)};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(1 - d) / (u -d)};

    commonType intrinsic {};
    commonType continuation {};

    std::vector<commonType> V_current (N + 1);
    for (int j {}; j < N + 1; ++j)
        {
            V_current[j] = std::max(data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>(j)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-j)) - data.strike , static_cast<commonType>(0));
        }

    for (int i {}; i < N ; ++i)
    {
        std::vector<commonType> V_new (N-i);
        for (int k {}; k < N-i; ++k)
        {
            continuation = exp_deltat * (p * V_current[k+1] + (1-p) * V_current[k]);
            intrinsic = std::max(data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>(k)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-i-k-1)) - data.strike , static_cast<commonType>(0));
            V_new[k] = std::max(continuation, intrinsic);
        }
        std::swap(V_current, V_new);
        
    }   
   
    return V_current[0];
}

template <typename F_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto americanBinomialTreeB76Put(const OptionDataB76<F_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.future)
{
    using commonType = decltype(data.future);
    
    auto deltat {data.maturity / static_cast<commonType>(N)};
    auto sqrt_deltat {std::sqrt(deltat)};
    auto exp_deltat {std::exp(-data.rate * deltat)};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(1 - d) / (u -d)};
    
    commonType intrinsic {};
    commonType continuation {};

    std::vector<commonType> V_current (N + 1);
    for (int j {}; j < N + 1; ++j)
        {
            V_current[j] = std::max(data.strike - data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>(j)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-j)) , static_cast<commonType>(0));
        }

    for (int i {}; i < N ; ++i)
    {
        std::vector<commonType> V_new (N-i);
        for (int k {}; k < N-i; ++k)
        {
            continuation = exp_deltat * (p * V_current[k+1] + (1-p) * V_current[k]);
            intrinsic = std::max(data.strike - data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>(k)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-i-k-1)), static_cast<commonType>(0));
            V_new[k] = std::max(continuation, intrinsic);
        }
        std::swap(V_current, V_new);
        
    }   
    return V_current[0];
}

int main()
{
    OptionDataB76 myOption {100, 100, 0.01, 0.3, 30};
    std::cout << "Call price: $" << americanBinomialTreeB76Call(myOption, 10000) << '\n';
    std::cout << "Put price: $" << americanBinomialTreeB76Put(myOption, 10000) << '\n';

    return 0;  
}
