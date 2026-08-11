#include <iostream>
#include <cmath>
#include <vector>
#include <type_traits>
#include <algorithm>
#include "../../Headers/black76_structs.hpp"

template <typename F_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto binomialTreeB76Call(const OptionDataB76<F_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.future)
{
    using commonType = decltype(data.future);
    
    auto sqrt_deltat {std::sqrt(data.maturity / static_cast<commonType>(N))};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(1 - d) / (u -d)};

    std::vector<commonType> V_current (N + 1);
    for (int j {}; j < N + 1; ++j)
        {
            V_current[j] = std::max(data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>(j)) * std::pow(static_cast<commonType>(d), static_cast<commonType>(N-j)) - data.strike , static_cast<commonType>(0));
        }<

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

template <typename F_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto binomialTreeB76Put(const OptionDataB76<F_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.future)
{
    using commonType = decltype(data.future);
    
    auto sqrt_deltat {std::sqrt(data.maturity / static_cast<commonType>(N))};

    auto u {std::exp(data.volatility * sqrt_deltat)};
    auto d {std::exp(-data.volatility * sqrt_deltat)};
    auto p {(1 - d) / (u -d)};

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
            V_new[k] = p * V_current[k+1] + (1-p) * V_current[k];
        }
        std::swap(V_current, V_new);
        
    }   
    return V_current[0] * std::exp(-data.rate * data.maturity);
}

int main()
{
    OptionDataB76 myOption {80, 100, 0.01, 0.3, 30};
    int iterations {10000};
    std::cout << "Call price: $" << binomialTreeB76Call(myOption, 10000) << '\n';
    std::cout << "Put price: $" << binomialTreeB76Put(myOption, 10000) << '\n';


    return 0;  
}
