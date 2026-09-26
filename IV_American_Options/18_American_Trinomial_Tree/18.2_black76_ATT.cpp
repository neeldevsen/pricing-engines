#include <iostream>
#include <cmath>
#include <vector>
#include <type_traits>
#include <algorithm>
#include "../../Headers/black76_structs.hpp"

template <typename F_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto americanTrinomialTreeB76Call(const OptionDataB76<F_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.future)
{
    using commonType = decltype(data.future);
    
    auto sqrt_2deltat {std::sqrt(2 * data.maturity / static_cast<commonType>(N))};
    auto exp_deltat {std::exp(-data.rate * data.maturity / N)};

    auto u {std::exp(data.volatility * sqrt_2deltat)};
    auto d {std::exp(-data.volatility * sqrt_2deltat)};
    auto p_u {std::pow((1 - std::exp(-data.volatility * sqrt_2deltat / static_cast<commonType>(2) )) / ((std::exp(data.volatility * sqrt_2deltat / static_cast<commonType>(2))) - std::exp(-data.volatility * sqrt_2deltat / static_cast<commonType>(2))), static_cast<commonType>(2))};
    auto p_d {std::pow((std::exp(data.volatility * sqrt_2deltat / static_cast<commonType>(2)) - 1) / ((std::exp(data.volatility * sqrt_2deltat / static_cast<commonType>(2))) - std::exp(-data.volatility * sqrt_2deltat / static_cast<commonType>(2))), static_cast<commonType>(2))};
    auto p_m {static_cast<commonType>(1) - p_u - p_d};

    commonType intrinsic {};
    commonType continuation {};

    std::vector<commonType> V_current (2 * N + 1);
    for (int j {}; j < 2 * N + 1; ++j)
        {
            V_current[j] = std::max(data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>(j- N)) - data.strike, static_cast<commonType>(0));
        }

    for (int i {}; i < N ; ++i)
    {
        std::vector<commonType> V_new (2 * N - 2 * i - 1);
        for (int k {}; k < 2* N-2*i -1; ++k)
        {
            continuation = exp_deltat *  (p_u * V_current[k+2] + p_m * V_current[k+1] + p_d * V_current[k]);
            intrinsic = std::max(data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>( k - N + i + 1)) - data.strike , static_cast<commonType>(0));
            V_new[k] = std::max(continuation, intrinsic);
        }
        std::swap(V_current, V_new);
    }   
   
    return V_current[0];
}

template <typename F_t, typename K_t, typename r_t, typename sigma_t, typename T_t>
auto americanTrinomialTreeB76Put(const OptionDataB76<F_t, K_t, r_t, sigma_t, T_t>& data, const int N) -> decltype(data.future)
{
    using commonType = decltype(data.future);
    
    auto sqrt_2deltat {std::sqrt(2 * data.maturity / static_cast<commonType>(N))};
    auto exp_deltat {std::exp(-data.rate * data.maturity / N)};

    auto u {std::exp(data.volatility * sqrt_2deltat)};
    auto d {std::exp(-data.volatility * sqrt_2deltat)};
    auto p_u {std::pow((1 - std::exp(-data.volatility * sqrt_2deltat / static_cast<commonType>(2) )) / ((std::exp(data.volatility * sqrt_2deltat / static_cast<commonType>(2))) - std::exp(-data.volatility * sqrt_2deltat / static_cast<commonType>(2))), static_cast<commonType>(2))};
    auto p_d {std::pow((std::exp(data.volatility * sqrt_2deltat / static_cast<commonType>(2)) - 1) / ((std::exp(data.volatility * sqrt_2deltat / static_cast<commonType>(2))) - std::exp(-data.volatility * sqrt_2deltat / static_cast<commonType>(2))), static_cast<commonType>(2))};
    auto p_m {static_cast<commonType>(1) - p_u - p_d};
    commonType intrinsic {};
    commonType continuation {};

    std::vector<commonType> V_current (2 * N + 1);
    for (int j {}; j < 2 * N + 1; ++j)
        {
            V_current[j] = std::max(data.strike - data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>(j- N)), static_cast<commonType>(0));
        }

    for (int i {}; i < N ; ++i)
    {
        std::vector<commonType> V_new (2 * N - 2 * i - 1);
        for (int k {}; k < 2* N-2*i -1; ++k)
        {
            continuation = exp_deltat *  (p_u * V_current[k+2] + p_m * V_current[k+1] + p_d * V_current[k]);
            intrinsic = std::max(data.strike - data.future * std::pow(static_cast<commonType>(u) , static_cast<commonType>( k - N + i + 1)), static_cast<commonType>(0));
            V_new[k] = std::max(continuation, intrinsic);
        }
        std::swap(V_current, V_new);
    }   
   
    return V_current[0];
}

int main()
{
    OptionDataB76 myOption {80, 100, 0.01, 0.3, 30};
    int iterations {10000};
    std::cout << "Call price: $" << americanTrinomialTreeB76Call(myOption, iterations) << '\n';
    std::cout << "Put price: $" << americanTrinomialTreeB76Put(myOption, iterations) << '\n';


    return 0;  
}
