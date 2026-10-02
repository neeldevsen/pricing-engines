#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include "barrier_struct.hpp"
#include <numeric>

template <typename S_t, typename K_t, typename B_t, typename r_t, typename sigma_t, typename T_t>
auto monteCarloBarrierUp(const OptionDataBarrier<S_t, K_t, B_t, r_t, sigma_t, T_t>& data, int simulations, int steps) -> OptionValueBarrierUp<decltype(data.spot)>
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    commonType payoffs_call_in {};
    commonType payoffs_call_out {};
    commonType payoffs_put_in {};
    commonType payoffs_put_out {};
    commonType deltat {data.maturity / static_cast<commonType>(steps)};
    OptionValueBarrierUp<commonType> values {};

    for (int i {}; i < simulations; ++i)
    {
        bool hit_barrier {false};
        commonType S_final {data.spot};
        for (int j {}; j < steps; ++j)
        {
            S_final *= std::exp((data.rate - 0.5 * data.volatility * data.volatility) * deltat + (data.volatility * std::sqrt(deltat) * W(mt)));
            if (S_final >= data.barrier && !hit_barrier)
            {
                hit_barrier = true;
            }    
        }

        if (hit_barrier)
            {
                payoffs_call_in += std::max(S_final - data.strike, static_cast<commonType>(0));
                payoffs_put_in += std::max(data.strike - S_final, static_cast<commonType>(0));
            }
        else
            {
                payoffs_call_out += std::max(S_final - data.strike, static_cast<commonType>(0));
                payoffs_put_out += std::max(data.strike - S_final, static_cast<commonType>(0));
            }
    }
    
    values.callUpIn = std::exp(-data.rate * data.maturity) * payoffs_call_in / static_cast<commonType>(simulations);
    values.callUpOut = std::exp(-data.rate * data.maturity) * payoffs_call_out / static_cast<commonType>(simulations);
    values.putUpIn = std::exp(-data.rate * data.maturity) * payoffs_put_in / static_cast<commonType>(simulations);
    values.putUpOut = std::exp(-data.rate * data.maturity)  * payoffs_put_out / static_cast<commonType>(simulations);

    return values;
}

template <typename S_t, typename K_t, typename B_t, typename r_t, typename sigma_t, typename T_t>
auto monteCarloBarrierDown(const OptionDataBarrier<S_t, K_t, B_t, r_t, sigma_t, T_t>& data, int simulations, int steps) -> OptionValueBarrierDown<decltype(data.spot)>
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    commonType payoffs_call_in {};
    commonType payoffs_call_out {};
    commonType payoffs_put_in {};
    commonType payoffs_put_out {};
    commonType deltat {data.maturity / static_cast<commonType>(steps)};
    OptionValueBarrierDown<commonType> values {};

    for (int i {}; i < simulations; ++i)
    {
        bool hit_barrier {false};
        commonType S_final {data.spot};
        for (int j {}; j < steps; ++j)
        {
            S_final *= std::exp((data.rate - 0.5 * data.volatility * data.volatility) * deltat + (data.volatility * std::sqrt(deltat) * W(mt)));
            if (S_final <= data.barrier && !hit_barrier)
            {
                hit_barrier = true;
            }    
        }
        
        if (hit_barrier)
            {
                payoffs_call_in += std::max(S_final - data.strike, static_cast<commonType>(0));
                payoffs_put_in += std::max(data.strike - S_final, static_cast<commonType>(0));
            }
        else
            {
                payoffs_call_out += std::max(S_final - data.strike, static_cast<commonType>(0));
                payoffs_put_out += std::max(data.strike - S_final, static_cast<commonType>(0));
            }
    }
    
    values.callDownIn = std::exp(-data.rate * data.maturity) * payoffs_call_in / static_cast<commonType>(simulations);
    values.callDownOut = std::exp(-data.rate * data.maturity) * payoffs_call_out / static_cast<commonType>(simulations);
    values.putDownIn = std::exp(-data.rate * data.maturity) * payoffs_put_in / static_cast<commonType>(simulations);
    values.putDownOut = std::exp(-data.rate * data.maturity) * payoffs_put_out / static_cast<commonType>(simulations);

    return values;
}


int main()
{
    int simulations {100000};
    int steps {1000};
    // Up options
    OptionDataBarrier myOptionUp {100, 100, 120, 0.05, 0.2, 1};
    OptionValueBarrierUp<double> valuesUp {monteCarloBarrierUp(myOptionUp, simulations, steps)};

    std::cout << "Up-Out Call Option price: $" << valuesUp.callUpOut << '\n';
    std::cout << "Up-In Call Option price: $" << valuesUp.callUpIn << '\n';
    std::cout << "Up-Out Put Option price: $" << valuesUp.putUpOut << '\n';
    std::cout << "Up-In Put Option price: $" << valuesUp.putUpIn << '\n';

    // Down options
    OptionDataBarrier myOptionDown {100, 100, 80, 0.05, 0.2, 1};
    OptionValueBarrierDown<double> valuesDown {monteCarloBarrierDown(myOptionDown, simulations, steps)};

    
    std::cout << "Down-Out Call Option price: $" << valuesDown.callDownOut << '\n';
    std::cout << "Down-In Call Option price: $" << valuesDown.callDownIn << '\n';
    std::cout << "Down-Out Put Option price: $" << valuesDown.putDownOut << '\n';
    std::cout << "Down-In Put Option price: $" << valuesDown.putDownIn << '\n';

    return 0;
}
