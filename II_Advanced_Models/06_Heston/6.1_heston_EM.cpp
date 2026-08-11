#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/heston_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename T_t, typename v_t, typename kappa_t, typename theta_t, typename xi_t, typename rho_t>
auto hestonNaiveEulerCall(OptionDataHeston<S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t>& data, int M, int N)
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W1 {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::normal_distribution W2 {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::vector<commonType> stock (M, data.spot);
    std::vector<commonType> payoff (M);
    commonType deltaT {data.maturity / static_cast<commonType>(N)};

    for (int i {0}; i < M; ++ i)
    {
        commonType variance {data.variance};
        for (int j {0}; j < N; ++j)
        {
            commonType Z1 {W1(mt)};
            commonType Z2 {W2(mt)};
            commonType Zv{data.rho * Z1 + std::sqrt(1- data.rho * data.rho) * Z2};
            stock[i] += data.rate * stock[i] * deltaT + stock[i] * std::sqrt(variance * deltaT) * Z1;
            variance += data.kappa * (data.theta - variance) * deltaT + data.xi * std::sqrt(variance * deltaT) * Zv;
        }
    }

    std::transform(stock.begin(), stock.end(), payoff.begin(), [&data](commonType S) {return std::max(S - data.strike, static_cast<commonType>(0));});

    return std::exp(-data.rate * data.maturity) * std::accumulate(payoff.begin(), payoff.end(), commonType{0}) / payoff.size(); 
}

template <typename S_t, typename K_t, typename r_t, typename T_t, typename v_t, typename kappa_t, typename theta_t, typename xi_t, typename rho_t>
auto hestonNaiveEulerPut(OptionDataHeston<S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t>& data, int M, int N)
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W1 {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::normal_distribution W2 {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::vector<commonType> stock (M, data.spot);
    std::vector<commonType> payoff (M);
    commonType deltaT {data.maturity / static_cast<commonType>(N)};

    for (int i {0}; i < M; ++ i)
    {
        commonType variance {data.variance};
        for (int j {0}; j < N; ++j)
        {
            commonType Z1 {W1(mt)};
            commonType Z2 {W2(mt)};
            commonType Zv{data.rho * Z1 + std::sqrt(1- data.rho * data.rho) * Z2};
            stock[i] += data.rate * stock[i] * deltaT + stock[i] * std::sqrt(variance * deltaT) * Z1;
            variance += data.kappa * (data.theta - variance) * deltaT + data.xi * std::sqrt(variance * deltaT) * Zv;
        }
    }

    std::transform(stock.begin(), stock.end(), payoff.begin(), [&data](commonType S) {return std::max(data.strike - S, static_cast<commonType>(0));});

    return std::exp(-data.rate * data.maturity) * std::accumulate(payoff.begin(), payoff.end(), commonType{0}) / payoff.size(); 
}
