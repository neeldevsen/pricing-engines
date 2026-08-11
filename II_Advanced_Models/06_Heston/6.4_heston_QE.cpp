#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/heston_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename T_t, typename v_t, typename kappa_t, typename theta_t, typename xi_t, typename rho_t>
auto hestonQECall(OptionDataHeston<S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t>& data, int M, int N, double threshold, double gamma1, double gamma2)
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::uniform_real_distribution U {static_cast<commonType>(0), static_cast<commonType>(1)};

    std::vector<commonType> X (M, std::log(data.spot));
    std::vector<commonType> payoff (M);

    commonType deltaT {data.maturity / static_cast<commonType>(N)};
    commonType gamma1_deltaT {static_cast<commonType>(gamma1) * deltaT};
    commonType gamma2_deltaT {static_cast<commonType>(gamma2) * deltaT};
    commonType exp_factor {std::exp(-data.kappa * deltaT)};

    commonType K0 {-data.rho * data.kappa * data.theta * deltaT / data.xi};
    commonType K1 {gamma1_deltaT * (data.kappa * data.rho / data.xi - static_cast<commonType>(0.5))- data.rho / data.xi};
    commonType K2 {gamma2_deltaT * (data.kappa * data.rho / data.xi - static_cast<commonType>(0.5))+ data.rho / data.xi};
    commonType K3 {gamma1_deltaT * (static_cast<commonType>(1) - data.rho * data.rho)};
    commonType K4 {gamma2_deltaT * (static_cast<commonType>(1) - data.rho * data.rho)};


    for (int i {0}; i < M; ++i)
    {
        commonType variance_n {data.variance};
        commonType variance_n_p1 {data.variance};
        for (int j {0}; j < N; ++j)
        {
            
            commonType s_squared {variance_n * data.xi * data.xi * exp_factor * (static_cast<commonType>(1) - exp_factor) / data.kappa + data.theta * std::pow(data.xi * (1-exp_factor), static_cast<commonType>(2.0)) / (2 * data.kappa) };
            commonType m {data.theta + (variance_n - data.theta) * exp_factor};
            commonType phi {s_squared / std::pow(m, static_cast<commonType>(2))};
            commonType Z1 {W(mt)};
            commonType Z2 {W(mt)};

            if (phi < static_cast<commonType>(threshold))
            {
                commonType b_squared {static_cast<commonType>(2) / phi - static_cast<commonType>(1) + std::sqrt((static_cast<commonType>(2) / phi) * (static_cast<commonType>(2) / phi - static_cast<commonType>(1) ))};
                commonType a {m / (static_cast<commonType>(1) + b_squared)};
                variance_n_p1 = a * std::pow((std::sqrt(b_squared) + Z1), static_cast<commonType>(2));
            }
            else
            {
                commonType p {(phi - static_cast<commonType>(1)) / (phi + static_cast<commonType>(1))};
                commonType beta {(static_cast<commonType>(1) - p) / m};
                commonType u {U(mt)};
                if (u <= p)
                {
                    variance_n_p1 = static_cast<commonType>(0);
                }
                else
                {
                    variance_n_p1 = -std::log((static_cast<commonType>(1) - p) / (static_cast<commonType>(1) - u)) / beta; 
                }
            }
            
            X[i] += data.rate * deltaT + K0 + K1 * variance_n + K2 * variance_n_p1 + std::sqrt(K3 * variance_n + K4 * variance_n_p1) * Z2;

            variance_n = variance_n_p1;
        }
    }

    std::transform(X.begin(), X.end(), payoff.begin(), [&data](commonType x) {return std::max(std::exp(x) - data.strike, static_cast<commonType>(0));});

    return std::exp(-data.rate * data.maturity) * std::accumulate(payoff.begin(), payoff.end(), commonType{0}) / payoff.size(); 
}

#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/heston_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename T_t, typename v_t, typename kappa_t, typename theta_t, typename xi_t, typename rho_t>
auto hestonQEPut(OptionDataHeston<S_t, K_t, r_t, T_t, v_t, kappa_t, theta_t, xi_t, rho_t>& data, int M, int N, double threshold, double gamma1, double gamma2)
{
    using commonType = decltype(data.spot);
    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::uniform_real_distribution U {static_cast<commonType>(0), static_cast<commonType>(1)};

    std::vector<commonType> X (M, std::log(data.spot));
    std::vector<commonType> payoff (M);

    commonType deltaT {data.maturity / static_cast<commonType>(N)};
    commonType gamma1_deltaT {static_cast<commonType>(gamma1) * deltaT};
    commonType gamma2_deltaT {static_cast<commonType>(gamma2) * deltaT};
    commonType exp_factor {std::exp(-data.kappa * deltaT)};

    commonType K0 {-data.rho * data.kappa * data.theta * deltaT / data.xi};
    commonType K1 {gamma1_deltaT * (data.kappa * data.rho / data.xi - static_cast<commonType>(0.5))- data.rho / data.xi};
    commonType K2 {gamma2_deltaT * (data.kappa * data.rho / data.xi - static_cast<commonType>(0.5))+ data.rho / data.xi};
    commonType K3 {gamma1_deltaT * (static_cast<commonType>(1) - data.rho * data.rho)};
    commonType K4 {gamma2_deltaT * (static_cast<commonType>(1) - data.rho * data.rho)};


    for (int i {0}; i < M; ++i)
    {
        commonType variance_n {data.variance};
        commonType variance_n_p1 {data.variance};
        for (int j {0}; j < N; ++j)
        {
            
            commonType s_squared {variance_n * data.xi * data.xi * exp_factor * (static_cast<commonType>(1) - exp_factor) / data.kappa + data.theta * std::pow(data.xi * (1-exp_factor), static_cast<commonType>(2.0)) / (2 * data.kappa) };
            commonType m {data.theta + (variance_n - data.theta) * exp_factor};
            commonType phi {s_squared / std::pow(m, static_cast<commonType>(2))};
            commonType Z1 {W(mt)};
            commonType Z2 {W(mt)};

            if (phi < static_cast<commonType>(threshold))
            {
                commonType b_squared {static_cast<commonType>(2) / phi - static_cast<commonType>(1) + std::sqrt((static_cast<commonType>(2) / phi) * (static_cast<commonType>(2) / phi - static_cast<commonType>(1) ))};
                commonType a {m / (static_cast<commonType>(1) + b_squared)};
                variance_n_p1 = a * std::pow((std::sqrt(b_squared) + Z1), static_cast<commonType>(2));
            }
            else
            {
                commonType p {(phi - static_cast<commonType>(1)) / (phi + static_cast<commonType>(1))};
                commonType beta {(static_cast<commonType>(1) - p) / m};
                commonType u {U(mt)};
                if (u <= p)
                {
                    variance_n_p1 = static_cast<commonType>(0);
                }
                else
                {
                    variance_n_p1 = -std::log((static_cast<commonType>(1) - p) / (static_cast<commonType>(1) - u)) / beta; 
                }
            }
            
            X[i] += data.rate * deltaT + K0 + K1 * variance_n + K2 * variance_n_p1 + std::sqrt(K3 * variance_n + K4 * variance_n_p1) * Z2;

            variance_n = variance_n_p1;
        }
    }

    std::transform(X.begin(), X.end(), payoff.begin(), [&data](commonType x) {return std::max(data.strike - std::exp(x), static_cast<commonType>(0));});

    return std::exp(-data.rate * data.maturity) * std::accumulate(payoff.begin(), payoff.end(), commonType{0}) / payoff.size(); 
}

