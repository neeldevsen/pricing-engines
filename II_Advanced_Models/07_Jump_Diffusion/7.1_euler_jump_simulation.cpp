#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/jump_diffusion_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionEulerSimCall(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int M, int N)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    commonType deltaT {data.maturity / static_cast<commonType>(N)};

    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::poisson_distribution P {data.jump_intensity * deltaT};
    
    std::vector<commonType> stocks (M, data.spot);
    std::vector<commonType> payoffs(M);    

    for (int i{0}; i < M; ++i)
    {
        for (int j{0}; j < N; ++j)
        {
            commonType deltaN {P(mt)};
            commonType Zn {W(mt)};
            commonType Zj {W(mt)};
            commonType J{std::exp(data.mean_log_jump + data.jump_size_volatility * Zj)};
            stocks[i] += (data.rate - data.jump_intensity * kappa_j) * stocks[i] * deltaT + data.volatility * stocks[i] * std::sqrt(deltaT) * Zn + stocks[i] * (J - static_cast<commonType>(1)) * deltaN;
        }
    }

    std::transform(stocks.begin(), stocks.end(), payoffs.begin(), [&data](commonType S){return std::max(S - data.strike, static_cast<commonType>(0));});
    return std::exp(-data.rate * data.maturity) * std::accumulate(payoffs.begin(), payoffs.end(), static_cast<commonType>(0)) /static_cast<commonType>(payoffs.size());
}

#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/jump_diffusion_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionEulerSimPut(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int M, int N)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    commonType deltaT {data.maturity / static_cast<commonType>(N)};

    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::poisson_distribution P {data.jump_intensity * deltaT}; 

    std::vector<commonType> stocks (M, data.spot);
    std::vector<commonType> payoffs(M);    

    for (int i{0}; i < M; ++i)
    {
        for (int j{0}; j < N; ++j)
        {
            commonType deltaN {P(mt)};
            commonType Zn {W(mt)};
            commonType Zj {W(mt)};
            commonType J{std::exp(data.mean_log_jump + data.jump_size_volatility * Zj)};
            stocks[i] += (data.rate - data.jump_intensity * kappa_j) * stocks[i] * deltaT + data.volatility * stocks[i] * std::sqrt(deltaT) * Zn + stocks[i] * (J - static_cast<commonType>(1)) * deltaN;
        }
    }

    std::transform(stocks.begin(), stocks.end(), payoffs.begin(), [&data](commonType S){return std::max(data.strike - S, static_cast<commonType>(0));});
    return std::exp(-data.rate * data.maturity) * std::accumulate(payoffs.begin(), payoffs.end(), static_cast<commonType>(0)) /static_cast<commonType>(payoffs.size());
}
