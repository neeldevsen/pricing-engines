#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/jump_diffusion_struct.hpp"
#include "../../Headers/functions.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionMertonFormulaCall(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int N_max)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    commonType K_0 {std::log(data.spot) + (data.rate - data.jump_intensity * kappa_j - static_cast<commonType>(0.5) * data.volatility * data.volatility) * data.maturity};
    commonType K_1 {data.volatility * data.volatility * data.maturity};
    commonType K_2 {data.jump_size_volatility * data.jump_size_volatility};



    std::vector<commonType> probabilities (N_max + 1);
    probabilities[0] = std::exp(-data.jump_intensity * data.maturity);

    for (int n{1}; n <= N_max; ++n)
    {
        probabilities[n] = probabilities[n-1] * data.jump_intensity * data.maturity / n;
    }

    commonType total_payoff {};

    for (int n{0}; n <= N_max; ++n)
    {
        commonType m {K_0 + data.mean_log_jump * n};
        commonType v {K_1 + K_2 * n};
        commonType d_2 {(m - std::log(data.strike)) / (std::sqrt(v))};
        commonType d_1 {d_2 + std::sqrt(v)};
        commonType payoff {std::exp(m + static_cast<commonType>(0.5) * v) * normalCDF(d_1) - data.strike * normalCDF(d_2)};

        total_payoff += payoff * probabilities[n];

    }

    return std::exp(-data.rate * data.maturity) * total_payoff;
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionMertonFormulaPut(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int N_max)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    commonType K_0 {std::log(data.spot) + (data.rate - data.jump_intensity * kappa_j - static_cast<commonType>(0.5) * data.volatility * data.volatility) * data.maturity};
    commonType K_1 {data.volatility * data.volatility * data.maturity};
    commonType K_2 {data.jump_size_volatility * data.jump_size_volatility};



    std::vector<commonType> probabilities (N_max + 1);
    probabilities[0] = std::exp(-data.jump_intensity * data.maturity);

    for (int n{1}; n <= N_max; ++n)
    {
        probabilities[n] = probabilities[n-1] * data.jump_intensity * data.maturity / n;
    }

    commonType total_payoff {};

    for (int n{0}; n <= N_max; ++n)
    {
        commonType m {K_0 + data.mean_log_jump * n};
        commonType v {K_1 + K_2 * n};
        commonType d_2 {(m - std::log(data.strike)) / (std::sqrt(v))};
        commonType d_1 {d_2 + std::sqrt(v)};
        commonType payoff {data.strike * normalCDF(-d_2) - std::exp(m + static_cast<commonType>(0.5) * v) * normalCDF(-d_1)};

        total_payoff += payoff * probabilities[n];

    }

    return std::exp(-data.rate * data.maturity) * total_payoff;
}

int main()
{
    OptionDataJD values {100, 100, 0.05 ,0.20, 1, 0.75, -0.10, 0.25};

    std::cout << jumpDiffusionMertonFormulaCall(values, 10000) << '\n';
    std::cout << jumpDiffusionMertonFormulaPut(values, 10000) << '\n';

}