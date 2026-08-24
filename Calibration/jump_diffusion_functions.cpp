#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../Headers/jump_diffusion_struct.hpp"
#include "../Headers/functions.hpp"
#include <pybind11/pybind11.h>

namespace py = pybind11;

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


//credits to ChatGPT for the code below 

PYBIND11_MODULE(jump_diffusion, m)
{
    m.doc() = "Merton Jump Diffusion series pricing engine";

    m.def(
        "call",
        [](double spot,
           double strike,
           double rate,
           double volatility,
           double maturity,
           double jump_intensity,
           double mean_log_jump,
           double jump_size_volatility,
           int N_max)
        {
            OptionDataJD<
                double,
                double,
                double,
                double,
                double,
                double,
                double,
                double
            > option{
                spot,
                strike,
                rate,
                volatility,
                maturity,
                jump_intensity,
                mean_log_jump,
                jump_size_volatility
            };

            return jumpDiffusionMertonFormulaCall(option, N_max);
        },
        py::arg("spot"),
        py::arg("strike"),
        py::arg("rate"),
        py::arg("volatility"),
        py::arg("maturity"),
        py::arg("jump_intensity"),
        py::arg("mean_log_jump"),
        py::arg("jump_size_volatility"),
        py::arg("N_max") = 50
    );

    m.def(
        "put",
        [](double spot,
           double strike,
           double rate,
           double volatility,
           double maturity,
           double jump_intensity,
           double mean_log_jump,
           double jump_size_volatility,
           int N_max)
        {
            OptionDataJD<
                double,
                double,
                double,
                double,
                double,
                double,
                double,
                double
            > option{
                spot,
                strike,
                rate,
                volatility,
                maturity,
                jump_intensity,
                mean_log_jump,
                jump_size_volatility
            };

            return jumpDiffusionMertonFormulaPut(option, N_max);
        },
        py::arg("spot"),
        py::arg("strike"),
        py::arg("rate"),
        py::arg("volatility"),
        py::arg("maturity"),
        py::arg("jump_intensity"),
        py::arg("mean_log_jump"),
        py::arg("jump_size_volatility"),
        py::arg("N_max") = 50
    );
}

