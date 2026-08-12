#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <vector>
#include <numeric>
#include "../../Headers/jump_diffusion_struct.hpp"

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionLogEulerSimCall(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int M, int N)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    commonType deltaT {data.maturity / static_cast<commonType>(N)};

    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::poisson_distribution P {data.jump_intensity * deltaT};
    
    std::vector<commonType> X (M, std::log(data.spot));
    std::vector<commonType> payoffs(M);    

    for (int i{0}; i < M; ++i)
    {
        for (int j{0}; j < N; ++j)
        {
            commonType deltaN {static_cast<commonType>(P(mt))}; 
            commonType Zj {0};   
            if (deltaN == 0)
            {
                Zj = 0;
            }
            else
            {
                std::normal_distribution J {deltaN * data.mean_log_jump, std::sqrt(deltaN) * data.jump_size_volatility };
                Zj = J(mt);
            }
            
            commonType Z {W(mt)};
            
            X[i] += (data.rate - data.jump_intensity * kappa_j - static_cast<commonType>(0.5) * data.volatility * data.volatility) * deltaT +  data.volatility * std::sqrt(deltaT) * Z + Zj;
        }
    }

    std::transform(X.begin(), X.end(), payoffs.begin(), [&data](commonType x){return std::max(std::exp(x) - data.strike, static_cast<commonType>(0));});
    return std::exp(-data.rate * data.maturity) * std::accumulate(payoffs.begin(), payoffs.end(), static_cast<commonType>(0)) /static_cast<commonType>(payoffs.size());
}

template <typename S_t, typename K_t, typename r_t, typename sigma_t, typename T_t, typename lambda_t, typename mu_t, typename j_t>
auto jumpDiffusionLogEulerSimPut(OptionDataJD<S_t, K_t, r_t, sigma_t, T_t, lambda_t, mu_t, j_t>& data, int M, int N)
{
    using commonType = decltype(data.spot);

    commonType kappa_j {std::exp(data.mean_log_jump + static_cast<commonType>(0.5) * data.jump_size_volatility * data.jump_size_volatility) - static_cast<commonType>(1)};
    commonType deltaT {data.maturity / static_cast<commonType>(N)};

    std::random_device rd{};
    std::mt19937 mt{rd()};
    std::normal_distribution W {static_cast<commonType>(0), static_cast<commonType>(1)};
    std::poisson_distribution P {data.jump_intensity * deltaT};
    
    std::vector<commonType> X (M, std::log(data.spot));
    std::vector<commonType> payoffs(M);    

    for (int i{0}; i < M; ++i)
    {
        for (int j{0}; j < N; ++j)
        {
            commonType deltaN {static_cast<commonType>(P(mt))};
            commonType Zj {0}; 

            if (deltaN == 0)
            {
                Zj = 0;
            }
            else
            {
                std::normal_distribution J {deltaN * data.mean_log_jump, std::sqrt(deltaN) * data.jump_size_volatility };
                Zj = J(mt);
            }
            
            commonType Z {W(mt)};
            
            X[i] += (data.rate - data.jump_intensity * kappa_j - static_cast<commonType>(0.5) * data.volatility * data.volatility) * deltaT +  data.volatility * std::sqrt(deltaT) * Z + Zj;
        }
    }

    std::transform(X.begin(), X.end(), payoffs.begin(), [&data](commonType x){return std::max(data.strike - std::exp(x), static_cast<commonType>(0));});
    return std::exp(-data.rate * data.maturity) * std::accumulate(payoffs.begin(), payoffs.end(), static_cast<commonType>(0)) / static_cast<commonType>(payoffs.size());
}

int main()
{
    OptionDataJD values {100, 100, 0.05 ,0.20, 1, 0.75, -0.10, 0.25};

    std::cout << jumpDiffusionLogEulerSimCall(values, 100000, 100) << '\n';
    std::cout << jumpDiffusionLogEulerSimPut(values, 100000, 100) << '\n';
}