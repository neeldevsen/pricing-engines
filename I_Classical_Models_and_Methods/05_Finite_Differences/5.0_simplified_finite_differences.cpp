#include <iostream>
#include <cmath>
#include <vector>
#include <type_traits>
#include <algorithm>
#include "../../Headers/black_scholes_struct.hpp"



std::vector<double> generateAlphas(double sigma, double S_max, int M, double r)
{
    double deltaS {S_max / M};
    std::vector<double> alphas(M-1);
    for (int i{0}; i < M-1; ++i)
    {
        double S_i {(i+1)*deltaS};
        alphas[i] = 0.5 * (std::pow(sigma*S_i, 2.0)) / (std::pow(deltaS, 2.0)) - (r*S_i) / (2*deltaS);
    }
    return alphas;
}

std::vector<double> generateBetas(double sigma, double S_max, int M, double r)
{
    double deltaS {S_max / M};
    std::vector<double> betas(M-1);
    for (int i{0}; i < M-1; ++i)
    {
        double S_i {(i+1)*deltaS};
        betas[i] = -(std::pow(sigma*S_i, 2.0)) / (std::pow(deltaS, 2.0)) - r;
    }
    return betas;
}

std::vector<double> generateGammas(double sigma, double S_max, int M, double r)
{
    double deltaS {S_max / M};
    std::vector<double> gammas(M-1);
    for (int i{0}; i < M-1; ++i)
    {
        double S_i {(i+1)*deltaS};
        gammas[i] = 0.5 * (std::pow(sigma*S_i, 2.0)) / (std::pow(deltaS, 2.0)) + (r*S_i) / (2*deltaS);

    }
    return gammas;
}

std::vector<double> generateALower(std::vector<double>& alphas, double deltaT, int M) // flag
{
    std::vector<double> A_lower (M-2);
    for (int i {}; i < M-2; ++i)
    {
        A_lower[i] = - 0.5 * deltaT * alphas[i+1];
    }
    return A_lower;
}

std::vector<double> generateAMain(std::vector<double>& betas, double deltaT, int M)
{
    std::vector<double> A_main (M-1);
    for (int i {}; i < M-1; ++i)
    {
        A_main[i] = 1 - 0.5 * deltaT * betas[i];
    }
    return A_main;
}

std::vector<double> generateAUpper(std::vector<double>& gammas, double deltaT, int M) // flag
{
    std::vector<double> A_upper (M-2);
    for (int i {}; i < M-2; ++i)
    {
        A_upper[i] = - 0.5 * deltaT * gammas[i];
    }
    return A_upper;
}


std::vector<double> generateRHS(std::vector<double>& alphas, std::vector<double>& betas, std::vector<double>& gammas, std::vector<double> V, double deltaT, int M, double K, int n, double r, double Smax)
{
    V[M] = Smax - K * std::exp(-r*n*deltaT); // changed for put
    V[0] = 0;
    std::vector<double> RHS (M-1);
    for (int i{}; i < M -1; ++i)
    {
        RHS[i] = 0.5*deltaT*alphas[i]*V[i] + (1 + 0.5*deltaT*betas[i])*V[i+1] + 0.5*deltaT * gammas[i]*V[i+2];
    }
    RHS[M - 2] -= (-0.5 * deltaT * gammas[M - 2]) * (Smax - K * std::exp(-r * (n + 1) * deltaT));
    return RHS;
}

std::vector<double> thomasAlgorithm(std::vector<double>& A_lower, std::vector<double>& A_main, std::vector<double>& A_upper, std::vector<double> RHS)
{
    int M {static_cast<int>(A_main.size()) + 1};
    std::vector<double> rho (M-2);
    std::vector<double> eta (M-1);
    std::vector<double> mu (M-1);

    for (int i{0}; i < M - 1; ++i)
    {
        mu[i] = RHS[i] / A_main[i];
        if (i == M - 2) {break;}
        eta[i+1] = A_lower[i] / A_main[i+1];
        rho[i] = A_upper[i] / A_main[i]; 
    }

    std::vector<double> rho_p (M-2);
    std::vector<double> mu_p (M-1);

    rho_p[0] = rho[0];
    mu_p[0] = mu[0];

    for (int i {1}; i < M - 1; ++i)
    {
        mu_p[i] = (mu[i] - eta[i] * mu_p[i-1]) / (1.0 - eta[i] * rho_p[i-1]);
        if (i == M-2) {break;}
        rho_p[i] = rho[i] / (1.0 - eta[i] * rho_p[i-1]);
    }
    
    std::vector<double> Vn (M);
    Vn[M-1] = mu_p[M-2];
    for (int i{M-2}; i >= 1; --i)
    {
        Vn[i] = mu_p[i-1] - rho_p[i-1] * Vn[i+1];
    }
    return Vn;
}

std::vector<double> stepCrankNicolson(std::vector<double>& A_lower, std::vector<double>& A_main, std::vector<double>& A_upper, std::vector<double>& alphas, std::vector<double>& betas, std::vector<double>& gammas, std::vector<double> V, double deltaT, int M, int n, double K, double r, double Smax)
{
    std::vector<double> RHS{generateRHS(alphas, betas, gammas, V, deltaT, M, K, n, r, Smax)};
    std::vector<double> V_next {thomasAlgorithm(A_lower, A_main, A_upper, RHS)};
    V_next[0] = 0;
    V_next.push_back(Smax - K * std::exp(-r*deltaT*(n+1)));
    return V_next;  
}

double solveEuroCallCrankNicolson(double St, double K, double r, double sigma, double T, int M, int N, double multiplier=4.0)
{
    double Smax {multiplier * std::max(St,K)};
    double deltaS {Smax / M};
    double deltaT {T / N};
    std::vector<double> S (M + 1);
    for (int i{}; i < M + 1; ++i)
    {
        S[i] = i * deltaS;
    }

    std::vector<double> V (M + 1);
    for (int i{}; i < M + 1; ++i)
    {
        V[i] = std::max(S[i] - K, static_cast<double>(0));
    }
    std::vector<double> alphas = generateAlphas(sigma, Smax, M, r);
    std::vector<double> betas = generateBetas(sigma, Smax, M, r);
    std::vector<double> gammas = generateGammas(sigma, Smax, M, r);

    std::vector<double> A_lower {generateALower(alphas, deltaT, M)};
    std::vector<double> A_main {generateAMain(betas, deltaT, M)};
    std::vector<double> A_upper {generateAUpper(gammas, deltaT, M)};

    for (int n{}; n < N; ++n)
    {
        std::vector<double> V_new {stepCrankNicolson(A_lower, A_main, A_upper, alphas, betas, gammas, V, deltaT, M, n, K, r, Smax)};   
        std::swap(V_new, V);
    }
    int j {static_cast<int>(St / deltaS)};
    j = std::max(static_cast<int>(0), static_cast<int>(std::min(M-1,static_cast<int>(j))));
    double price {V[j] + (St - S[j]) * (V[j+1] - V[j]) / (S[j+1] - S[j])};
    return price;
}

int main()
{
    double St {100};
    double Smult {40};
    double K {100};
    double rate {0.01};
    double volatility {0.3};
    double maturity {30};
    int M {5000};
    int N {5000};
    std::cout << "Call price: " << solveEuroCallCrankNicolson(St, K, rate, volatility, maturity, M, N, Smult);

    return 0;
}

