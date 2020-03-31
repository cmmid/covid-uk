// randomizer.h
// (C) 2013-2020 Nicholas G Davies

#ifndef RANDOMIZER_H
#define RANDOMIZER_H

#include <vector>
#include <limits>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>
#include <stdexcept>
#include <string>

class Randomizer
{
public:
    Randomizer(unsigned long int seed = 0);
    ~Randomizer();

    void Reset();

    double Uniform(double min = 0.0, double max = 1.0);
    double RoundedUniform(double min = 0.0, double max = 1.0, double shoulder = 0.01);
    double Normal(double mean = 0.0, double sd = 1.0);
    double Normal(double mean, double sd, double clamp);
    double Cauchy(double x0 = 0.0, double gamma = 1.0);
    double LogNormal(double zeta = 0.0, double sd = 1.0);
    double Exponential(double rate = 1.0);
    double Gamma(double shape, double scale);
    double Beta(double alpha, double beta);
    unsigned int Discrete(unsigned int size);
    int Discrete(int min, int max);
    // int Discrete(std::vector<unsigned int>& cumulative_weights);
    // int Discrete(std::vector<double>& cumulative_weights);
    void Multinomial(unsigned int N, std::vector<double>& p, std::vector<unsigned int>& n);
    bool Bernoulli(double p);
    unsigned int Binomial(unsigned int n, double p);
    unsigned int BetaBinomial(unsigned int n, double p, double a_plus_b);
    int Poisson(double mean);
    int Geometric(double p);
    int Round(double x);

    template <typename T>
    void Shuffle(std::vector<T>& vec);

    inline gsl_rng* GSL_RNG() { return r; }

private:
    unsigned long int seed;
    gsl_rng* r;
};

Randomizer::Randomizer(unsigned long int s)
 : seed(s), r(gsl_rng_alloc(gsl_rng_mt19937))
{
    Reset();
}

Randomizer::~Randomizer()
{
    gsl_rng_free(r);
}

void Randomizer::Reset()
{
    gsl_rng_set(r, seed);
}

double Randomizer::Uniform(double min, double max)
{
    return min + gsl_rng_uniform(r) * (max - min);
}

double Randomizer::RoundedUniform(double min, double max, double shoulder)
{
    if (min >= max)
        return min;
    double z = Uniform();
    double sd = shoulder * (max - min) / ((1 - shoulder) * 2.50662827463);
    if (z < shoulder / 2)
        return min - abs(Normal(0, sd));
    else if (z < shoulder)
        return max + abs(Normal(0, sd));
    else
        return Uniform(min, max);
}

double Randomizer::Normal(double mean, double sd)
{
    return mean + gsl_ran_gaussian_ziggurat(r, sd);
}

double Randomizer::Normal(double mean, double sd, double clamp)
{
    double n;
    do n = Normal(mean, sd); while (std::fabs(n - mean) > clamp);
    return n;
}

double Randomizer::LogNormal(double zeta, double sd)
{
    return gsl_ran_lognormal(r, zeta, sd);
}

double Randomizer::Cauchy(double x0, double gamma)
{
    return x0 + gsl_ran_cauchy(r, gamma);
}

double Randomizer::Exponential(double rate)
{
    return gsl_ran_exponential(r, 1. / rate);
}

double Randomizer::Gamma(double shape, double scale)
{
    return gsl_ran_gamma(r, shape, scale);
}

double Randomizer::Beta(double alpha, double beta)
{
    return gsl_ran_beta(r, alpha, beta);
}

unsigned int Randomizer::Discrete(unsigned int size)
{
    // TODO improve...
    if (size > gsl_rng_max(r))
        throw std::runtime_error("Generator cannot produce integers larger than " + std::to_string(gsl_rng_max(r)));
    return gsl_rng_get(r) % size;
}

int Randomizer::Discrete(int min, int max)
{
    return min + gsl_rng_uniform_int(r, max - min + 1);
}
// ///
// int Randomizer::Discrete(std::vector<unsigned int>& cumulative_weights)
// {

//     if (cumulative_weights.back() == 0)
//         return Discrete(cumulative_weights.size());
//     return std::distance(cumulative_weights.begin(),
//                 std::lower_bound(cumulative_weights.begin(),
//                     cumulative_weights.end(),
//                     Discrete(cumulative_weights.back())));
// }
// ///
// int Randomizer::Discrete(std::vector<double>& cumulative_weights)
// {
//     if (cumulative_weights.back() == 0)
//         return Discrete(cumulative_weights.size());
//     return std::distance(cumulative_weights.begin(),
//                 std::lower_bound(cumulative_weights.begin(),
//                     cumulative_weights.end(),
//                     Uniform(0.0, cumulative_weights.back())));
// }

void Randomizer::Multinomial(unsigned int N, std::vector<double>& p, std::vector<unsigned int>& n)
{
    gsl_ran_multinomial(r, p.size(), N, &p[0], &n[0]);
}

bool Randomizer::Bernoulli(double p)
{
    if (p <= 0) return false;
    if (p >= 1) return true;
    return gsl_rng_uniform(r) < p;
}

unsigned int Randomizer::Binomial(unsigned int n, double p)
{
    if (p <= 0) return 0;
    return gsl_ran_binomial(r, p, n);
}

unsigned int Randomizer::BetaBinomial(unsigned int n, double p, double a_plus_b)
{
    if (a_plus_b > 0)
        p = gsl_ran_beta(r, a_plus_b * p, a_plus_b * (1 - p));
    return gsl_ran_binomial(r, p, n);
}

int Randomizer::Poisson(double mean)
{
    if (mean <= 0) return 0;
    return gsl_ran_poisson(r, mean);
}

int Randomizer::Geometric(double p)
{
    if (p <= 0) return 0;
    return gsl_ran_geometric(r, p);
}

int Randomizer::Round(double x)
{
    int sign = x < 0 ? -1 : 1;
    double intpart, fracpart;
    fracpart = std::modf(std::fabs(x), &intpart);
    return sign * (intpart + Bernoulli(fracpart));
}

template <typename T>
void Randomizer::Shuffle(std::vector<T>& vec)
{
    gsl_ran_shuffle(r, &vec[0], vec.size(), sizeof(vec[0]));
}


#endif