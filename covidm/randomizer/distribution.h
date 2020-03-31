// distribution.h

#ifndef DISTRIBUTION_H
#define DISTRIBUTION_H

#include "randomizer.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_cdf.h>
#include <algorithm>
#include <string>
#include <iostream>

class Distribution
{
public:
    enum Type { Cauchy, LogNormal, Normal, Uniform, Beta, Exponential, Gamma, Rounded };
    static const std::string Types[];

    Distribution(Type t, std::vector<double> parameters, std::vector<double> trunc = {}, std::vector<double> shift = {}, std::vector<double> init = {});
    Distribution(std::string code);
    Distribution(std::vector<std::string> subcode, std::string code = "");

    double Random(Randomizer& R);
    double RandomInit(Randomizer& R);
    double LogProbability(double x);
    double LowerBound();
    double UpperBound();

    friend std::ostream& operator << (std::ostream& out, const Distribution& d)
    {
        out << Types[d.type] << "(" << d.a << ", " << d.b << ")";
        if (d.t_lp_adj != 0)
            out << " T [" << d.t0 << ", " << d.t1 << "] adj " << d.t_lp_adj;
        if (d.s0 != 0 || d.s1 != 1)
            out << " * " << d.s1 << " + " << d.s0;
        return out;
    }

private:
    static constexpr double ShoulderArea = 1e-20;

    void Init(Type t, std::vector<double> parameters, std::vector<double> trunc = {}, std::vector<double> shift = {}, std::vector<double> init = {});
    void SetTrunc(double xmin, double xmax);
    void SetFromCode(std::vector<std::string> subcode, std::string code);
    static std::vector<std::string> Unserialize(std::string s);
    static double RoundedUniformCDF(double min, double max, double shoulder, double x);

    Type type;
    double a, b;
    double t0, t1, t_lp_adj;
    double s0, s1;
    double i0, i1;
};

const std::string Distribution::Types[] = { "Cauchy", "LogNormal", "Normal", "Uniform", "Beta", "Exponential", "Gamma", "Rounded" };

Distribution::Distribution(Type t, std::vector<double> parameters, std::vector<double> trunc, std::vector<double> shift, std::vector<double> init)
{
    Init(t, parameters, trunc, shift, init);
}

Distribution::Distribution(std::string code)
{
    SetFromCode(Unserialize(code), code);
}

Distribution::Distribution(std::vector<std::string> subcode, std::string code)
{
    SetFromCode(subcode, code);
}

double Distribution::Random(Randomizer& R)
{
    double x;
    do
    {
        switch (type)
        {
            case Cauchy:
                x = R.Cauchy(a, b); break;
            case LogNormal:
                x = R.LogNormal(a, b); break;
            case Normal:
                x = R.Normal(a, b); break;
            case Uniform:
                x = R.Uniform(a, b); break;
            case Beta:
                x = R.Beta(a, b); break;
            case Exponential:
                x = R.Exponential(a); break;
            case Gamma:
                x = R.Gamma(a, b); break;
            case Rounded:
                x = R.RoundedUniform(a, b, ShoulderArea); break;
            default:
                throw std::runtime_error("No distribution type in Distribution::Random.");
        }
        x = s0 + x * s1;
    } while (x < t0 || x > t1);

    return x;
}

double Distribution::RandomInit(Randomizer& R)
{
    if (i1 < i0)
        return Random(R);
    else
        return R.Uniform(i0, i1);
}

double Distribution::LogProbability(double x)
{
    using std::log;

    if (x < t0 || x > t1)
        return -std::numeric_limits<double>::infinity();

    // TODO this can quite conceivably fail when x is close to the limit of a bounded distribution which is scaled or translated.
    // ideally, x would be slightly moved inside so long as this was all valid.
    x = (x - s0) / s1;
    double d;

    switch (type)
    {
        case Cauchy:
            d = -log(b * (1. + ((x - a) / b) * ((x - a) / b))) - 1.1447298858494;
            break;
        case LogNormal:
            d = -(log(x) - a) * (log(x) - a) / (2 * b * b) - log(x * b) - 0.918938533204673;
            break;
        case Normal:
            d = -(x - a) * (x - a) / (2 * b * b) - log(b) - 0.918938533204673;
            break;
        case Uniform:
            if (b == a) d = 0;
            else d = -log(b - a);
            break;
        case Beta:
            if (a == 1 && b == 1) d = 0; // needed because at x = 0 or x = 1 when alpha = beta = 1, below formula fails
            else d = gsl_sf_lngamma(a + b) - gsl_sf_lngamma(a) - gsl_sf_lngamma(b) + (a - 1) * log(x) + (b - 1) * log(1 - x);
            break;
        case Exponential:
            d = -log(a) - a * x;
            break;
        case Gamma:
            d = -gsl_sf_lngamma(a) - a * log(b) + (a - 1) * log(x) - x / b;
            break;
        case Rounded:
            if (x < a)      { double sd = ShoulderArea * (b - a) / ((1 - ShoulderArea) * 2.50662827463); d = log(ShoulderArea) - (x - a) * (x - a) / (2 * sd * sd) - log(sd) - 0.918938533204673; }
            else if (x < b) { d = log((1 - ShoulderArea) / (b - a)); }
            else            { double sd = ShoulderArea * (b - a) / ((1 - ShoulderArea) * 2.50662827463); d = log(ShoulderArea) - (x - b) * (x - b) / (2 * sd * sd) - log(sd) - 0.918938533204673; }
            break;
        default:
            throw std::runtime_error("No distribution type in Distribution::LogProbability.");
    }

    return d + t_lp_adj;
}

double Distribution::LowerBound()
{
    return t0;
}

double Distribution::UpperBound()
{
    return t1;
}

void Distribution::Init(Type t, std::vector<double> parameters, std::vector<double> trunc, std::vector<double> shift, std::vector<double> init)
{
    type = t;

    auto set = [](double& a, double& b, double A, double B, std::vector<double> v) {
        a = v.size() > 0 ? v[0] : A;
        b = v.size() > 1 ? v[1] : B;
    };

    set(a, b, t <= Uniform ? 0 : 1, 1, parameters);
    set(t0, t1, std::numeric_limits<double>::lowest(), std::numeric_limits<double>::max(), trunc);
    set(s0, s1, 0., 1., shift);
    set(i0, i1, 0., -1., init);

    SetTrunc(t0, t1);
}

void Distribution::SetTrunc(double xmin, double xmax)
{
    // First, adjust xmin and xmax to be compatible with the underlying distribution
    double natural_xmin, natural_xmax, sign = s1 < 0 ? -1 : 1;
    switch (type)
    {
        case Cauchy:
        case Normal:
        case Rounded:
            natural_xmin = std::numeric_limits<double>::lowest();
            natural_xmax = std::numeric_limits<double>::max();
            break;
        case LogNormal:
            natural_xmin = s0 + sign * std::numeric_limits<double>::min();
            natural_xmax = sign * std::numeric_limits<double>::max();
            break;
        case Exponential:
        case Gamma:
            natural_xmin = s0;
            natural_xmax = sign * std::numeric_limits<double>::max();
            break;
        case Uniform:
            natural_xmin = s0 + a * s1;
            natural_xmax = s0 + b * s1;
            break;
        case Beta:
            natural_xmin = s0 + 0 * s1;
            natural_xmax = s0 + 1 * s1;
            break;

        default:
            throw std::runtime_error("No distribution type in Distribution::SetTrunc.");
    }

    if (natural_xmin > natural_xmax)
        std::swap(natural_xmin, natural_xmax);

    t0 = std::max(xmin, natural_xmin);
    t1 = std::min(xmax, natural_xmax);

    // Only calculate CDF adjustment if the underlying distribution is being cut off
    if (xmin <= natural_xmin && xmax >= natural_xmax)
    {
        t_lp_adj = 0;
    }

    // Calculate the log-probability adjustment according to how much of the underlying distribution is being cut off
    else
    {
        double q0 = (t0 - s0) / s1, q1 = (t1 - s0) / s1;
        if (q0 > q1)
            std::swap(q0, q1);
        double cdf0, cdf1;

        switch (type)
        {
            case Cauchy:
                cdf0 = gsl_cdf_cauchy_P(q0 - a, b);
                cdf1 = gsl_cdf_cauchy_P(q1 - a, b);
                break;
            case LogNormal:
                cdf0 = gsl_cdf_lognormal_P(q0, a, b);
                cdf1 = gsl_cdf_lognormal_P(q1, a, b);
                break;
            case Normal:
                cdf0 = gsl_cdf_gaussian_P(q0 - a, b);
                cdf1 = gsl_cdf_gaussian_P(q1 - a, b);
                break;
            case Uniform:
                cdf0 = gsl_cdf_flat_P(q0, a, b);
                cdf1 = gsl_cdf_flat_P(q1, a, b);
                break;
            case Beta:
                cdf0 = gsl_cdf_beta_P(q0, a, b);
                cdf1 = gsl_cdf_beta_P(q1, a, b);
                break;
            case Exponential:
                cdf0 = gsl_cdf_exponential_P(q0, 1. / a);
                cdf1 = gsl_cdf_exponential_P(q1, 1. / a);
                break;
            case Gamma:
                cdf0 = gsl_cdf_gamma_P(q0, a, b);
                cdf1 = gsl_cdf_gamma_P(q1, a, b);
                break;
            case Rounded:
                cdf0 = RoundedUniformCDF(a, b, ShoulderArea, q0);
                cdf1 = RoundedUniformCDF(a, b, ShoulderArea, q1);
                break;
            default:
                throw std::runtime_error("No distribution type in Distribution::SetTrunc.");
        }
        t_lp_adj = -std::log(cdf1 - cdf0);
    }
}

void Distribution::SetFromCode(std::vector<std::string> subcode, std::string code)
{
    auto trunc_start = std::find(subcode.begin(), subcode.end(), std::string("T"));
    auto scale_start = std::find(subcode.begin(), subcode.end(), std::string("S"));
    auto init_start = std::find(subcode.begin(), subcode.end(), std::string("I"));

    if (subcode.size() < 1)
        throw std::runtime_error("Could not interpret code " + code + " in Distribution::SetFromCode.");
    else if (subcode.size() == 1)
    {
        try {
            double value = stod(subcode[0]);
            Init(Uniform, { value, value + 1e-6 }, {}, {}, {});
        } catch (...) {
            throw std::runtime_error("Could not interpret single-value code " + code + " as double.");
        }
        return;
    }

    Type dist;
    if      (subcode[0] == "C") dist = Cauchy;
    else if (subcode[0] == "L") dist = LogNormal;
    else if (subcode[0] == "N") dist = Normal;
    else if (subcode[0] == "U") dist = Uniform;
    else if (subcode[0] == "B") dist = Beta;
    else if (subcode[0] == "E") dist = Exponential;
    else if (subcode[0] == "G") dist = Gamma;
    else if (subcode[0] == "R") dist = Rounded;
    else throw std::runtime_error("Unrecognised distribution code " + subcode[0] + " in Distribution::SetFromCode.");

    std::vector<double> parameters, trunc, scale, init;

    for (auto i = subcode.begin() + 1; i < subcode.end() && i != trunc_start && i != scale_start && i != init_start; ++i)
        parameters.push_back(std::stod(*i));
    for (auto i = trunc_start + 1; i < subcode.end() && i != scale_start && i != init_start; ++i)
        trunc.push_back(std::stod(*i));
    for (auto i = scale_start + 1; i < subcode.end() && i != trunc_start && i != init_start; ++i)
        scale.push_back(std::stod(*i));
    for (auto i = init_start + 1; i < subcode.end() && i != trunc_start && i != scale_start; ++i)
        init.push_back(std::stod(*i));

    Init(dist, parameters, trunc, scale, init);
}

std::vector<std::string> Distribution::Unserialize(std::string s)
{
    std::string sep = " \f\n\r\t\v";
    std::vector<std::string> ret;
    std::string::size_type n0 = s.find_first_not_of(sep), n1;
    while ((n1 = s.find_first_of(sep, n0)) != std::string::npos)
    {
        ret.push_back(s.substr(n0, n1 - n0));
        n0 = s.find_first_not_of(sep, n1);
    }
    if (n0 != std::string::npos)
        ret.push_back(s.substr(n0));
    return ret;
}

double Distribution::RoundedUniformCDF(double min, double max, double shoulder, double x)
{
    double sd = shoulder * (max - min) / ((1 - shoulder) * 2.50662827463);

    if (x < min)        return gsl_cdf_gaussian_P(x - min, sd) * shoulder;
    else if (x < max)   return shoulder / 2 + (x - min) / (max - min) * (1 - shoulder);
    else                return gsl_cdf_gaussian_P(x - max, sd) * shoulder + (1 - shoulder);
}

#endif