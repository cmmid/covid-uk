#include <cmath>
#include <iostream>
#include <iomanip>
#include <vector>
#include <limits>
#include <string>
#include <stdexcept>
#include <nlopt.h>
//#include <omp.h>

template <typename Observables, typename LikelihoodFunc, typename ReportFunc>
void DEMCMC_Priors(Randomizer& R, LikelihoodFunc likelihood, ReportFunc report,
    int burn_in, int iterations, int n_chains, std::vector<Distribution>& priors,
    bool verbose = true, std::vector<std::string> param_names = std::vector<std::string>(),
    bool reeval_likelihood = false, bool in_parallel = false, int n_threads = -1, 
    bool classic_gamma = false, int R_skip = -1,
    std::vector<std::vector<double>> init = std::vector<std::vector<double>>(), int init_iter = 1)
{
    using namespace std;

// #ifdef _OPENMP
//     if (in_parallel && n_threads > 0)
//         omp_set_num_threads(n_threads);
// #endif

    if (n_chains < 3)
        throw runtime_error("Cannot use DE-MCMC with fewer than 3 chains.");

    // Store acceptance rates
    unsigned int ar_size = 1000;
    unsigned int ar_i = 0;
    bool show_ar = false;
    vector<bool> acceptances(ar_size, false);

    // Storage for chains and settings
    int n_theta = priors.size();
    vector<vector<double>> chains(n_chains, vector<double>(n_theta, 0.0));
    vector<Observables> obs(n_chains, Observables());
    vector<double> p(n_chains, 0.0);    // log probability for each chain
    vector<double> l(n_chains, 0.0);    // log likelihood for each chain
    double b = 0.001;

    // Storage for calls to Randomizer - to make thread-safe
    vector<vector<double>> random_perturbations(n_chains, vector<double>(n_theta, 0.0));
    vector<double> random_tests(n_chains, 0.0);
    vector<double> random_gammas(n_chains, 0.0);
    vector<vector<int>> random_chains(n_chains, vector<int>(2, 0));

    // Storage for calculation of Gelman-Rubin-Brooks R diagnostic
    int running_count = 0;
    vector<vector<double>> running_mean(n_chains, vector<double>(n_theta, 0.0));
    vector<vector<double>> running_M2(n_chains, vector<double>(n_theta, 0.0));
    vector<vector<double>> running_s2(n_chains, vector<double>(n_theta, 0.0));
    vector<double> R_hat(n_theta, 0.0);
    const int n_between_R = 100;

    // Assemble target func
    auto target = [&](vector<double>& theta, Observables& obs, double& l)
    {
        double p = 0;
        for (int d = 0; d < n_theta; ++d)
        {
            double pd = priors[d].LogProbability(theta[d]);
            if (pd == -std::numeric_limits<double>::infinity())
            {
                l = -std::numeric_limits<double>::infinity();
                return -std::numeric_limits<double>::infinity();
            }
            p += pd;
        }

        l = likelihood(theta, obs);

        return l + p;
    };

    // Initialize chains . . .
    if (init.empty())
    {
        // . . . from prior
        for (int c = 0; c < n_chains; ++c)
            for (int d = 0; d < n_theta; ++d)
                chains[c][d] = priors[d].RandomInit(R);
    }
    else
    {
        // . . . from initial values supplied
        if ((int)init.size() != n_chains || (int)init[0].size() != n_theta)
            throw runtime_error("init vector supplied is not the right size.");
        chains = init;
    }

    // Set initial probabilities and observables
    if (verbose)
        cout << "Initializing chains...\n";
    for (int c = 0; c < n_chains; ++c)
    {
        p[c] = target(chains[c], obs[c], l[c]);
        if (verbose)
            cout << "." << flush;
    }

    // Do iterations
    if (verbose)
        cout << "\nIterating...\n" << flush;
    for (int i = init_iter; i < burn_in + iterations; ++i)
    {
        // If requested, re-evaluate likelihood for next iteration
        if (reeval_likelihood)
        {
            // #pragma omp parallel for if(in_parallel) schedule(dynamic)
            for (int c = 0; c < n_chains; ++c)
                p[c] = target(chains[c], obs[c], l[c]);
        }

        // Prepare storage and random variates
        bool migration = i < burn_in * 0.75 ? R.Bernoulli(0.1) : false;
        vector<int> migration_indices(n_chains, 0);
        if (migration)
        {
            for (int c = 0; c < n_chains; ++c)
                migration_indices[c] = c;
            R.Shuffle(migration_indices);
        }
        for (int c = 0; c < n_chains; ++c)
        {
            for (int d = 0; d < n_theta; ++d)
                random_perturbations[c][d] = R.Uniform(-b, b);
            random_tests[c] = R.Uniform();
            if (!migration)
            {
                if (classic_gamma)
                    random_gammas[c] = (i % 10 == 0 ? 1.0 : 2.38 / sqrt(2 * n_theta));
                else
                    random_gammas[c] = R.Uniform(0.5, 1.0);
                do random_chains[c][0] = R.Discrete(n_chains); while (random_chains[c][0] == c);
                do random_chains[c][1] = R.Discrete(n_chains); while (random_chains[c][1] == c || random_chains[c][1] == random_chains[c][0]);
            }
        }
        auto saved_chains = chains;
        vector<int> accept(n_chains, 0);

        // #pragma omp parallel for if(in_parallel) schedule(dynamic)
        for (int c = 0; c < n_chains; ++c)
        {
            vector<double> theta_p = chains[c];
            int c_from = c;

            // Generate proposal, either by migration...
            if (migration)
            {
                c_from = migration_indices[c];
                theta_p = saved_chains[migration_indices[(c + 1) % n_chains]];
                for (int d = 0; d < n_theta; ++d)
                    theta_p[d] += random_perturbations[c][d];
            }
            else // ... or by directed mutation
            {
                for (int d = 0; d < n_theta; ++d)
                    theta_p[d] += random_gammas[c] * (saved_chains[random_chains[c][1]][d] - saved_chains[random_chains[c][0]][d]) + random_perturbations[c][d];
            }

            // Calculate log-probability and accept or reject
            Observables obs_p;
            double l_p = 0;
            double p_p = target(theta_p, obs_p, l_p);
            if ( (p_p == -std::numeric_limits<double>::infinity() && p[c_from] == -std::numeric_limits<double>::infinity() && random_tests[c] < 0.5)
                || (p_p > -std::numeric_limits<double>::infinity() && random_tests[c] < exp(p_p - p[c_from])) )
            {
                accept[c_from] = 1;
                chains[c_from] = theta_p;
                obs[c_from] = obs_p;
                p[c_from] = p_p;
                l[c_from] = l_p;
            }
        }

        // Update acceptances
        for (int c = 0; c < n_chains; ++c)
        {
            if (ar_i == ar_size - 1) show_ar = true;
            acceptances[ar_i] = accept[c];
            ar_i = (ar_i + 1) % ar_size;
        }

        // Update Gelman-Rubin-Brooks R
        bool R_all_ok = true;
        if (R_skip > 0)
        {
            ++running_count;
            for (int c = 0; c < n_chains; ++c)
            {
                for (int d = 0; d < n_theta; ++d)
                {
                    double delta = chains[c][d] - running_mean[c][d];
                    running_mean[c][d] += delta / running_count;
                    double delta2 = chains[c][d] - running_mean[c][d];
                    running_M2[c][d] += delta * delta2;
                }
            }

            // Calculate R every n_between_R generations
            if (i % n_between_R == 0)
            {
                // Finalise running mean and variance
                for (int c = 0; c < n_chains; ++c)
                    for (int d = 0; d < n_theta; ++d)
                        running_s2[c][d] = running_M2[c][d] / (running_count - 1);

                // Calculate statistic for each parameter
                for (int d = 0; d < n_theta; ++d)
                {
                    double M = n_chains;
                    double N = running_count;
                    double W = 0, X = 0, B = 0;
                    for (int c = 0; c < n_chains; ++c)
                    {
                        W += running_s2[c][d];
                        X += running_mean[c][d];
                    }
                    W /= M;
                    X /= M;

                    for (int c = 0; c < n_chains; ++c)
                        B += (running_mean[c][d] - X) * (running_mean[c][d] - X);
                    B *= N / (M - 1);

                    double var = ((N - 1) / N) * W + B / N;
                    R_hat[d] = std::sqrt(var / W);

                    if (R_hat[d] > 1.05)
                        R_all_ok = false;
                }
            }
        }

        // Report results of this iteration
        for (int c = 0; c < n_chains; ++c)
            report(i - burn_in, p[c], c, l[c], chains[c], obs[c]);

        // Print progress
        if (verbose)
        {
            cout << "." << flush;
            if (i % 100 == 0)
            {
                cout << "\n" << (i < burn_in ? "burn-in" : "main") << " iteration " << i - burn_in << ":";
                if (!param_names.empty())
                {
                    cout << "\n         " << setw(12) << right << "log (P)" << setw(12) << right << "log (L)";
                    for (auto n : param_names)
                        cout << setw(12) << right << n;
                }
                for (int c = 0; c < n_chains; ++c)
                {
                    cout << "\nchain" << setw(4) << right << c << setw(12) << right << p[c] << setw(12) << right << l[c];
                    for (int d = 0; d < n_theta; ++d)
                        cout << setw(12) << right << chains[c][d];
                }
                if (R_skip > 0)
                {
                    cout << "\nGelman-Rubin-Brooks diagnostic R ";
                    for (int d = 0; d < n_theta; ++d)
                        cout << setw(12) << right << R_hat[d];
                }

                double acceptance_rate = show_ar ? (double)count(acceptances.begin(), acceptances.end(), true) / ar_size : -1;
                cout << "\nacceptance rate: " << acceptance_rate << "\n\n" << flush;
            }
        }

        // Skip past burn-in if R OK
        if (R_skip > 0 && i > R_skip && i < burn_in && R_all_ok)
        {
            if (verbose)
                cout << "\n\nSkipping to iterations (R < 1.05 reached).\n\n";
            i = burn_in - 1;
        }

    }
    if (verbose)
        cout << "\n";
}

template <typename Observables, typename PFunc>
struct OptimizeData
{
    PFunc target;
    Observables obs;
    std::vector<Distribution> priors;
    bool verbose;
    unsigned int count;

    OptimizeData(PFunc targ)
     : target(targ) { }

    static double Test(unsigned int n, const double* thet, double* grad, void* f_data)
    {
        (void) grad;
        auto& THIS = *(OptimizeData<Observables, PFunc>*)(f_data);
        ++THIS.count;
        std::vector<double> theta(thet, thet + n);

        // Assemble prior probability
        double p = 0;
        for (unsigned int d = 0; d < theta.size(); ++d)
        {
            double pd = THIS.priors[d].LogProbability(theta[d]);
            if (pd == std::numeric_limits<double>::lowest())
                return std::numeric_limits<double>::lowest();
            p += pd;
        }

        // Get target probability (i.e. likelihood)
        double t = THIS.target(theta, THIS.obs);

        if (THIS.verbose)
        {
            if (THIS.count % 100 == 0)
            {
                std::cout << "Optimize step " << THIS.count << ": ";
                for (auto th : theta)
                    std::cout << " " << th;
                std::cout << ", ll = " << t << ", lp = " << t + p << ".\n";
            }
        }

        return t + p;
    }
};


template <typename Observables, typename LikelihoodFunc, typename ReportFunc>
void Optimize_Priors(Randomizer& R, LikelihoodFunc likelihood, ReportFunc report, std::vector<Distribution>& priors,
    bool global, std::string global_algorithm, unsigned int global_maxeval, double global_ftol_abs,
    bool local, std::string local_algorithm, unsigned int local_maxeval, double local_ftol_abs, bool verbose,
    std::vector<double>* initial = 0, double stepsize = -1)
{
    using namespace std;

    nlopt_srand(R.Discrete(std::numeric_limits<unsigned int>::max()));

    int n_theta = priors.size();
    vector<double> lb, ub;
    for (unsigned int i = 0; i < priors.size(); ++i)
    {
        lb.push_back(priors[i].LowerBound());
        ub.push_back(priors[i].UpperBound());
    }

    OptimizeData<Observables, LikelihoodFunc> f_data(likelihood);
    f_data.priors = priors;
    f_data.verbose = verbose;
    f_data.count = 0;
    vector<double> theta(n_theta, 0.0);

    if (initial)
        theta = *initial;
    else for (int d = 0; d < n_theta; ++d)
    {
        theta[d] = priors[d].RandomInit(R);
    }

    double p;
    int error = 0;

    if (global)
    {
        if (verbose)
            cout << "Global optimization:\n";

        nlopt_opt go;
        if (global_algorithm == "DIRECT")
            go = nlopt_create(NLOPT_GN_DIRECT, n_theta);
        else if (global_algorithm == "DIRECT_L")
            go = nlopt_create(NLOPT_GN_DIRECT_L, n_theta);
        else if (global_algorithm == "DIRECT_L_RAND")
            go = nlopt_create(NLOPT_GN_DIRECT_L_RAND, n_theta);
        else if (global_algorithm == "DIRECT_NOSCAL")
            go = nlopt_create(NLOPT_GN_DIRECT_NOSCAL, n_theta);
        else if (global_algorithm == "DIRECT_L_NOSCAL")
            go = nlopt_create(NLOPT_GN_DIRECT_L_NOSCAL, n_theta);
        else if (global_algorithm == "DIRECT_L_RAND_NOSCAL")
            go = nlopt_create(NLOPT_GN_DIRECT_L_RAND_NOSCAL, n_theta);
        else if (global_algorithm == "ORIG_DIRECT")
            go = nlopt_create(NLOPT_GN_ORIG_DIRECT, n_theta);
        else if (global_algorithm == "ORIG_DIRECT_L")
            go = nlopt_create(NLOPT_GN_ORIG_DIRECT_L, n_theta);
        else if (global_algorithm == "CRS2_LM")
            go = nlopt_create(NLOPT_GN_CRS2_LM, n_theta);
        else if (global_algorithm == "ISRES")
            go = nlopt_create(NLOPT_GN_ISRES, n_theta);
        else if (global_algorithm == "ESCH")
            go = nlopt_create(NLOPT_GN_ESCH, n_theta);
        else
            throw runtime_error("global_algorithm must be DIRECT, DIRECT_L, DIRECT_L_RAND, DIRECT_NOSCAL, DIRECT_L_NOSCAL, DIRECT_L_RAND_NOSCAL, ORIG_DIRECT, ORIG_DIRECT_L, CRS2_LM, ISRES, or ESCH. See NLopt documentation online for details.");

        nlopt_set_max_objective(go, OptimizeData<Observables, LikelihoodFunc>::Test, &f_data);
        nlopt_set_lower_bounds(go, &lb[0]);
        nlopt_set_upper_bounds(go, &ub[0]);

        nlopt_set_maxeval(go, global_maxeval);
        nlopt_set_ftol_abs(go, global_ftol_abs);

        if (nlopt_optimize(go, &theta[0], &p) < 0)
            cout << "NLOPT global optimization failed.\n";
        nlopt_destroy(go);
    }

    if (local)
    {
        if (verbose)
            cout << "Local refinement:\n";

        nlopt_opt lo;
        if (local_algorithm == "COBYLA")
            lo = nlopt_create(NLOPT_LN_COBYLA, n_theta);
        else if (local_algorithm == "BOBYQA")
            lo = nlopt_create(NLOPT_LN_BOBYQA, n_theta);
        else if (local_algorithm == "NEWUOA_BOUND")
            lo = nlopt_create(NLOPT_LN_NEWUOA_BOUND, n_theta);
        else if (local_algorithm == "PRAXIS")
            lo = nlopt_create(NLOPT_LN_PRAXIS, n_theta);
        else if (local_algorithm == "NELDERMEAD")
            lo = nlopt_create(NLOPT_LN_NELDERMEAD, n_theta);
        else if (local_algorithm == "SBPLX")
            lo = nlopt_create(NLOPT_LN_SBPLX, n_theta);
        else
            throw runtime_error("local_algorithm must be COBYLA, BOBYQA, NEWUOA_BOUND, PRAXIS, NELDERMEAD, or SBPLX. See NLopt documentation online for details.");

        nlopt_set_max_objective(lo, OptimizeData<Observables, LikelihoodFunc>::Test, &f_data);
        nlopt_set_lower_bounds(lo, &lb[0]);
        nlopt_set_upper_bounds(lo, &ub[0]);

        nlopt_set_maxeval(lo, local_maxeval);
        nlopt_set_ftol_abs(lo, local_ftol_abs);

        if (stepsize > 0)
            nlopt_set_initial_step1(lo, stepsize);

        do {
            error = 0;
            if (nlopt_optimize(lo, &theta[0], &p) < 0)
            {
                cout << "NLOPT local optimization failed.\n";
                error = -1;
                cout << "Trying different step size.\n";
                nlopt_set_initial_step1(lo, 0.001 * R.Exponential(1.0));
            }
        } while (error < 0);
    }

    double lprior = 0;
    for (unsigned int d = 0; d < theta.size(); ++d)
        lprior += priors[d].LogProbability(theta[d]);

    report(-f_data.count, p, error, p - lprior, theta, f_data.obs);
}

