// Rcpp_interface.h

// [[Rcpp::export]]
Rcpp::List cm_backend_simulate(Rcpp::List parameters, unsigned int n_run = 1, unsigned long int seed = 0)
{
    // Initialise parameters for this simulation
    Randomizer Rand(seed);
    Parameters covidm_parameters;
    SetParameters(covidm_parameters, parameters, Rand);

    Rcpp::List dynamics;
    Rcpp::List csvs;

    for (unsigned int r = 0; r < n_run; ++r)
    {
        // Run the simulation
        Parameters P = covidm_parameters;
        Reporter rep = RunSimulation(P, Rand);

        dynamics.push_back(rep.dynamics_df);
        csvs.push_back(rep.csv);
    }

    return Rcpp::List::create(
        Rcpp::Named("dynamics") = dynamics,
        Rcpp::Named("csv") = csvs
    );
}

// [[Rcpp::export]]
Rcpp::DataFrame cm_evaluate_distribution(string dist_code, unsigned int steps = 101, double xmin = 0, double xmax = -1)
{
    Distribution dist(dist_code);
    if (xmax < xmin)
    {
        xmin = dist.LowerBound();
        xmax = dist.UpperBound();
    }

    vector<double> x(steps, 0.);
    vector<double> p(steps, 0.);

    for (unsigned int s = 0; s < steps; ++s)
    {
        x[s] = xmin + ((xmax - xmin) / (steps - 1.)) * s;
        p[s] = exp(dist.LogProbability(x[s]));
    }

    Rcpp::DataFrame results = Rcpp::DataFrame::create(
        Rcpp::Named("x") = x,
        Rcpp::Named("p") = p
    );

    return results;
}

// ]]Rcpp::export("cm_backend_particle_likelihood")[[
// double ParticleLikelihood_R(Rcpp::List parameters, unsigned int n_particles)
// {
//     Randomizer Rand;
//     Parameters P;
//     SetParameters(P, parameters, Rand);
//     auto l = ParticleLikelihood(P, Rand, n_particles);
//     return l;
// }

// double ParticleLikelihood(Parameters& P, Randomizer& Rand, unsigned int n_particles)
// {
//     // Metapopulations as particles
//     vector<Metapopulation> mpp(n_particles, Metapopulation(P));
//     vector<Metapopulation> mpp_other(n_particles, Metapopulation(P));
//     vector<Reporter> report(n_particles, Reporter(P, false));
//     vector<Reporter> report_other(n_particles, Reporter(P, false));
//     vector<double> weights(n_particles);
//     double ll = 0;
//     double max_lwt = 0;
    
//     vector<double> n_00_06(n_particles, 0.0);
//     vector<double> n_07_17(n_particles, 0.0);
//     vector<double> n_18_24(n_particles, 0.0);
//     vector<double> n_25_49(n_particles, 0.0);
//     vector<double> n_50_64(n_particles, 0.0);
//     vector<double> n_65plu(n_particles, 0.0);
//     vector<unsigned int> picks(n_particles, 0);

//     auto weight_agedist = [&](const unsigned int observed[])
//     {
//         max_lwt = std::numeric_limits<double>::lowest();

//         for (unsigned int i = 0; i < n_particles; ++i)
//         {
//             auto cas = report[i].true_cases[0];
//             double n_00_06 = cas[0] + 0.4 * cas[1];
//             double n_07_17 = 0.6 * cas[1] + cas[2] + 0.6 * cas[3];
//             double n_18_24 = 0.4 * cas[3] + cas[4];
//             double n_25_49 = cas[5] + cas[6] + cas[7] + cas[8] + cas[9];
//             double n_50_64 = cas[10] + cas[11] + cas[12];
//             double n_65plu = cas[13] + cas[14] + cas[15] + cas[16] + cas[17];
            
//             double freq[] = { 
//                 n_00_06 + 1, 
//                 n_07_17 + 1, 
//                 n_18_24 + 1, 
//                 n_25_49 + 1, 
//                 n_50_64 + 1, 
//                 n_65plu + 1
//              };

//             double lwt = gsl_ran_multinomial_lnpdf(6, freq, observed);

//             ll += lwt;
//             max_lwt = max(lwt, max_lwt);
//             weights[i] = lwt;
//             report[i].true_cases[0].assign(18, 0.0);
//         }
//     };

//     auto weight_ncases = [&](const double observed)
//     {
//         max_lwt = std::numeric_limits<double>::lowest();

//         // weight
//         for (unsigned int i = 0; i < n_particles; ++i)
//         {
//             auto cases = max(1.0, accumulate(report[i].rep_cases[0].begin(), report[i].rep_cases[0].end(), 0.0));
//             double lwt = observed * log(cases) - cases - gsl_sf_lngamma(observed + 1);

//             ll += lwt;
//             max_lwt = max(lwt, max_lwt);
//             weights[i] = lwt;
//         }
//     };

//     auto resample = [&]()
//     {
//         // normalise weights such that max weight is 1
//         for (unsigned int i = 0; i < n_particles; ++i)
//             weights[i] = exp(weights[i] - max_lwt);

//         // resample
//         picks.assign(n_particles, 0);
//         Rand.Multinomial(n_particles, weights, picks);

//         unsigned int newind = 0;
//         for (unsigned int i = 0; i < n_particles; ++i)
//         {
//             for (unsigned int j = 0; j < picks[i]; ++j)
//             {
//                 mpp_other[newind] = mpp[i];
//                 report_other[newind] = report[i];
//                 ++newind;
//             }
//         }
//         swap(mpp, mpp_other);
//         swap(report, report_other);
//     };

//     unsigned int time_steps = (P.time1 - P.time0) / P.time_step;
//     const unsigned int counts_43[] = { 0, 1, 1, 21, 20, 26 };
//     const unsigned int counts_51[] = { 0, 1, 14, 119, 71, 32 };
//     const unsigned int counts_64[] = { 4, 14, 55, 322, 104, 38 };
//     for (unsigned int ts = 0; ts < time_steps; ++ts)
//     {
//         // Update metapopulations
//         double t = P.time0 + ts * P.time_step;
//         for (unsigned int m = 0; m < mpp.size(); ++m)
//             mpp[m].Tick(P, Rand, t, ts, report[m]);

//         // HACK
//         if (t == 43) {
//             weight_agedist(counts_43);
//             resample();
//         } else if (t == 51) {
//             weight_agedist(counts_51);
//             resample();
//         } else if (t == 55) {
//             weight_ncases(571);
//             resample();
//         } else if (t == 62) {
//             weight_ncases(7736);
//             resample();
//         } else if (t == 64) {
//             weight_agedist(counts_64);
//             resample();
//         } else if (t == 69) {
//             weight_ncases(28060);
//             resample();
//         }
//     }

//     return ll;
// }

/* part filter algo

1. For each particle j = 1 … J
    2. initialise the state of particle j
    3. initialise the weight of particle j
4. For each observation time t = 1 … T
    5. resample particles
    6. For each particle j = 1 … J
        7. propagate particle j to next observation time
        8. weight particle j

*/
