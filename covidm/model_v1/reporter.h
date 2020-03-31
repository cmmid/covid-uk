// reporter.h

// For reporting results
class Reporter
{
public:
    Reporter(Parameters& P);

    // Access / modify data
    double& operator()(double t, unsigned int p, unsigned int a, unsigned int c)
    {
        unsigned int row = (unsigned int)(t - t0) * n_populations * n_age_groups + p * n_age_groups + a;
        return data[c][row];
    }

//private:
    double t0;
    unsigned int n_times;
    unsigned int n_populations;
    unsigned int n_age_groups;
    vector<string> col_names;

    vector<Rcpp::NumericVector> data;
    Rcpp::DataFrame dynamics_df;
    string csv;
};

Reporter::Reporter(Parameters& P)
 : t0(P.time0),
   n_times((unsigned int)(P.time1 - P.time0) + 1),
   n_populations(P.pop.size()),
   n_age_groups(P.pop[0].size.size()),
   col_names({ "S", "E", "Ip", "Is", "Ia", "R", "cases", "cases_reported", "subclinical" })
{
    if (P.time_step != 1. / P.report_every)
        throw logic_error("Reporter requires P.time_step = 1 / P.report_every.");

    // Create times
    Rcpp::NumericVector t(n_times * n_populations * n_age_groups, 0.);
    for (unsigned int it = 0; it < n_times; ++it)
        for (unsigned int j = 0; j < n_populations * n_age_groups; ++j)
            t[it * n_populations * n_age_groups + j] = P.time0 + it * P.time_step * P.report_every;

    // Create identifier columns
    dynamics_df = Rcpp::DataFrame::create(
        Rcpp::Named("t") = t,
        Rcpp::Named("population") = Rcpp::rep(Rcpp::rep_each(Rcpp::seq(1, n_populations), n_age_groups), n_times),
        Rcpp::Named("group") = Rcpp::rep(Rcpp::seq(1, n_age_groups), n_times * n_populations)
    );

    // Create space for built-in compartments
    for (unsigned int c = 0; c < col_names.size(); ++c)
        data.push_back(Rcpp::NumericVector(n_times * n_populations * n_age_groups, 0.));

    // User-defined compartments
    for (auto& pr : P.processes)
    {
        for (unsigned int j = 0; j < pr.names.size(); ++j)
        {
            for (unsigned int k = 0; k < pr.report[j].size(); ++k)
            {
                col_names.push_back(pr.names[j] + "_" + pr.report[j][k]);
                data.push_back(Rcpp::NumericVector(n_times * n_populations * n_age_groups, 0.));

                if (pr.report[j][k] == 'p') {
                    pr.p_cols.push_back(data.size() - 1);
                    pr.p_ids.push_back(pr.ids[j]);
                } else if (pr.report[j][k] == 'i') {
                    pr.i_cols.push_back(data.size() - 1);
                    pr.i_ids.push_back(pr.ids[j]);
                } else if (pr.report[j][k] == 'o') {
                    pr.o_cols.push_back(data.size() - 1);
                    pr.o_ids.push_back(pr.ids[j]);
                } else {
                    throw runtime_error("Unrecognized process report type " + string(1, pr.report[j][k]) + ".");
                }
            }
        }
    }

    // Allocate all columns to the dataframe
    for (unsigned int c = 0; c < col_names.size(); ++c)
    {
        dynamics_df.push_back(data[c], col_names[c]);
    }

    // Set dataframe as a data.table
    Rcpp::Function setDT("setDT"); 
    setDT(dynamics_df);
}
