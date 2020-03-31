// parameters.h

#include <Rcpp.h>

class Reporter;
struct PopulationParameters;
struct Parameters;

// Observer function
struct Observer
{
    Observer() : func("rnorm"), null(true) { } // TODO CHANGE THIS TO STH ELSE
    Observer& operator=(Rcpp::Function f) { func = f; null = false; return *this; }
    bool operator()(Parameters* parent, PopulationParameters& params, double t, Reporter& rep);
    Rcpp::Function func;
    bool null;
};

// Schedule entry
struct ScheduleEntry
{
    double t;
    string variable;
    Rcpp::RObject value;
};

// Helpers to set parameters
#define _CheckSet(variable) else if (name == #variable) { ParamSet(variable, value); }
#define _CheckSetParent(P, variable) else if (P != 0 && name == #variable) { ParamSet(P->variable, value); }

void ParamSet(Discrete& variable, Rcpp::RObject& value)
{
    variable = Rcpp::as<vector<double>>(value);
}

void ParamSet(vector<double>& variable, Rcpp::RObject& value)
{
    variable = Rcpp::as<vector<double>>(value);
}

void ParamSet(Matrix& variable, Rcpp::RObject& value)
{
    Rcpp::NumericMatrix rhs = Rcpp::as<Rcpp::NumericMatrix>(value);
    variable = Matrix(0.0, rhs.nrow(), rhs.ncol());
    for (int r = 0; r < rhs.nrow(); ++r)
        for (int c = 0; c < rhs.ncol(); ++c)
            variable(r, c) = rhs(r, c);
}

void ParamSet(vector<Matrix>& variable, Rcpp::RObject& value)
{
    Rcpp::List ml = Rcpp::as<Rcpp::List>(value);
    for (unsigned int i = 0; i < ml.size(); ++i)
    {
        Rcpp::RObject mat = Rcpp::as<Rcpp::RObject>(ml[i]);
        ParamSet(variable[i], mat);
    }
}

void ParamSet(vector<ProcessSpec>& variable, Rcpp::RObject& value)
{
    Rcpp::List pl = Rcpp::as<Rcpp::List>(value);
    unsigned int pc_id = 0;
    vector<string> pc_names;
    for (unsigned int i = 0; i < pl.size(); ++i)
    {
        Rcpp::List pli = Rcpp::as<Rcpp::List>(pl[i]);

        ProcessSpec process;
        process.source_name = Rcpp::as<string>(pli["source"]);
        process.type = Rcpp::as<string>(pli["type"]);

        process.names = Rcpp::as<vector<string>>(pli["names"]);
        process.ids = vector<unsigned int>(process.names.size(), 0);
        for (unsigned int j = 0; j < process.ids.size(); ++j)
        {
            if (process.names[j] == "null")
            {
                process.ids[j] = Null;
            }
            else 
            {
                process.ids[j] = pc_id++;
                pc_names.push_back(process.names[j]);
            }
        }
        process.report = Rcpp::as<vector<string>>(pli["report"]);

        Matrix m_prob, m_delays;
        Rcpp::RObject r_prob = Rcpp::as<Rcpp::RObject>(pli["prob"]);
        Rcpp::RObject r_delays = Rcpp::as<Rcpp::RObject>(pli["delays"]);
        ParamSet(m_prob, r_prob);
        ParamSet(m_delays, r_delays);

        for (unsigned int c = 0; c < m_prob.NCol(); ++c)
        {
            process.prob.push_back(vector<double>(m_prob.NRow(), 0.));
            for (unsigned int r = 0; r < m_prob.NRow(); ++r)
                process.prob[c][r] = m_prob(r, c);
        }

        for (unsigned int r = 0; r < m_delays.NRow(); ++r)
        {
            process.delays.push_back(Discrete());
            std::vector<double> uw(m_delays.NCol(), 0.);
            for (unsigned int c = 0; c < m_delays.NCol(); ++c)
                uw[c] = m_delays(r, c);
            process.delays.back() = uw;
        }

        variable.push_back(process);
    }

    // Set source_ids
    for (auto& pr : variable)
    {
        auto sn = std::find(pc_names.begin(), pc_names.end(), pr.source_name);
        if (sn == pc_names.end())
        {
            if (pr.source_name == "S")
                pr.source_id = srcS;
            else if (pr.source_name == "E")
                pr.source_id = srcE;
            else if (pr.source_name == "E_Ip")
                pr.source_id = srcEp;
            else if (pr.source_name == "E_Ia")
                pr.source_id = srcEa;
            else if (pr.source_name == "Ip")
                pr.source_id = srcIp;
            else if (pr.source_name == "Is")
                pr.source_id = srcIs;
            else if (pr.source_name == "H")
                pr.source_id = srcH;
            else if (pr.source_name == "Ia")
                pr.source_id = srcIa;
            else if (pr.source_name == "I")
                pr.source_id = srcI;
            else
                throw logic_error("Unrecognized process source name " + pr.source_name);
        }
        else
        {
            pr.source_id = (unsigned int)(sn - pc_names.begin());
        }
    }
}

// Population-level parameters for batch of simulations
struct PopulationParameters
{
public:
    PopulationParameters() : needs_recalc(true) {}

    bool needs_recalc;  // does contact matrix need recalculation?
    Matrix cm;          // contact matrix

    Discrete dE;
    Discrete dIp;       // TODO: any need for these to be age-specific?
    Discrete dIa;
    Discrete dIs;
    Discrete dH;
    Discrete dC;

    vector<double> size;
    vector<Matrix> matrices;
    vector<double> contact;
    vector<double> contact_mult;
    vector<double> contact_lowerto;
    vector<double> u;
    vector<double> fIp;
    vector<double> fIa;
    vector<double> fIs;
    vector<double> y;
    vector<double> rho;
    vector<double> tau;

    vector<double> seed_times;
    Discrete dist_seed_ages;

    Observer observer;
    vector<ScheduleEntry> schedule;

    string name;
    vector<string> group_names;

    void Set(Parameters* parent, string& name, Rcpp::RObject& value);

    void Recalculate()
    {
        if (needs_recalc)
        {
            if (matrices.empty())
                throw logic_error("No contact matrices defined.");
            if (matrices.size() != contact.size())
                throw logic_error("Number of contact components not equal to number of matrices.");

            unsigned int ncol = matrices[0].nc;
            unsigned int nrow = matrices[0].x.size() / ncol;

            auto c_mult = [&](int m) { if (!contact_mult.empty()) return contact_mult[m]; return 1.0; };
            auto c_lowerto = [&](int m) { if (!contact_lowerto.empty()) return contact_lowerto[m]; return std::numeric_limits<double>::max(); };

            cm = Matrix(0, nrow, ncol);

            for (unsigned int r = 0; r < nrow; ++r)
                for (unsigned int c = 0; c < ncol; ++c)
                    for (unsigned int m = 0; m < matrices.size(); ++m)
                        cm(r, c) += matrices[m](r, c) * min(contact[m] * c_mult(m), c_lowerto(m));

            needs_recalc = false;
        }
    }
};

struct Parameters
{
public:
    double time_step;
    double time0;
    double time1;
    unsigned int report_every;
    bool fast_multinomial;
    bool deterministic;

    vector<PopulationParameters> pop;

    vector<ProcessSpec> processes;    
    Matrix travel;
};

void PopulationParameters::Set(Parameters* parent, string& name, Rcpp::RObject& value)
{
    if (name == "contact" || name == "contact_mult" || name == "contact_lowerto" || name == "matrices")
        needs_recalc = true;

    if (false) {}
    _CheckSet(dE)
    _CheckSet(dIp)
    _CheckSet(dIa)
    _CheckSet(dIs)
    _CheckSet(dH)
    _CheckSet(dC)
    _CheckSet(size)
    _CheckSet(matrices)
    _CheckSet(contact)
    _CheckSet(contact_mult)
    _CheckSet(contact_lowerto)
    _CheckSet(u)
    _CheckSet(fIp)
    _CheckSet(fIa)
    _CheckSet(fIs)
    _CheckSet(y)
    _CheckSet(rho)
    _CheckSet(tau)
    _CheckSetParent(parent, travel)
    else
    {
        throw logic_error("Unrecognised parameter " + name + ".");
    }
}


// Helpers to set parameters
#define ParamAssign(t, v)               if (list.containsElementNamed(#v)) P.v = Rcpp::as<t>(list[#v]);
#define ParamMatrixAssign(v)            if (list.containsElementNamed(#v)) SetMatrix(P.v, Rcpp::as<Rcpp::NumericMatrix>(list[#v]));
#define ParamPopAssign(t, v, i)         if (popi.containsElementNamed(#v)) P.pop[i].v = Rcpp::as<t>(popi[#v]);
#define ParamPopMatrixAssign(v, i)      if (popi.containsElementNamed(#v)) SetMatrix(P.pop[i].v, Rcpp::as<Rcpp::NumericMatrix>(popi[#v]));

void SetMatrix(Matrix& mat, const Rcpp::NumericMatrix& rhs)
{
    mat = Matrix(0.0, rhs.nrow(), rhs.ncol());
    for (int r = 0; r < rhs.nrow(); ++r)
        for (int c = 0; c < rhs.ncol(); ++c)
            mat(r, c) = rhs(r, c);
}

void SetParameters(Parameters& P, Rcpp::List list, Randomizer& Rand)
{
    // TODO clean up and use namespace Rcpp
    // TODO use the ParamSet functions above for all of these parameter assignments
    ParamAssign(double, time_step);
    ParamAssign(double, time0);
    ParamAssign(double, time1);
    ParamAssign(unsigned int, report_every);
    ParamAssign(bool, fast_multinomial);
    ParamAssign(bool, deterministic);

    if (P.report_every != 1/P.time_step)
        throw("report_every must be the reciprocal of time_step.");

    if (list.containsElementNamed("pop"))
    {
        Rcpp::List populations = Rcpp::as<Rcpp::List>(list["pop"]);
        unsigned int np = populations.size();
        P.pop.assign(np, PopulationParameters());

        for (unsigned int i = 0; i < np; ++i)
        {
            Rcpp::List popi = Rcpp::as<Rcpp::List>(populations[i]);

            string type = Rcpp::as<string>(popi["type"]);
            if (type != "SEI3R")
                throw runtime_error("Population type must be SEI3R.");

            ParamPopAssign(vector<double>, dE, i);
            ParamPopAssign(vector<double>, dIp, i);
            ParamPopAssign(vector<double>, dIa, i);
            ParamPopAssign(vector<double>, dIs, i);
            ParamPopAssign(vector<double>, dH, i);
            ParamPopAssign(vector<double>, dC, i);

            if (P.fast_multinomial)
            {
                P.pop[i].dE.mn_approx.Set(P, Rand, P.pop[i].dE.weights);
                P.pop[i].dIp.mn_approx.Set(P, Rand, P.pop[i].dIp.weights);
                P.pop[i].dIa.mn_approx.Set(P, Rand, P.pop[i].dIa.weights);
                P.pop[i].dIs.mn_approx.Set(P, Rand, P.pop[i].dIs.weights);
                P.pop[i].dH.mn_approx.Set(P, Rand, P.pop[i].dH.weights);
                P.pop[i].dC.mn_approx.Set(P, Rand, P.pop[i].dC.weights);
            }

            ParamPopAssign(vector<double>, size, i);

            if (popi.containsElementNamed("matrices"))
            {
                Rcpp::List ml = Rcpp::as<Rcpp::List>(popi["matrices"]);
                for (unsigned int m = 0; m < ml.size(); ++m)
                {
                    Matrix mat;
                    SetMatrix(mat, Rcpp::as<Rcpp::NumericMatrix>(ml[m]));
                    P.pop[i].matrices.push_back(mat);
                }
            }    

            ParamPopAssign(vector<double>, contact, i);
            ParamPopAssign(vector<double>, contact_mult, i);
            ParamPopAssign(vector<double>, contact_lowerto, i);
            ParamPopAssign(vector<double>, u, i);
            ParamPopAssign(vector<double>, fIp, i);
            ParamPopAssign(vector<double>, fIa, i);
            ParamPopAssign(vector<double>, fIs, i);
            ParamPopAssign(vector<double>, y, i);
            ParamPopAssign(vector<double>, rho, i);
            ParamPopAssign(vector<double>, tau, i);

            ParamPopAssign(vector<double>, seed_times, i);
            ParamPopAssign(vector<double>, dist_seed_ages, i);

            // Read in observer
            if (popi.containsElementNamed("observer"))
            {
                Rcpp::RObject obs = Rcpp::as<Rcpp::RObject>(popi["observer"]);
                if (obs.isNULL())
                    P.pop[i].observer.null = true;
                else
                    P.pop[i].observer = Rcpp::as<Rcpp::Function>(obs);
            }

            // Read in scheduled changes
            Rcpp::List sched = Rcpp::as<Rcpp::List>(popi["schedule"]);
            for (unsigned int j = 0; j < sched.size(); ++j)
            {
                Rcpp::List se = Rcpp::as<Rcpp::List>(sched[j]);
                vector<ScheduleEntry> entries;
                double t = 0;

                for (unsigned int k = 0; k < se.size(); ++k)
                {
                    string name = Rcpp::as<string>(Rcpp::as<Rcpp::CharacterVector>(se.names())[k]);
                    Rcpp::RObject value = Rcpp::as<Rcpp::RObject>(se[k]);
                    if (name == "t")
                        t = Rcpp::as<double>(value);
                    else
                        entries.push_back(ScheduleEntry { -1, name, value });
                }

                for (unsigned int k = 0; k < entries.size(); ++k)
                {
                    entries[k].t = t;
                    P.pop[i].schedule.push_back(entries[k]);
                }
            }

            // Names
            ParamPopAssign(string, name, i);
            ParamPopAssign(vector<string>, group_names, i);

            // Calculate
            P.pop[i].Recalculate();
        }
    }

    // Read in processes
    if (list.containsElementNamed("processes"))
    {
        Rcpp::RObject r_processes = Rcpp::as<Rcpp::RObject>(list["processes"]);
        ParamSet(P.processes, r_processes);
    }
    
    ParamMatrixAssign(travel);
}
