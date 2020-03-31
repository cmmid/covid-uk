// corona.cpp
// SEI3HR dynamics

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppGSL)]]

#include <vector>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <ctime>
#include <limits>
#include "../randomizer/randomizer.h"
#include "../randomizer/distribution.h"
using namespace std;
#include "helper.h"
#include "process.h"
#include "parameters.h"
#include "compartment.h"
#include "reporter.h"

bool Observer::operator()(Parameters* parent, PopulationParameters& pp, double t, Reporter& rep)
{
    if (null)
        return true;
    
    using namespace Rcpp;
    RObject observer_return = func(t, rep.dynamics_df);

    if (is<List>(observer_return))
    {
        List ret = as<List>(observer_return);
        if (ret.containsElementNamed("changes"))
        {
            if (is<List>(ret["changes"]))
            {
                List changes = as<List>(ret["changes"]);
                for (unsigned int i = 0; i < changes.size(); ++i)
                {
                    string name = as<string>(as<CharacterVector>(changes.names())[i]);
                    RObject value = as<RObject>(changes[i]);
                    pp.Set(parent, name, value);
                }
                pp.Recalculate();
            }
            else
            {
                throw logic_error("changes must be a list with named entries.");
            }
        }
        if (ret.containsElementNamed("print"))
        {
            string message = as<string>(ret["print"]);
            cout << message << '\n';
        }
        if (ret.containsElementNamed("csv"))
        {
            string csv = as<string>(ret["csv"]);
            rep.csv += csv;
            if (csv.back() != '\n')
                rep.csv += '\n';
        }
        if (ret.containsElementNamed("halt"))
        {
            bool halt = as<bool>(ret["halt"]);
            return !halt;
        }
    }
    else if (!observer_return.isNULL())
    {
        throw runtime_error("Observer function must return either a list or NULL.");
    }
    return true;
}


//
// MODEL DYNAMICS
//

// A population of individuals, with SEI3HR dynamics.
class Population
{
public:
    // Construct a population with the specified size by age group; initially all uninfected
    Population(Parameters& P, unsigned int pindex)
     : seed_row(0), schedule_row(0), p(pindex)
    {
        // Set up built-in compartments
        S  = P.pop[p].size;
        E  = vector<Compartment>(S.size());
        Ip = vector<Compartment>(S.size());
        Ia = vector<Compartment>(S.size());
        Is = vector<Compartment>(S.size());
        H  = vector<Compartment>(S.size());
        C  = vector<Compartment>(S.size());
        R  = vector<double>(S.size(), 0.);

        // Set up user-specified processes
        unsigned int n_pc = 0;
        for (auto& p : P.processes)
            n_pc += p.ids.size();
        pc = vector<vector<Compartment>>(n_pc, vector<Compartment>(S.size()));
        pci = vector<double>(pc.size(), 0.);
        pco = vector<double>(pc.size(), 0.);
    }

    // Do seeding and calculate contagiousness
    void Contagiousness(Parameters& P, Randomizer& Rand, double t, vector<double>& contag)
    {
        auto add = [&](unsigned int age, double n)
        {
            if (S[age] < n)
                throw logic_error("Not enough unexposed individuals to seed new infections.");
            S[age] -= n;
            E[age].Add(P, Rand, n, P.pop[p].dE);
        };

        // Do seeding
        while (seed_row < P.pop[p].seed_times.size() && t >= P.pop[p].seed_times[seed_row])
        {
            if (P.deterministic)
            {
                for (unsigned int a = 0; a < S.size(); ++a)
                    add(a, P.pop[p].dist_seed_ages.weights[a]);
            }
            else
            {
                Rand.Multinomial(1, P.pop[p].dist_seed_ages.weights, P.pop[p].dist_seed_ages.storage);
                for (unsigned int a = 0; a < S.size(); ++a)
                {
                    if (P.pop[p].dist_seed_ages.storage[a] == 1)
                    {
                        add(a, 1);
                        break;
                    }
                }
            }
            ++seed_row;
        }

        // Do scheduled changes
        while (schedule_row < P.pop[p].schedule.size() && t >= P.pop[p].schedule[schedule_row].t)
        {
            P.pop[p].Set(&P, P.pop[p].schedule[schedule_row].variable, P.pop[p].schedule[schedule_row].value);
            ++schedule_row;
        }
        P.pop[p].Recalculate();

        // Calculate contagiousness from this population
        for (unsigned int a = 0; a < contag.size(); ++a)
            contag[a] = (P.pop[p].size[a] == 0) ? 0 : (P.pop[p].fIp[a] * Ip[a].Size() + P.pop[p].fIa[a] * Ia[a].Size() + P.pop[p].fIs[a] * Is[a].Size()) / P.pop[p].size[a];
    }

    // Execute one time step's events
    bool Tick(Parameters& P, Randomizer& Rand, double t, vector<double>& infec, Reporter& rep)
    {
        (void) t;

        // Calculate force of infection in this compartment
        lambda.assign(infec.size(), 0.0);
        for (unsigned int a = 0; a < lambda.size(); ++a)
            for (unsigned int b = 0; b < lambda.size(); ++b)
                lambda[a] += P.pop[p].u[a] * P.pop[p].cm(a,b) * infec[b];

        // Helpers
        auto multinomial = [&](double n, vector<double>& p, vector<double>& nd_out, vector<unsigned int>& ni_out) {
            nd_out.resize(p.size(), 0.);
            if (P.deterministic)
            {
                for (unsigned int i = 0; i < p.size(); ++i)
                    nd_out[i] = n * p[i];
            }
            else
            {
                ni_out.resize(p.size(), 0);
                Rand.Multinomial(n, p, ni_out);
                for (unsigned int i = 0; i < p.size(); ++i)
                    nd_out[i] = ni_out[i];
            }
        };

        auto binomial = [&](double n, double p) {
            if (P.deterministic)
                return n * p;
            else
                return (double)Rand.Binomial(n, p);
        };

        // Do state transitions and reporting for each age group
        for (unsigned int a = 0; a < lambda.size(); ++a)
        {
            // 0. Report prevalences
            if (t == (int)t)
            {
                // Built-in states
                rep(t, p, a, 0) = S[a];
                rep(t, p, a, 1) = E[a].Size();
                rep(t, p, a, 2) = Ip[a].Size();
                rep(t, p, a, 3) = Is[a].Size();
                rep(t, p, a, 4) = Ia[a].Size();
                rep(t, p, a, 5) = R[a];

                // User-specified processes
                for (auto& process : P.processes)
                {
                    for (unsigned int i = 0; i < process.p_cols.size(); ++i)
                        rep(t, p, a, process.p_cols[i]) = pc[process.p_ids[i]][a].Size();
                }

            }

            // 1. Built-in states

            // S -> E
            double nS_E = binomial(S[a], 1.0 - exp(-lambda[a] * P.time_step));
            S[a] -= nS_E;
            E[a].Add(P, Rand, nS_E, P.pop[p].dE);

            // E -> Ip/Ia
            double nE_Ipa = E[a].Mature();
            double nE_Ip = binomial(nE_Ipa, P.pop[p].y[a]);
            double nE_Ia = nE_Ipa - nE_Ip;
            Ip[a].Add(P, Rand, nE_Ip, P.pop[p].dIp);
            Ia[a].Add(P, Rand, nE_Ia, P.pop[p].dIa);

            // Ip -> Is -- also, true case onsets
            double nIp_Is = Ip[a].Mature();
            Is[a].Add(P, Rand, nIp_Is, P.pop[p].dIs);

            // Reported cases
            double n_to_report = binomial(nIp_Is, P.pop[p].rho[a]);
            C[a].Add(P, Rand, n_to_report, P.pop[p].dC);
            double n_reported = C[a].Mature();

            // Is -> H
            double nIs_H = Is[a].Mature();
            H[a].Add(P, Rand, nIs_H, P.pop[p].dH);

            // H -> R
            double nH_R = H[a].Mature();
            R[a] += nH_R;

            // Ia -> R
            double nIa_R = Ia[a].Mature();
            R[a] += nIa_R;

            // 2. User-specified processes

            // Mature process compartments
            for (auto& process : P.processes)
                for (auto compartment_id : process.ids)
                    if (compartment_id != Null)
                        pco[compartment_id] = pc[compartment_id][a].Mature();

            // Seed processes
            for (auto& process : P.processes)
            {
                // Determine number of individuals entering the process
                double n_entering = 0.;
                switch (process.source_id)
                {
                    case srcS:
                        n_entering = nS_E; break;
                    case srcE:
                        n_entering = nE_Ipa; break;
                    case srcEp:
                        n_entering = nE_Ip; break;
                    case srcEa:
                        n_entering = nE_Ia; break;
                    case srcIp:
                        n_entering = nIp_Is; break;
                    case srcIs:
                        n_entering = nIs_H; break;
                    case srcH:
                        n_entering = nH_R; break;
                    case srcIa:
                        n_entering = nIa_R; break;
                    case srcI:
                        n_entering = nH_R + nIa_R; break;
                    default:
                        n_entering = pco[process.source_id]; break;
                }

                multinomial(n_entering, process.prob[a], nd_out, ni_out);

                // Seed this process's compartments
                unsigned int c = 0;
                for (unsigned int compartment_id : process.ids)
                {
                    if (compartment_id != Null)
                    {
                        pc[compartment_id][a].Add(P, Rand, nd_out[c], process.delays[c]);
                        pci[compartment_id] = nd_out[c];
                    }
                    ++c;
                }
            }

            // 3. Report incidence / outcidence
            // Built-in states
            rep(t, p, a, 6) += nIp_Is;
            rep(t, p, a, 7) += n_reported;
            rep(t, p, a, 8) += nE_Ia;

            // User-specified processes
            for (auto& process : P.processes)
            {
                for (unsigned int i = 0; i < process.i_cols.size(); ++i)
                    rep(t, p, a, process.i_cols[i]) += pci[process.i_ids[i]];
                for (unsigned int i = 0; i < process.o_cols.size(); ++i)
                    rep(t, p, a, process.o_cols[i]) += pco[process.o_ids[i]];
            }

        }

        // Run observer at the last time step of each day.
        if (t + P.time_step == int(t + P.time_step))
            return P.pop[p].observer(&P, P.pop[p], (int)t, rep);
        else
            return true;
    }

//private:
    vector<double> lambda;
    vector<double> S, R;                        // Susceptible, recovered
    vector<Compartment> E, Ip, Ia, Is, H, C;    // Exposed, presymptomatic, asymptomatic, symptomatic, hospitalised, cases (reported)
    unsigned int seed_row;                      // Which seed event is next
    unsigned int schedule_row;                  // Which schedule event is next
    unsigned int p;                             // Which population this is
    vector<vector<Compartment>> pc;             // User-specified process compartments, indexed by process id, then group
    vector<unsigned int> ni_out;                // Temporary storage
    vector<double> nd_out;                      // Temporary storage
    vector<double> pci;                         // Temporary storage
    vector<double> pco;                         // Temporary storage
};

// A metapopulation, containing multiple subpopulations.
class Metapopulation
{
public:
    Metapopulation(Parameters& P)
    {
        for (unsigned int i = 0; i < P.pop.size(); ++i)
            pops.push_back(Population(P, i));
    }

    // Execute one time step's events
    bool Tick(Parameters& P, Randomizer& Rand, double t, unsigned int ts, Reporter& rep)
    {
        ///Report(P, t, ts, report, false);

        unsigned int n_ages = P.pop[0].size.size();

        // Calculate contagiousness from each population
        // NOTE -- 'contag' subscripted first by j, then by a.
        // It's the effective number of infectious individuals FROM subpop j of age a.
        contag.assign(pops.size(), vector<double>(n_ages, 0.0));
        for (unsigned int j = 0; j < pops.size(); ++j)
            pops[j].Contagiousness(P, Rand, t, contag[j]);

        // note -- 'infec' subscripted first by i, then by a
        // It's the effective number of infectious individuals who are CURRENTLY IN subpop i of age a.
        infec.assign(pops.size(), vector<double>(n_ages, 0.0)); 
        for (unsigned int i = 0; i < pops.size(); ++i)
            for (unsigned int j = 0; j < pops.size(); ++j)
                for (unsigned int a = 0; a < n_ages; ++a)
                    infec[i][a] += P.travel(j, i) * contag[j][a] * (j != i ? P.pop[j].tau[a] : 1.0);

        // Update populations
        bool keep_going = true;
        for (unsigned int i = 0; i < pops.size(); ++i)
            keep_going = keep_going && pops[i].Tick(P, Rand, t, infec[i], rep);

        return keep_going;
    }

//private:
    vector<vector<double>> contag;
    vector<vector<double>> infec;
    vector<Population> pops;
};

Reporter RunSimulation(Parameters& P, Randomizer& Rand)
{
    Metapopulation mp(P);
    Reporter rep(P);

    // Run simulation
    unsigned int time_steps = (1 + P.time1 - P.time0) / P.time_step;
    bool ended_early = false;
    for (unsigned int ts = 0; ts < time_steps; ++ts)
    {
        if (!mp.Tick(P, Rand, P.time0 + ts * P.time_step, ts, rep))
        {
            ended_early = true;
            break;
        }
    }

    return rep;
}

#ifndef CORONA_CPP
#include "Rcpp_interface.h"
#endif

int main()
{
    return 0;
}