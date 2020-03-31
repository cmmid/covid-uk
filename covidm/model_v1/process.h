//
// PROCESS
//

enum SourceID
{
    srcS = 1000000,
    srcE,
    srcEp,
    srcEa,
    srcIp,
    srcIs,
    srcH,
    srcIa,
    srcI
};

const unsigned int Null = 999999;

struct ProcessSpec
{
    unsigned int source_id; // identifier for the compartment which, upon exiting, individuals enter this process
    string source_name;     // name of the above
    string type;            // ignored for now - multinomial or dirichlet multinomial

    vector<string> names;       // names of sub-processes
    vector<unsigned int> ids;   // identifiers of sub-process compartments
    vector<string> report;      // reporting mode of sub-processes: empty, "i", "o", "p", or a combination of these
    vector<vector<double>> prob;// probability by group of entering each sub-process from the source above; indexed by group then by subprocess
    vector<Discrete> delays;    // delays for each sub-process

    vector<unsigned int> p_cols;
    vector<unsigned int> p_ids;
    vector<unsigned int> i_cols;
    vector<unsigned int> i_ids;
    vector<unsigned int> o_cols;
    vector<unsigned int> o_ids;     // prevalence, incidence, and outcidence data columns and sub-process identifiers for this process
};
