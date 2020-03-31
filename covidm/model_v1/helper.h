#include <gsl/gsl_sf_gamma.h>

//
// TIMING
//

// Clock functions
static vector<double> ClockTimes;
static double C0;

double Clock()  { return double(clock()) / CLOCKS_PER_SEC; }

void StartClocking()
{
    ClockTimes.clear();
    C0 = Clock();
}

void ClockCheckpoint(unsigned int cp)
{
    double C1 = Clock();
    if (cp >= ClockTimes.size()) ClockTimes.resize(cp + 1, 0.0);
    ClockTimes[cp] += C1 - C0;
    C0 = C1;
}

void ShowClockInfo()
{
    double total = 0;
    for (unsigned int i = 0; i < ClockTimes.size(); ++i)
    {
        cout << "Checkpoint " << i << ": " << ClockTimes[i] << "\n";
        total += ClockTimes[i];
    }
    cout << "Total time: " << total << " seconds.\n";
}

//
// DISCRETE DISTRIBUTION
//

struct Parameters;

class MNApprox
{
public:
    static const unsigned int NVariants = 32;
    static const unsigned int VariantMask = 31;

    void Set(Parameters& P, Randomizer& Rand, vector<double>& p);

    void operator()(unsigned int N, vector<unsigned int>& out)
    {
        out.assign(out.size(), 0);
        for (unsigned int shift = 0; shift < 32; ++shift)
        {
            if ((N >> shift) & 1)
            {
                for (unsigned int i = 0; i < out.size(); ++i)
                {
                    out[i] += x[shift][cycle[shift]][i];
                }
                cycle[shift] = (cycle[shift] + 1) & VariantMask;
            }
        }
    }

private:
    vector<vector<vector<unsigned int>>> x;
    vector<unsigned int> cycle;
};


// An arbitrary discrete distribution, used for maturation times in compartments
class Discrete
{
public:
    // Create distribution from set of unnormalised weights
    void operator=(std::vector<double> uw)
    {
        weights = uw;

        double total = accumulate(weights.begin(), weights.end(), 0.0);
        for (auto& w : weights)
            w /= total;

        storage.assign(weights.size(), 0);
    }

    // Normalised weights
    std::vector<double> weights;

    // Storage for draws of integers
    std::vector<unsigned int> storage;

    // Fast approximate multinomial draws
    MNApprox mn_approx;
};


//
// MATRIX
//

struct Matrix
{
    Matrix()
     : nc(0)
    { }

    Matrix(double X, unsigned int nrow, unsigned int ncol)
     : x(nrow * ncol, X), nc(ncol)
    { }

    Matrix(vector<double> X, unsigned int nrow, unsigned int ncol)
     : x(X), nc(ncol)
    {
        if (x.size() != nrow * ncol)
            throw runtime_error("Improperly sized matrix.");
    }

    double& operator()(unsigned int i, unsigned int j)
    {
        return x[i * nc + j];
    }

    unsigned int NCol() { return nc; }
    unsigned int NRow() { return x.size() / nc; }

    vector<double> x;
    unsigned int nc;
};


void MNApprox::Set(Parameters& P, Randomizer& Rand, vector<double>& p)
{
    x.assign(32, vector<vector<unsigned int>>(NVariants, vector<unsigned int>(p.size(), 0)));
    cycle.assign(NVariants, 0);
    for (unsigned int shift = 0; shift < 32; ++shift)
        for (unsigned int variant = 0; variant < NVariants; ++variant)
            gsl_ran_multinomial(Rand.GSL_RNG(), p.size(), 1 << shift, &p[0], &x[shift][variant][0]);
}
