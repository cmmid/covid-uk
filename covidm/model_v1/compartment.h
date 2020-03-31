//
// COMPARTMENT
//

// A compartment containing individuals
class Compartment
{
public:
    // Construct the compartment
    Compartment()
     : size(0) { }

    // Add n individuals to the compartment, with maturation times t
    // Note that t[0] means individuals who will be in the compartment for 0 time steps
    void Add(Parameters& P, Randomizer& Rand, double n, Discrete& mat)
    {
        // Expand compartment span if needed
        if (contents.size() < mat.weights.size())
            contents.resize(mat.weights.size(), 0);

        // Seed compartment
        size += n;
        if (P.deterministic)
        {
            for (unsigned int i = 0; i < mat.weights.size(); ++i)
                contents[i] += mat.weights[i] * n;
        }
        else
        {
            if (P.fast_multinomial)
                mat.mn_approx(n, mat.storage);
            else
                Rand.Multinomial(n, mat.weights, mat.storage);
            for (unsigned int i = 0; i < mat.weights.size(); ++i)
                contents[i] += mat.storage[i];
        }
    }

    // Mature the compartment by one time step, returning the number of individuals who have left the compartment.
    double Mature()
    {
        if (contents.empty())
            return 0.0;

        auto m = contents.front();
        size -= m;
        contents.erase(contents.begin());
        contents.push_back(0);
        return m;
    }

    // Return the total number of individuals in the compartment
    double Size()
    {
        return size;
    }
    
private:
    // Individuals in the compartment indexed by the number of days to maturation
    vector<double> contents;

    // Total number of individuals in the compartment
    double size;
};