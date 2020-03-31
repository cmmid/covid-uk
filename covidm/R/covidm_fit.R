# covidm_fit.R
# framework for model fitting

# Automatically set parameters if they start with a "."
cm_parameters_func_ = function(parameters, x)
{
    default = function(x, d) { if (is.na(x)) return(d) else return(as.numeric(x)) }
    
    direct = (substr(names(x), 1, 1) == ".");
    n = strsplit(names(x)[direct], ".", fixed = T);
    v = x[direct];
    
    for (i in seq_along(n)) {
        m = n[[i]];
        name = m[2];
        if (length(m) == 2) {
            j = 1;
            k0 = 1;
            k1 = length(parameters$pop[[j]][[name]]);
        } else if (length(m) == 3) {
            j = as.numeric(m[3]);
            k0 = 1;
            k1 = length(parameters$pop[[j]][[name]]);
        } else if (length(m) == 4) {
            j = as.numeric(m[3]);
            k0 = k1 = as.numeric(m[4]);
        } else {
            j = as.numeric(m[3]);
            k0 = as.numeric(m[4]);
            k1 = as.numeric(m[5]);
        }
        parameters$pop[[j]][[name]][k0:k1] = v[i];
    }
    
    return (parameters)
}


# Internal - likelihood function for MCMC backend
# x is a posterior sample, details contains extra stuff needed to evaluate the likelihood (see cm_fit).
cm_mcmc_likelihood_ = function(x, details, i)
{
    # Set parameters
    names(x) = names(details$priors);
    parameters = cm_parameters_func_(details$base_parameters, x);
    if (!is.null(details$parameters_func)) {
        parameters = details$parameters_func(parameters, x);
    }
    parameters = cm_translate_parameters(parameters);
    cm_check_parameters(parameters);
    
    # Run simulation 
    # TODO what to do about model seed?
    dynamics = cm_backend_simulate(parameters, details$fit_sampling, i)$dynamics;###details$model_seed);
    dynamics = cm_package_dynamics(dynamics);
    dynamics = cm_annotate_dynamics(dynamics, parameters);

    # Calculate and return likelihood
    if (details$fit_sampling == 1) {
        return (details$likelihood_func(parameters, dynamics, details$data, x))
    } else {
        dynamics = split(dynamics, by = "run");
        return (sum(sapply(dynamics, function(d) details$likelihood_func(parameters, d, details$data, x), simplify = T)) / details$fit_sampling);
    }
}

# Fit the model
cm_fit = function(base_parameters, priors, parameters_func, likelihood_func, data,
                  fit_method = "mcmc", fit_version = "latest", fit_seed = 0, fit_sampling = 1,
                  model_version = "latest", model_seed = 0,
                  mcmc_burn_in = 1000, mcmc_chains = 2 * length(priors), mcmc_samples = NA, mcmc_iter = 1000,
                  mcmc_init_opt = F, mcmc_reeval_likelihood = F,
                  opt_global = T, opt_global_algorithm = "ISRES", opt_global_maxeval = 1000, opt_global_ftol_abs = 1e-6,
                  opt_local = T, opt_local_algorithm = "NELDERMEAD", opt_local_maxeval = 1000, opt_local_ftol_abs = 1e-6,
                  verbose = T, options = list(
                      fit_method = fit_method, fit_version = fit_version, fit_seed = fit_seed, fit_sampling = fit_sampling,
                      model_version = model_version, model_seed = model_seed, 
                      mcmc_burn_in = mcmc_burn_in, mcmc_chains = mcmc_chains, 
                      mcmc_samples = mcmc_samples, mcmc_iter = mcmc_iter, 
                      mcmc_init_opt = mcmc_init_opt, mcmc_reeval_likelihood = mcmc_reeval_likelihood,
                      opt_global = opt_global, opt_global_algorithm = opt_global_algorithm, 
                      opt_global_maxeval = opt_global_maxeval, opt_global_ftol_abs = opt_global_ftol_abs, 
                      opt_local = opt_local, opt_local_algorithm = opt_local_algorithm, 
                      opt_local_maxeval = opt_local_maxeval, opt_local_ftol_abs = opt_local_ftol_abs,
                      verbose = verbose),
                  debug = F)
{
    # Check parameters
    # TODO implement these.
    if (options$fit_version != "latest")   { stop("fit_verstion must be \"latest\"."); }
    if (options$model_version != "latest") { stop("model_verstion must be \"latest\"."); }
    
    if (options$fit_sampling > 1 & base_parameters$deterministic) {
        message("cm_fit integrating over multiple runs per parameter combination, but running the deterministic model. Is this what you intended?");
    }

    # Create details list
    details = list(
        base_parameters = base_parameters,
        priors = priors,
        parameters_func = parameters_func,
        likelihood_func = likelihood_func,
        data = data,
        fit_method = options$fit_method,
        fit_version = options$fit_version,
        fit_seed = options$fit_seed,
        fit_sampling = options$fit_sampling,
        model_version = options$model_version,
        model_seed = options$model_seed
    );

    # Debug mode tests.
    if (debug) {
        cat("Debug mode.\n");
        cat("Sampling from prior.\n");
        psamp = cm_backend_prior_sample(priors);
        print(psamp);
        readline(prompt = "Press <Return> to run parameters_func . . .");
        
        parameters = details$parameters_func(details$base_parameters, psamp);
        ans = readline(prompt = "View results? (y/n) ");
        if (ans == "y") {
            View(parameters);
        }
        readline(prompt = "Press <Return> to run cm_translate_parameters . . .");
        
        parameters = cm_translate_parameters(parameters);
        ans = readline(prompt = "View results? (y/n) ");
        if (ans == "y") {
            View(parameters);
        }
        readline(prompt = "Press <Return> to run cm_check_parameters . . .");

        cm_check_parameters(parameters);
        readline(prompt = "Press <Return> to run simulation . . .");

        dynamics = cm_backend_simulate(parameters, 1, details$model_seed)$dynamics;
        dynamics = cm_package_dynamics(dynamics);
        dynamics = cm_annotate_dynamics(dynamics, parameters);

        ans = readline(prompt = "View results? (y/n) ");
        if (ans == "y") {
            View(dynamics);
        }
        readline(prompt = "Press <Return> to evaluate likelihood . . .");

        ll = details$likelihood_func(parameters, dynamics, details$data, psamp);
        print(ll);

        return (0);
    }

    # Run fitting
    message("Fitting model...");
    if (fit_method == "mcmc") {
        # Determine number of post-burn-in samples
        iter = options$mcmc_iter;
        if (!is.na(options$mcmc_samples)) {
            iter = ceiling(options$mcmc_samples / options$mcmc_chains);
        }
        
        # Optimize first?
        if (options$mcmc_init_opt) {
            # Get initial values for fitting from optimization
            initial = cm_backend_optimize(cm_mcmc_likelihood_, details, priors,
                options$fit_seed,
                options$opt_global, options$opt_global_algorithm, options$opt_global_maxeval, options$opt_global_ftol_abs,
                options$opt_local, options$opt_local_algorithm, options$opt_local_maxeval, options$opt_local_ftol_abs,
                options$verbose);
            setDT(initial);
            initial = unname(unlist(initial[, 5:ncol(initial)]));
            initial = matrix(initial, nrow = options$mcmc_chains, ncol = length(initial), byrow = T);
            
            # Perturb these slightly for the different chains
            perturbed = initial * matrix(rnorm(nrow(initial) * ncol(initial), 1.0, 1e-2), nrow = nrow(initial), ncol = ncol(initial));
            
            # Sample posterior with MCMC
            posterior = cm_backend_mcmc_init(cm_mcmc_likelihood_, details, priors, perturbed,
                options$fit_seed, options$mcmc_burn_in, options$mcmc_chains, iter, options$verbose, options$mcmc_reeval_likelihood, F, -1);
        } else {
            # Sample posterior with MCMC
            posterior = cm_backend_mcmc(cm_mcmc_likelihood_, details, priors,
                options$fit_seed, options$mcmc_burn_in, options$mcmc_chains, iter, options$verbose, options$mcmc_reeval_likelihood, F, -1);
        }
        setDT(posterior);
    } else if (fit_method == "optimize") {
        # Estimate posterior maximum using numerical optimization
        posterior = cm_backend_optimize(cm_mcmc_likelihood_, details, priors,
            options$fit_seed,
            options$opt_global, options$opt_global_algorithm, options$opt_global_maxeval, options$opt_global_ftol_abs,
            options$opt_local, options$opt_local_algorithm, options$opt_local_maxeval, options$opt_local_ftol_abs,
            options$verbose);
        setDT(posterior);
    } else {
        stop("fit_method must be \"mcmc\" or \"optimize\".");
    }
    
    # Prepare results
    results = list(
        base_parameters = base_parameters,
        priors = priors,
        parameters_func = parameters_func,
        likelihood_func = likelihood_func,
        data = data,
        posterior = posterior,
        options = options
    );
    class(results) = c("cm.fit", "list");

    return (results)
}

# Return the DIC for an MCMC model fit.
cm_DIC = function(fit)
{
    # Ensure this fit is of the right type
    if (fit$options$fit_method != "mcmc") {
        stop("cm_DIC needs a model fitted with MCMC to calculate DIC.");
    }
    
    # Get posterior mean
    Ep = unlist(colMeans(fit$posterior[, 5:ncol(fit$posterior)]));
    
    # Run simulation with posterior mean
    parameters = cm_translate_parameters(fit$parameters_func(fit$base_parameters, Ep));
    run = cm_simulate(parameters)$dynamics; # TODO ADD SEED!!!...?

    # Calculate DIC
    D_mean = -2 * fit$likelihood_func(parameters, run$dynamics, fit$data, Ep);
    mean_D = fit$posterior[, mean(-2 * lp)];
    return (c(DIC = 2 * mean_D - D_mean));
}

# Return the AIC for an MCMC model fit.
cm_AIC = function(fit)
{
    # Get estimated maximum likelihood
    L = fit$posterior[, max(ll)];

    # Calculate AIC
    return (c(AIC = -2 * L + 2 * (ncol(fit$posterior) - 4)));
}

