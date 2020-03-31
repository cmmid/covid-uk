# covidm_run.R
# functions to run the model

# package dynamics returned by cm_backend_simulate into a single data.table
cm_package_dynamics = function(dynamics_list)
{
    dynamics = NULL;
    for (run in seq_along(dynamics_list))
    {
        setDT(dynamics_list[[run]]);
        dynamics = rbind(dynamics,
            cbind(run = run, dynamics_list[[run]])
        );
    }
    
    melt(dynamics, id.vars = c("run", "t", "population", "group"), variable.name = "compartment", value.name = "value");
}

# TODO -- in progress... interface likely to change.
cm_sample_parameters = function(fit, n)
{
    row = sample.int(nrow(fit$posterior), size = 1, replace = T);
    # Set parameters
    x = unlist(fit$posterior[row, ]);
    parameters = cm_translate_parameters(fit$parameters_func(fit$base_parameters, x));
    return (parameters)
}

# TODO -- in progress... interface likely to change.
cm_sample_fit = function(fit, n)
{
    dynamics = NULL;
    rows = sample.int(nrow(fit$posterior), size = n, replace = T);
    
    for (i in 1:n)
    {
        # Set parameters
        x = unlist(fit$posterior[rows[i],]);
        parameters = cm_translate_parameters(fit$parameters_func(fit$base_parameters, x));
    
        # Run simulation
        dyn = cm_backend_simulate(parameters, 1, fit$options$model_seed)$dynamics;
        dyn = cm_package_dynamics(dyn);
        dyn[, run := i];
        dyn = cm_annotate_dynamics(dyn, parameters);
        
        dynamics = rbind(dynamics, dyn);
    }
    
    return (dynamics)
}

# alter any duplicated names
cm_unduplicate_names = function(names)
{
    if (anyDuplicated(names)) {
        return (paste0(names, "__", 1:length(names)))
    }
    return (unlist(names))
}

# Annotate returned dynamics with more informative names for population and group columns.
cm_annotate_dynamics = function(dynamics, parameters)
{
    pnames = list()
    for (p in 1:length(parameters$pop)) {
        pnames[[p]] = parameters$pop[[p]]$name;
        if (is.null(pnames[[p]])) {
            pnames[[p]] = p;
        }
        group_names = parameters$pop[[p]]$group_names;
        if (!is.null(group_names)) {
            group_names = cm_unduplicate_names(group_names);
            dynamics[population == p, group_name := group_names[group]];
        } else {
            dynamics[population == p, group_name := as.character(group)];
        }
    }
    
    pnames = cm_unduplicate_names(pnames);
    dynamics[, population := pnames[population]];
    
    dynamics[, population := factor(population, levels = unique(population))];
    dynamics[, group := factor(group_name, levels = unique(group_name))];
    dynamics[, group_name := NULL];
    
    return (dynamics);
}

# TODO -- in progress... interface likely to change.
cm_simulate = function(parameters, n = 1, model_seed = 0)
{
    results = cm_backend_simulate(cm_translate_parameters(parameters), n, model_seed);
    dyn = cm_package_dynamics(results$dynamics);
    dyn = cm_annotate_dynamics(dyn, parameters);

    run = list(
        base_parameters = parameters,
        dynamics = dyn,
        csv = results$csv,
        changes = NULL, # some analogue of posterior here
        parameters_func = NULL, # some analogue of parameters_func here as well.
        options = list( # fill out
            model_seed = model_seed, # implement
            model_version = 1 # likewise
        )
    );
    
    class(run) = c("cm.run", "list");
    
    return (run)
}
