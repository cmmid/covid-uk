# covidm_params.R
# parameters for the model

# return translated parameters to work with the backend,
# i.e. fix any times expressed in dates to be expressed in days since date0.
cm_translate_parameters = function(p)
{
    translate_time = function(t) {
        if (is.numeric(t)) {
            return (t)
        } else {
            return (as.numeric(ymd(t) - ymd(p$date0)));
        }
    }
    
    p$time0 = translate_time(p$time0);
    p$time1 = translate_time(p$time1);
    
    for (pi in seq_along(p$pop)) {
        p$pop[[pi]]$seed_times = sapply(p$pop[[pi]]$seed_times, translate_time);
        for (si in seq_along(p$pop[[pi]]$schedule)) {
            p$pop[[pi]]$schedule[[si]]$t = translate_time(p$pop[[pi]]$schedule[[si]]$t);
        }
    }
    
    return (p);
}

# check parameters for validity.
cm_check_parameters = function(parameters)
{
    req = function(v, x)
    {
        if (!exists(v, parameters)) {
            stop(paste0("Parameter ", v, " required, but not found."));
        } else if (!eval(parse(text = x), parameters)) {
            stop(paste0("Parameters check: ", x, " failed."));
        }
    }
    
    reqp = function(i, v, x)
    {
        if (!exists(v, parameters$pop[[i]])) {
            stop(paste0("Population parameter ", v, " required, but not found in population ", i, "."));
        } else if (!eval(parse(text = x), parameters$pop[[i]])) {
            stop(paste0("Parameters check: ", x, " failed in population ", i, "."));
        }
    }

    req("time_step",        "is.numeric(time_step) & time_step > 0");
    req("date0",            "is.Date(ymd(date0))");
    req("time0",            "is.numeric(time0)");
    req("time1",            "is.numeric(time1)");
    req("report_every",     "is.numeric(report_every) & report_every == 1. / time_step");
    req("fast_multinomial", "is.logical(fast_multinomial)");
    req("deterministic",    "is.logical(deterministic)");
    req("pop",              "is.list(pop) & length(pop) > 0");
    req("travel",           "is.matrix(travel) & all(travel >= 0) & nrow(travel) == length(parameters$pop) & ncol(travel) == nrow(travel)");
    
    for (i in 1:length(parameters$pop))
    {
        reqp(i, "type", "type == 'SEI3R'");
        reqp(i, "dE",   "is.numeric(dE) & all(dE >= 0) & any(dE > 0)");
        reqp(i, "dIp",  "is.numeric(dIp) & all(dIp >= 0) & any(dIp > 0)");
        reqp(i, "dIs",  "is.numeric(dIs) & all(dIs >= 0) & any(dIs > 0)");
        reqp(i, "dIa",  "is.numeric(dIa) & all(dIa >= 0) & any(dIa > 0)");
        reqp(i, "dH",   "is.numeric(dH) & all(dH >= 0) & any(dH > 0)");
        reqp(i, "dC",   "is.numeric(dC) & all(dC >= 0) & any(dC > 0)");
        reqp(i, "size", "is.numeric(size) & all(size >= 0) & any(size > 0)");
        reqp(i, "matrices", "is.list(matrices) & length(matrices) > 0");
        for (m in 1:length(parameters$pop[[i]]$matrices))
        {
            reqp(i, "matrices", sprintf("is.matrix(matrices[[%d]]) & all(matrices[[%d]] >= 0) & nrow(matrices[[%d]]) == length(size) & ncol(matrices[[%d]]) == nrow(matrices[[%d]])", m, m, m, m, m));
        }
        reqp(i, "contact",         "is.numeric(contact) & length(contact) == length(matrices) & all(contact >= 0)");
        reqp(i, "contact_mult",    "(is.numeric(contact_mult) & length(contact_mult) == length(matrices) & all(contact_mult >= 0)) | length(contact_mult) == 0");
        reqp(i, "contact_lowerto", "(is.numeric(contact_lowerto) & length(contact_lowerto) == length(matrices) & all(contact_lowerto >= 0)) | length(contact_lowerto) == 0");
        
        reqp(i, "u",       "is.numeric(u) & length(u) == length(size) & all(u >= 0)");
        reqp(i, "y",       "is.numeric(y) & length(y) == length(size) & all(y >= 0)");
        reqp(i, "fIp",     "is.numeric(fIp) & length(fIp) == length(size) & all(fIp >= 0)");
        reqp(i, "fIa",     "is.numeric(fIa) & length(fIa) == length(size) & all(fIa >= 0)");
        reqp(i, "fIs",     "is.numeric(fIs) & length(fIs) == length(size) & all(fIs >= 0)");
        reqp(i, "rho",     "is.numeric(rho) & length(rho) == length(size) & all(rho >= 0)");
        reqp(i, "tau",     "is.numeric(tau) & length(tau) == length(size) & all(tau >= 0)");
        
        reqp(i, "seed_times",      "is.numeric(seed_times) & !is.unsorted(seed_times)");
        reqp(i, "dist_seed_ages",  "is.numeric(dist_seed_ages) & length(dist_seed_ages) == length(size)");
        
        if (!is.null(parameters$pop[[i]]$observer)) {
            if (!(is.function(parameters$pop[[i]]$observer) & length(formals(parameters$pop[[i]]$observer) == 4))) {
                stop(paste0("observer has to be either NULL or a function taking 4 arguments, but is not in population", i));
            }
        }
        reqp(i, "schedule", "is.list(schedule)");
        schedule_times = sapply(parameters$pop[[i]]$schedule, function(x) x$t);
        if (is.unsorted(schedule_times)) {
            stop(paste0("elements t of schedule_times need to be ordered, but are not in population ", i));
        }
    }
}

# Get demographics for a given location, with error checking.
cm_get_demographics = function(dem_location, n_groups = NULL)
{
    # Get demographics.
    demographics = cm_populations[name == dem_location];
    if (nrow(demographics) == 0) {
        message(paste0("Could not find demographics for dem_location ", dem_location, "."));
        answer = readline(prompt = "View options? (y/n) ");
        if (answer != "n") {
            print(cm_populations[, unique(name)]);
        }
        stop();
    }
    if (nrow(demographics) != demographics[, uniqueN(age)]) {
        stop(paste0("Age not unique in cm_population[name == \"", dem_location, "\"]. This means cm_populations is misspecified for this location."));
    }
    
    # Adjust number of age groups if needed.
    if (!is.null(n_groups)) {
        if (n_groups > nrow(demographics)) {
            stop(sprintf("Requested %d age groups for model (up to %d+), but demographic data only goes up to %s.",
                         n_groups, (n_groups - 1) * 5, demographics[.N, age]));
        } else if (n_groups < nrow(demographics)) {
            demographics[n_groups]$f = demographics[n_groups:.N, sum(f)];
            demographics[n_groups]$m = demographics[n_groups:.N, sum(m)];
            demographics[n_groups]$age = demographics[n_groups, sub("-[0-9]+", "\\+", age)];
            demographics = demographics[1:n_groups];
        }
    }
    return (demographics);
}

# Get matrices for a given location, with error checking.
cm_get_matrices = function(mat_location, dem_location = NULL)
{
    # If requested, guess mat_location from dem_location.
    guess = F;
    if (mat_location == "guess") {
        if (is.null(dem_location)) {
            stop("cm_get_matrices needs a dem_location to guess the matrices location from.");
        }
        mat_location = dem_location;
        guess = T;
    }
    
    # Try to find matrices for mat_location.
    if (!mat_location %in% names(cm_matrices)) {
        message(paste0("Could not find matrices for mat_location ", mat_location, "."));
        answer = readline(prompt = "View options? (y/n) ");
        if (answer != "n") {
            print(names(cm_matrices));
        }
        stop();
    }
    if (sum(mat_location %in% names(cm_matrices) > 1)) {
        stop(paste0("Duplicate entries for ", mat_location, " in cm_matrices. This means cm_matrices is misspecified."));
    }
    mat = cm_matrices[[mat_location]];
    if (is.list(mat) & length(mat) > 0 & all(sapply(mat, is.matrix))) {
        return (mat);
    } else {
        if (guess) {
            stop(paste0("Could not guess mat_location for dem_location ", dem_location, "."));
        }
        stop(paste0("No valid entry in cm_matrices for matrix location ", mat_location, "."));
    }
}

# Split matrices
# ex_in: bounds separate into contacts *ex*clusively between lower groups, and *in*clusively between higher groups.
cm_split_matrices_ex_in = function(parameters, bounds)
{
    for (pi in seq_along(parameters$pop))
    {
        ng = nrow(parameters$pop[[pi]]$matrices[[1]]);
        if (any(bounds < 1 | bounds > ng)) {
            stop("Bounds must lie within [1, nrow(mat)] for splitting contact matrices.");
        }
        nmat0 = length(parameters$pop[[pi]]$matrices);
        parameters$pop[[pi]]$matrices = rep(parameters$pop[[pi]]$matrices, length(bounds) + 1);
        
        for (b in seq_along(bounds))
        {
            lb = floor(bounds[b]);
            fb = bounds[b] %% 1;
            mask1 = matrix(1, nrow = ng, ncol = ng);
            if (lb > 1) {
                mask1[1:(lb - 1), 1:(lb - 1)] = 0;
            }
            mask1[lb, 1:lb] = 1 - fb;
            mask1[1:lb, lb] = 1 - fb;
            mask0 = 1 - mask1;
            
            for (m in seq_len(nmat0)) {
                names(parameters$pop[[pi]]$matrices)[m + b * nmat0] = paste0(names(parameters$pop[[pi]]$matrices)[m + b * nmat0], b + 1);
                parameters$pop[[pi]]$matrices[[m + (b - 1) * nmat0]] = mask0 * parameters$pop[[pi]]$matrices[[m + (b - 1) * nmat0]];
                parameters$pop[[pi]]$matrices[[m +       b * nmat0]] = mask1 * parameters$pop[[pi]]$matrices[[m +       b * nmat0]];
            }
        }
        parameters$pop[[pi]]$contact = rep_len(parameters$pop[[pi]]$contact, nmat0 * (length(bounds) + 1));
    }
    
    return (parameters)
}

# TODO
# cm_split_matrices_in_ex
# cm_split_matrices_custom

# Get default population parameters, SEI3R model
cm_base_pop_SEI3R = function(n_groups)
{
    list(
        type = "SEI3R",
        dE  = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p, # Derived from Backer et al Eurosurveillance
        dIp = cm_delay_gamma(2.4, 4.0, t_max = 60, t_step = 0.25)$p, # Derived from Backer et al Eurosurveillance
        dIa = cm_delay_gamma(7.0, 4.0, t_max = 60, t_step = 0.25)$p, # Assumed 7 days subclinical shedding
        dIs = cm_delay_gamma(3.2, 3.7, t_max = 60, t_step = 0.25)$p, # Zhang et al 2020
        dH = 1, # hospitalization ignored
        dC = 1, # no case reporting delay
        
        size = rep(1000, n_groups),
        matrices = list(base = diag(n_groups) * 0.5 + 0.5/n_groups),
        contact = 1,
        contact_mult = numeric(),
        contact_lowerto = numeric(),
        u = rep(0.08, n_groups),
        y = rep(0.5, n_groups),
        fIp = rep(1, n_groups),
        fIs = rep(1, n_groups),
        fIa = rep(0.5, n_groups),
        rho = rep(1, n_groups),
        tau = rep(1, n_groups),
        
        seed_times = 1,
        dist_seed_ages = rep(1, n_groups),
        
        schedule = list(),
        observer = NULL
    )
}

# Build parameters for a single location, SEI3R model
cm_build_pop_SEI3R = function(dem_location, mat_location = "guess",
    dE = NULL, dIp = NULL, dIs = NULL, dIa = NULL, dH = NULL, dC = NULL, contact = NULL,
    u = NULL, y = NULL, fIp = NULL, fIa = NULL, fIs = NULL, rho = NULL, tau = NULL,
    seed_times = NULL, dist_seed_ages = NULL, observer = NULL, schedule = NULL)
{
    # Get desired demographics and matrices.
    matrices = cm_get_matrices(mat_location, dem_location);
    n_groups = nrow(matrices[[1]]);
    demographics = cm_get_demographics(dem_location, n_groups);
    
    # Get base population parameters.
    pop = cm_base_pop_SEI3R(n_groups);
    
    # Set population parameters.
    assign = function(pop, name, value) {
        if (!is.null(value)) {
            pop[[name]] = value;
        }
        return (pop);
    }
    
    assign_g = function(pop, name, value, n_groups) {
        if (!is.null(value)) {
            if (length(value) == 1 | length(value) == n_groups) {
                pop[[name]] = rep_len(value, n_groups)
            }
        }
        return (pop);
    }
    
    pop$name = dem_location;
    pop$group_names = colnames(matrices[[1]]);
    
    pop = assign(pop, "dE", dE);
    pop = assign(pop, "dIp", dIp);
    pop = assign(pop, "dIs", dIs);
    pop = assign(pop, "dIa", dIa);
    pop = assign(pop, "dH", dH);
    pop = assign(pop, "dC", dC);
    pop$size = demographics[, round((f + m) * 1000)];
    pop$matrices = matrices;
    pop$contact = rep(1, length(matrices));
    pop = assign(pop, "contact", contact);
    
    pop = assign_g(pop, "u", u, n_groups);
    pop = assign_g(pop, "y", y, n_groups);
    pop = assign_g(pop, "fIp", fIp, n_groups);
    pop = assign_g(pop, "fIs", fIs, n_groups);
    pop = assign_g(pop, "fIa", fIa, n_groups);
    pop = assign_g(pop, "rho", rho, n_groups);
    pop = assign_g(pop, "tau", tau, n_groups);

    pop = assign(pop, "seed_times", seed_times);
    pop = assign(pop, "dist_seed_ages", dist_seed_ages);
    pop = assign(pop, "observer", observer);
    pop = assign(pop, "schedule", schedule);
    
    return (pop)
}

# Get default simulation parameters, SEI3R model
cm_base_parameters_SEI3R = function(n_groups = 1, pop = cm_base_pop_SEI3R(n_groups))
{
    # If just a single population, rather than a list of populations, has been passed to this function, rectify that.
    if (is.character(pop$type)) {
        pop = list(pop);
    }
    
    list(
        time_step = 0.25,
        date0 = "2020-01-01",
        time0 = 0,
        time1 = 365,
        report_every = 4,
        fast_multinomial = F,
        deterministic = T,
        pop = pop,
        travel = diag(length(pop)),
        processes = NULL
    )
}

# Build parameters for one or several locations, SEI3R model
cm_parameters_SEI3R = function(dem_locations, mat_locations = "guess", date_start = "2020-03-01", date_end = "2021-03-01", deterministic = T, processes = NULL,
    dE = NULL, dIp = NULL, dIs = NULL, dIa = NULL, dH = NULL, dC = NULL, contact = NULL,
    u = NULL, y = NULL, fIp = NULL, fIa = NULL, fIs = NULL, rho = NULL, tau = NULL,
    seed_times = NULL, dist_seed_ages = NULL, observer = NULL, schedule = NULL)
{
    # Check parameters
    if (length(mat_locations) != length(dem_locations)) {
        if (length(mat_locations) == 1) {
            mat_locations = rep_len(mat_locations, length(dem_locations));
        } else {
            stop("dem_locations and mat_locations must have the same length; or mat_locations can have length 1 and will be used for all dem_locations.");
        }
    }
    
    # Get population parameters.
    pop = list();
    for (i in seq_along(dem_locations))
    {
        pop[[i]] = cm_build_pop_SEI3R(dem_locations[i], mat_locations[i],
            dE = dE, dIp = dIp, dIs = dIs, dIa = dIa, dH = dH, dC = dC, contact = contact,
            u = u, y = y, fIp = fIp, fIa = fIa, fIs = fIs, rho = rho, tau = tau,
            seed_times = seed_times, dist_seed_ages = dist_seed_ages, observer = observer, schedule = schedule);
    }

    # Build simulation parameters around this population.
    parameters = cm_base_parameters_SEI3R(n_groups, pop);
    parameters$date0 = date_start;
    parameters$time0 = 0;
    parameters$time1 = date_end;
    parameters$deterministic = deterministic;
    parameters$processes = processes;
    
    return (parameters)
}

# Get regions for the UK.
cm_uk_locations = function(country, level) {
    # Check country code
    country = toupper(country);
    if (country == "UK") { 
        country = "EWSN";
    }
    if (!country %like% "^(UK|[EWSN]+)$") {
        stop("country must be UK, or a combination of E, W, S, and/or N.");
    }
    
    # Interpret level
    level = as.integer(level);
    if (level < 0 | level > 4) {
        stop("level must be 0, 1, 2, 3, or 4");
    }
    
    if (level == 0) {
        if (country != "EWSN") {
            stop("For level 0, country must be UK.");
        }
        return ("UK | UNITED KINGDOM");
    } else if (level == 1) {
        gE = "Country";
        gW = "Country";
        gS = "Country";
        gN = "Country";
    } else if (level == 2) {
        gE = "Region";
        gW = "Country";
        gS = "Country";
        gN = "Country";
    } else if (level == 3) {
        gE = c("Metropolitan County", "County", "Unitary Authority", "London Borough");
        gW = "Unitary Authority";
        gS = "Council Area";
        gN = "Local Government District";
    } else if (level == 4) {
        gE = c("Metropolitan District", "Non-metropolitan District", "Unitary Authority", "London Borough");
        gW = "Unitary Authority";
        gS = "Council Area";
        gN = "Local Government District";
    }
    
    # Extract locations
    locs = NULL;
    if (country %like% "E") { locs = c(locs, cm_structure_UK[Code %like% "^E" & Geography1 %in% gE, Name]); }
    if (country %like% "W") { locs = c(locs, cm_structure_UK[Code %like% "^W" & Geography1 %in% gW, Name]); }
    if (country %like% "S") { locs = c(locs, cm_structure_UK[Code %like% "^S" & Geography1 %in% gS, Name]); }
    if (country %like% "N") { locs = c(locs, cm_structure_UK[Code %like% "^N" & Geography1 %in% gN, Name]); }
    return (paste0("UK | ", locs));
}

