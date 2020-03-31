# covidm_misc.R
# various useful functions

# Wrappers around qs to save/load an R object.
cm_save = function(x, file) {
    qsave(x, file, "balanced");
}

cm_load = function(file) {
    qread(file);
}

# Safe sampling (unlike sample, works even if x is only one number)
cm_resample = function(x, ...) x[sample.int(length(x), ...)]

# Return the mean and estimated highest density interval from a set of samples x.
cm_mean_hdi = function(x, credMass = 0.95)
{
    h = hdi(x, credMass);
    m = mean(x);
    return (list(mean = m, lower = h[[1]], upper = h[[2]]))
}

# construct a delay distribution following a gamma distribution with mean mu and shape parameter shape.
cm_delay_gamma = function(mu, shape, t_max, t_step)
{
    scale = mu / shape;
    t_points = seq(0, t_max, by = t_step);
    heights = pgamma(t_points + t_step/2, shape, scale = scale) - 
        pgamma(pmax(0, t_points - t_step/2), shape, scale = scale);
    return (data.table(t = t_points, p = heights / sum(heights)))
}

# construct a delay distribution that effectively skips the compartment
cm_delay_skip = function(t_max, t_step)
{
    t_points = seq(0, t_max, by = t_step);
    return (data.table(t = t_points, p = c(1, rep(0, t_max / t_step))))
}

# smoothly interpolate between points (x0, y0) and (x1, y1) using cosine interpolation.
# for x < x0, returns y0; for x > x1, returns y1; for x0 < x < x1, returns the cosine interpolation between y0 and y1
cm_interpolate_cos = function(x, x0, y0, x1, y1)
{
    ifelse(x < x0, y0, ifelse(x > x1, y1, y0 + (y1 - y0) * (0.5 - 0.5 * cos(pi * (x - x0) / (x1 - x0)))))
}

# for a set of age group bounds age_bounds (including the lower and upper bound as well as intermediate bounds),
# return a vector giving length(age_bounds) - 1 numerical entries between 0 and 1 inclusive, depending on how much
# each age group lies within the range age_min to age_max.
cm_age_coefficients = function(age_min, age_max, age_bounds)
{
    x = rep(0, length(age_bounds) - 1);
    for (ag in 1:(length(age_bounds) - 1)) {
        ag_lower = age_bounds[ag];
        ag_upper = age_bounds[ag + 1];
        ag_0 = max(age_min, ag_lower);
        ag_1 = min(age_max, ag_upper);
        x[ag] = max(0, (ag_1 - ag_0) / (ag_upper - ag_lower));
    }
    return (x);
}


# look up pop_region and return vector with the number of people in different age bands (5 year size)
cm_make_population = function(pop_region, n_age_groups)
{
    pop_size = cm_populations[name == pop_region, .(pop = round(sum(f, m) * 1000)), by = age]$pop;
    
    if (length(pop_size) == 0) {
        stop(paste("Could not find region", pop_region, "in cm_populations."));
    } else if (anyDuplicated(cm_populations[name == pop_region, age])) {
        stop(paste("Duplicate entries for region", pop_region, "in cm_populations."));
    }
    
    if (n_age_groups == 1) {
        pop_size = sum(pop_size);
    } else {
        pop_size = c(pop_size[1:(n_age_groups - 1)], sum(pop_size[n_age_groups:length(pop_size)]));
    }
    return (pop_size);
}


# expects dynamic results from a single run; the start date of the simulation; 
# the date upon which measurement starts; dates marking the end of each measurement period (i.e. n bounds for n measurements);
# and age group outer bounds (n+1 bounds for n measurements)
cm_case_distribution = function(z, date_simulation_start, date_measurement_start, dates_measurement_end, age_bounds, compart = "cases")
{
    # do date magic
    dss = ymd(date_simulation_start);
    t_ref = as.numeric(dss);
    t_bounds = c(as.numeric(ymd(date_measurement_start)) - t_ref,
                 as.numeric(ymd(dates_measurement_end)) - t_ref);

    # do age magic
    n_age_groups = z[, max(group)];
    age_dists = list();
    age_names = c();
    for (a in 1:(length(age_bounds) - 1)) {
        age_dists[[a]] = cm_age_coefficients(age_bounds[a], age_bounds[a + 1], 5 * (0:n_age_groups));
        age_names = c(age_names, paste(age_bounds[a], "-", age_bounds[a + 1]));
    }

    # cast incidence into a usable form
    incidence = data.table(dcast(z[compartment == compart], t ~ group, fun.aggregate = sum, value.var = "value"))

    # slice up cases by date
    dist = NULL;
    for (slice in 1:length(dates_measurement_end)) {
        inc = unname(incidence[t >= t_bounds[slice] & t < t_bounds[slice + 1], colSums(.SD), .SDcols = -1]);
        t_name = paste(dss + t_bounds[slice] + min(1, slice - 1), "-", dss + t_bounds[slice + 1]);
        dist = rbind(dist, 
            data.table(
                date = rep(t_name, length(age_dists)),
                age = age_names,
                cases = sapply(age_dists, function(age_dist) sum(inc * age_dist))
            )
        )
    }
    
    # also report fraction of cases by age group
    dist[, fcases := cases / sum(cases), by = date];
    
    return (dist)
}

# Calculate R0 for a given population
cm_calc_R0 = function(parameters, population) {
    po = parameters$pop[[population]];
    dIp = sum(po$dIp * seq(0, by = parameters$time_step, length.out = length(po$dIp)));
    dIs = sum(po$dIs * seq(0, by = parameters$time_step, length.out = length(po$dIs)));
    dIa = sum(po$dIa * seq(0, by = parameters$time_step, length.out = length(po$dIa)));
    
    cm = Reduce('+', mapply(function(c, m) c * m, po$contact, po$matrices, SIMPLIFY = F));

    ngm = po$u * t(t(cm) * (
        po$y * (po$fIp * dIp + po$fIs * dIs) + 
        (1 - po$y) * po$fIa * dIa)
    )
    abs(eigen(ngm)$values[1])
}

