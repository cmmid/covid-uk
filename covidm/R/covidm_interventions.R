# covidm_interventions.R
# setting up intervention scenarios

# set up interventions
# TODO remove reliance on cm_translate_parameters, for efficiency?
# TODO I don't like that this has to set contact and fIs.
cm_iv_build = function(parameters)
{
    p = cm_translate_parameters(parameters)
    data.table(date = (ymd(p$date0) + p$time0):(ymd(p$date0) + p$time1));
}

cm_iv_checkset = function(iv, what, change)
{
    if (!(what %in% names(iv))) {
        iv[, (what) := list(rep(1, length(change)))];
    }
}

# set up school breaks
cm_iv_school_breaks = function(iv, ymd_break_start, ymd_break_end, sf = 0)
{
    cm_iv_set(iv, ymd_break_start, ymd_break_end, contact = c(1, 1, sf, 1));
}

# set up other contact interventions
# TODO add a trace component
# TODO min option instead of *
cm_iv_contact = function(iv, ymd_iv_first_day, ymd_iv_last_day, cf)
{
    cm_iv_set(iv, ymd_iv_first_day, ymd_iv_last_day, contact = cf);
}

# apply arbitrary interventions
cm_iv_set = function(iv, ymd_iv_first_day, ymd_iv_last_day, ...)
{
    changes = list(...);
    if (length(changes) == 1 & is.list(changes[[1]]) & is.null(names(changes))) {
        changes = changes[[1]];
    }
    
    # Apply each intervention
    for (i in seq_along(changes))
    {
        # Check column is present, initialize if it is not
        col = names(changes)[i];
        cm_iv_checkset(iv, col, changes[[i]]);
        
        # Check start & end dates are same length
        if (length(ymd_iv_first_day) != length(ymd_iv_last_day)) {
            stop("length(ymd_iv_first_day) != length(ymd_iv_last_day)");
        }
        
        # Apply intervention for each time period
        for (j in seq_along(ymd_iv_first_day))
        {
            # Convert start and end dates
            t0 = ymd(ymd_iv_first_day[j]);
            t1 = ymd(ymd_iv_last_day[j]);
            
            # Apply intervention
            iv[, (col) := as(get(col), "list")];
            iv[date >= t0 & date <= t1, (col) := lapply(get(col), function(x) x * changes[[i]])];
        }
    }
}

# apply interventions to a parameter set
# TODO merge with existing schedule
cm_iv_apply = function(parameters, iv, populations = -1)
{
    # Make schedule of interventions
    schedule = list();
    for (ivr in seq_len(nrow(iv)))
    {
        # Skip this row if nothing has changed
        if (ivr == 1) {
        } else if (any(unlist(iv[ivr, .SD, .SDcols = -"date"]) != unlist(iv[ivr - 1, .SD, .SDcols = -"date"]))) {
        } else next;
        
        # Set up changes
        change = list(t = as.numeric(as_date(iv[ivr, date]) - ymd(parameters$date0)));
        for (iv_col in 2:ncol(iv))
        {
            if (!(names(iv)[iv_col] %like% "^trace")) {
                change[[length(change) + 1]] = iv[ivr, ..iv_col][[1]][[1]];
                names(change)[length(change)] = names(iv)[iv_col];
            }
        }
        
        # Add to schedule
        schedule[[length(schedule) + 1]] = change;
    }

    # Apply to each population requested
    if (populations[1] < 0) {
        populations = seq_along(parameters$pop);
    }
    for (pi in populations)
    {
        if (length(parameters$pop[[pi]]$schedule) == 0) {
            parameters$pop[[pi]]$schedule = schedule;
        } else {
            stop(paste0("Schedule already set for population ", pi, ". (Merging of schedules not yet implemented.)"));
        }
    }
    
    return (parameters);
}
