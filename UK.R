# - - - - - - - - - - - - - - - - - - - - - - - 
# UK model: load data and analyse scenarios
# - - - - - - - - - - - - - - - - - - - - - - - 

library(rlang)
library(stringr)

# Set path
# Set this path to the base directory of the repository.
covid_uk_path = "~/Dropbox/COVID-UK"

# covidm options
cm_path = paste0(covid_uk_path, "/covidm/");
if (grepl(Sys.info()["user"], pattern = "^adamkuchars(ki)?$")) { cm_path = "~/Documents/GitHub/covid-uk/covidm/" }
source(paste0(cm_path, "/R/covidm.R"))

# build parameters for entire UK, for setting R0.
parametersUK1 = cm_parameters_SEI3R(cm_uk_locations("UK", 0), 
                                    dE  = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p,
                                    dIp = cm_delay_gamma(1.5, 4.0, t_max = 60, t_step = 0.25)$p,
                                    dIs = cm_delay_gamma(3.5, 4.0, t_max = 60, t_step = 0.25)$p,
                                    dIa = cm_delay_gamma(5.0, 4.0, t_max = 60, t_step = 0.25)$p,
                                    deterministic = F);

# build parameters for regions of UK, down to the county level (level 3).
locations = cm_uk_locations("UK", 3);
parameters = cm_parameters_SEI3R(locations, date_start = "2020-01-29", date_end = "2021-12-31",
                                 dE  = cm_delay_gamma(4.0, 4.0, t_max = 60, t_step = 0.25)$p, # 6.5 day serial interval.
                                 dIp = cm_delay_gamma(1.5, 4.0, t_max = 60, t_step = 0.25)$p, # 1.5 days w/o symptoms
                                 dIs = cm_delay_gamma(3.5, 4.0, t_max = 60, t_step = 0.25)$p, # 5 days total of infectiousness
                                 dIa = cm_delay_gamma(5.0, 4.0, t_max = 60, t_step = 0.25)$p, # 5 days total of infectiousness here as well.
                                 deterministic = F);

# Split off the elderly (70+, age groups 15 and 16) so their contact matrices can be manipulated separately
parameters = cm_split_matrices_ex_in(parameters, 15);

# Create additional matrix for child-elderly contacts
for (j in seq_along(parameters$pop))
{
  # Recover home/other contact matrix
  mat_ref = parameters$pop[[j]]$matrices[[1]] + parameters$pop[[j]]$matrices[[4]] + 
    parameters$pop[[j]]$matrices[[5]] + parameters$pop[[j]]$matrices[[8]];
  
  gran = 5/7; # adjustment for weekdays only.
  N = nrow(mat_ref);
  popsize = parameters$pop[[j]]$size;
  mat = matrix(0, ncol = N, nrow = N);
  
  # Add child-grandparent contacts: under 15s to 55+s
  for (a in 1:3) {
    dist = c(rep(0, 10 + a), mat_ref[a, (11 + a):N]);
    dist = dist/sum(dist);
    mat[a, ] = mat[a, ] + gran * dist;
    mat[, a] = mat[, a] + (gran * dist) * (popsize[a] / popsize);
  }
  
  # Add child-grandparent contact matrix to population
  parameters$pop[[j]]$matrices$gran = mat;
  parameters$pop[[j]]$contact = c(parameters$pop[[j]]$contact, 0);
}


# Health burden processes
probs = fread(
  "Age,Prop_symptomatic,IFR,Prop_inf_hosp,Prop_inf_critical,Prop_critical_fatal,Prop_noncritical_fatal,Prop_symp_hospitalised,Prop_hospitalised_critical
  10,0.66,8.59E-05,0.002361009,6.44E-05,0.5,0,0,0.3
  20,0.66,0.000122561,0.003370421,9.19E-05,0.5,9.47E-04,0.007615301,0.3
  30,0.66,0.000382331,0.010514103,0.000286748,0.5,0.001005803,0.008086654,0.3
  40,0.66,0.000851765,0.023423527,0.000638823,0.5,0.001231579,0.009901895,0.3
  50,0.66,0.001489873,0.0394717,0.001117404,0.5,0.002305449,0.018535807,0.3
  60,0.66,0.006933589,0.098113786,0.005200192,0.5,0.006754596,0.054306954,0.3
  70,0.66,0.022120421,0.224965092,0.016590316,0.5,0.018720727,0.150514645,0.3
  80,0.66,0.059223786,0.362002579,0.04441784,0.5,0.041408882,0.332927412,0.3
  100,0.66,0.087585558,0.437927788,0.065689168,0.5,0.076818182,0.617618182,0.3")

reformat = function(P)
{
  # 70-74,3388.488  75-79,2442.147  80-84,1736.567  85-89,1077.555  90-94,490.577  95-99,130.083  100+,15.834
  x = c(P[1:7], weighted.mean(c(P[8], P[9]), c(3388.488 + 2442.147, 1736.567 + 1077.555 + 490.577 + 130.083 + 15.834)));
  return (rep(x, each = 2))
}

P.icu_symp     = reformat(probs[, Prop_symp_hospitalised * Prop_hospitalised_critical]);
P.nonicu_symp  = reformat(probs[, Prop_symp_hospitalised * (1 - Prop_hospitalised_critical)]);
P.death_icu    = reformat(probs[, Prop_critical_fatal]);
P.death_nonicu = reformat(probs[, Prop_noncritical_fatal]);


burden_processes = list(
  list(source = "Ip", type = "multinomial", names = c("to_icu", "to_nonicu", "null"), report = c("", "", ""),
       prob = matrix(c(P.icu_symp, P.nonicu_symp, 1 - P.icu_symp - P.nonicu_symp), nrow = 3, ncol = 16, byrow = T),
       delays = matrix(c(cm_delay_gamma(7, 7, 60, 0.25)$p, cm_delay_gamma(7, 7, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 3, byrow = T)),
  
  list(source = "to_icu", type = "multinomial", names = "icu", report = "p",
       prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_gamma(10, 10, 60, 0.25)$p, nrow = 1, byrow = T)),
  
  list(source = "to_nonicu", type = "multinomial", names = "nonicu", report = "p",
       prob = matrix(1, nrow = 1, ncol = 16, byrow = T),
       delays = matrix(cm_delay_gamma(8, 8, 60, 0.25)$p, nrow = 1, byrow = T)),
  
  list(source = "Ip", type = "multinomial", names = c("death", "null"), report = c("o", ""),
       prob = matrix(c(P.death_nonicu, 1 - P.death_nonicu), nrow = 2, ncol = 16, byrow = T),
       delays = matrix(c(cm_delay_gamma(22, 22, 60, 0.25)$p, cm_delay_skip(60, 0.25)$p), nrow = 2, byrow = T))
)
parameters$processes = burden_processes

# Observer for lockdown scenarios
observer_lockdown = function(lockdown_trigger) function(time, dynamics)
{
  # Get current icu prevalence
  icu_prevalence = dynamics[t == time, sum(icu_p)];
  
  # Determine lockdown trigger
  trigger = lockdown_trigger;
  
  # If ICU prevalence exceeds a threshold, turn on lockdown
  if (icu_prevalence >= trigger) {
    return (list(csv = paste(time, "trace_lockdown", "All", 2, sep = ","),
                 changes = list(contact_lowerto = c(1, 0.1, 0.1, 0.1,  1, 0.1, 0.1, 0.1,  1))));
  } else  {
    return (list(csv = paste(time, "trace_lockdown", "All", 1, sep = ","),
                 changes = list(contact_lowerto = c(1, 1, 1, 1, 1, 1, 1, 1, 1))));
  }
  return (list(csv = paste(time, "trace_lockdown", "All", 1, sep = ",")))
}

# Load age-varying symptomatic rate
covid_scenario = qread(paste0(covid_uk_path, "/data/2-linelist_symp_fit_fIa0.5.qs"));

# Identify London boroughs for early seeding, and regions of each country for time courses
london = cm_structure_UK[match(str_sub(locations, 6), Name), Geography1 %like% "London"]
england = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^E" & !(Geography1 %like% "London")]
wales = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^W"]
scotland = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^S"]
nireland = cm_structure_UK[match(str_sub(locations, 6), Name), Code %like% "^N"]
westmid = cm_structure_UK[match(str_sub(locations, 6), Name), Name == "West Midlands (Met County)"]
cumbria = cm_structure_UK[match(str_sub(locations, 6), Name), Name == "Cumbria"]

save = function(run)
{
  # if (analysis == 3) {
  #     filename = paste0("~/Dropbox/COVID-UK Storage/", run$dynamics$scenario[1], "-", run$dynamics$run[1], ".qs");
  #     cm_save(run, filename);
  # }
}

add_totals = function(run, totals)
{
  regions = run$dynamics[, unique(population)];
  
  # totals by age
  totals0 = run$dynamics[, .(total = sum(value)), by = .(scenario, run, compartment, group)];
  return (rbind(totals, totals0))
}

add_dynamics = function(run, dynamics, iv)
{
  regions = run$dynamics[, unique(population)];
  
  interv = data.table(scenario = run$dynamics$scenario[1], run = run$dynamics$run[1], t = unique(run$dynamics$t), 
                      compartment = "trace_school", region = "All", value = unlist(iv$trace_school));
  if (!is.null(iv$trace_intervention)) {
    interv = rbind(interv,
                   data.table(scenario = run$dynamics$scenario[1], run = run$dynamics$run[1], t = unique(run$dynamics$t), 
                              compartment = "trace_intervention", region = "All", value = unlist(iv$trace_intervention)));
  } else {
    interv = rbind(interv,
                   data.table(scenario = run$dynamics$scenario[1], run = run$dynamics$run[1], t = unique(run$dynamics$t), 
                              compartment = "trace_intervention", region = "All", value = 1));
  }
  
  csvlines = NULL;
  if (nchar(run$csv[[1]]) > 0) {
    csvlines = fread(run$csv[[1]], header = F);
    csvlines = cbind(run$dynamics$scenario[1], run$dynamics$run[1], csvlines);
    names(csvlines) = c("scenario", "run", "t", "compartment", "region", "value");
    csvlines = unique(csvlines);
  }
  
  # time courses
  return (rbind(dynamics,
                run$dynamics[population %in% locations[westmid],  .(region = "West Midlands",    value = sum(value)), by = .(scenario, run, t, compartment)],
                run$dynamics[population %in% locations[cumbria],  .(region = "Cumbria",          value = sum(value)), by = .(scenario, run, t, compartment)],
                run$dynamics[population %in% locations[london],   .(region = "London",           value = sum(value)), by = .(scenario, run, t, compartment)],
                run$dynamics[population %in% locations[england],  .(region = "England",          value = sum(value)), by = .(scenario, run, t, compartment)],
                run$dynamics[population %in% locations[wales],    .(region = "Wales",            value = sum(value)), by = .(scenario, run, t, compartment)],
                run$dynamics[population %in% locations[scotland], .(region = "Scotland",         value = sum(value)), by = .(scenario, run, t, compartment)],
                run$dynamics[population %in% locations[nireland], .(region = "Northern Ireland", value = sum(value)), by = .(scenario, run, t, compartment)],
                run$dynamics[,                                    .(region = "United Kingdom",   value = sum(value)), by = .(scenario, run, t, compartment)],
                interv,
                csvlines
  ))
}

#############
# MAIN CODE #
#############

argv = commandArgs(trailingOnly = T);
argc = length(argv);
if (argc != 2) {
  stop("Must provide two arguments: analysis set and number of runs.");
}
analysis = as.numeric(argv[argc-1]);
n_runs = as.numeric(argv[argc]);

if (analysis == 1) {
  # Define school terms, base versus intervention (both same here)
  school_close_b =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_b = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  school_close_i =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_i = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  
  # Define interventions to be used
  interventions = list(
    `School Closures`   = list(contact = c(1.0, 1.0, 0.0, 1.0,  1.0, 1.0, 0.0, 1.0,  0)),
    `Social Distancing` = list(contact = c(1.0, 0.5, 1.0, 0.5,  1.0, 0.5, 1.0, 0.5,  0)),
    `Elderly Shielding` = list(contact = c(1.0, 1.0, 1.0, 1.0,  1.0, 0.25, 1.0, 0.25,  0)),
    `Self-Isolation`    = list(fIs = rep(0.65, 16)),
    `Combination`       = list(contact = c(1.0, 0.5, 0.0, 0.5,  1.0, 0.25, 0.0, 0.25,  0), fIs = rep(0.65, 16))
  );
  
  # Set options
  option.trigger = "national";
  option.duration = 7 * 12;
  option.lockdown = NA;
  option.intervention_shift = 0;
} else if (analysis == 2.1) {
  # Define school terms, base versus intervention (both same here)
  school_close_b =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_b = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  school_close_i =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_i = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  
  # Define interventions to be used
  interventions = list(
    `Combination`       = list(contact = c(1.0, 0.5, 0.0, 0.5,  1.0, 0.25, 0.0, 0.25,  0), fIs = rep(0.65, 16))
  );
  
  # Set options
  option.trigger = "national";
  option.duration = 7 * 12;
  option.lockdown = NA;
  option.intervention_shift = c(0, 14, 28, 56);
} else if (analysis == 2.2) {
  # Define school terms, base versus intervention (both same here)
  school_close_b =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_b = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  school_close_i =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_i = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  
  # Define interventions to be used
  interventions = list(
    `Combination`       = list(contact = c(1.0, 0.5, 0.0, 0.5,  1.0, 0.25, 0.0, 0.25,  0), fIs = rep(0.65, 16))
  );
  
  # Set options
  option.trigger = "local";
  option.duration = 7 * 12;
  option.lockdown = NA;
  option.intervention_shift = c(0, 14, 28, 56);
} else if (analysis == 3) {
  # Define school terms, base versus intervention (schools close from 23 March)
  school_close_b =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_b = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  school_close_i =  c("2020-2-16", "2020-3-23",                           "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_i = c("2020-2-22",                           "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  
  # Define interventions to be used
  interventions = list(
    `Intensive Interventions` = list(contact = c(1, 0.655, 1, 0.59155,  1, 0.25, 1, 0.157375,  0), fIs = rep(0.65, 16))
  );
  
  # Set options
  option.trigger = "2020-03-17";
  option.duration = 364;
  option.lockdown = c(NA, 1000, 2000, 5000);
  option.intervention_shift = 0;
} else if (analysis == 4) {
  # Define school terms, base versus intervention (both same here)
  school_close_b =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_b = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  school_close_i =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_i = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  
  # Define interventions to be used
  interventions = list(
    `Intensive`                    = list(contact = c(1, 0.655, 1, 0.59155,  1, 0.25, 1, 0.157375,  0.0), fIs = rep(0.65, 16)),
    `Intensive + School`           = list(contact = c(1, 0.655, 0, 0.59155,  1, 0.25, 0, 0.157375,  0.0), fIs = rep(0.65, 16)),
    `Intensive + School + G20`     = list(contact = c(1, 0.655, 0, 0.59155,  1, 0.25, 0, 0.157375,  0.2), fIs = rep(0.65, 16)),
    `Intensive + School + G50`     = list(contact = c(1, 0.655, 0, 0.59155,  1, 0.25, 0, 0.157375,  0.5), fIs = rep(0.65, 16)),
    `Intensive + School + G100`    = list(contact = c(1, 0.655, 0, 0.59155,  1, 0.25, 0, 0.157375,  1.0), fIs = rep(0.65, 16))
  );
  
  # Set options
  option.trigger = "2020-03-17";
  option.duration = 125;
  option.lockdown = NA;
  option.intervention_shift = 0;
  parameters$time1 = "2020-07-20";
} else if (analysis == 6) {
  # Define school terms, base versus intervention (both same here)
  school_close_b =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_b = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  school_close_i =  c("2020-2-16", "2020-4-05", "2020-5-24", "2020-7-22", "2020-10-25", "2020-12-20", "2021-02-14", "2021-04-01", "2021-05-30", "2021-07-25");
  school_reopen_i = c("2020-2-22", "2020-4-18", "2020-5-30", "2020-9-01", "2020-10-31", "2021-01-02", "2021-02-20", "2021-04-17", "2021-06-05", "2021-09-01");
  
  # Define interventions to be used
  interventions = list(
    `Background`               = list(contact = c(1, 1, 0, 1,          1, 0.25, 0, 0.25,  0), fIs = rep(0.65, 16)),
    `Background + 0% Sports`   = list(contact = c(1, 1, 0, 1 - 0.0041, 1, 0.25, 0, 0.25,  0), fIs = rep(0.65, 16)),
    `Background + 25% Leisure` = list(contact = c(1, 1, 0, 1 - 0.362,  1, 0.25, 0, 0.25,  0), fIs = rep(0.65, 16))
  );
  
  # Set options
  option.trigger = "2020-03-17";
  option.duration = 168;
  option.lockdown = NA;
  option.intervention_shift = 0;
  parameters$time1 = "2020-09-01";
}

# Pick R0s 
set.seed(9876);
R0s = rnorm(n_runs, mean = 2.675739, sd = 0.5719293)

# Do runs
dynamics = data.table()
totals = data.table()
print(Sys.time())
set.seed(1234567);

for (r in 1:n_runs) {
  cat(paste0(r, ": R0 = ", R0s[r], "\n"));
  
  # 1. Pick age-varying symptomatic rate
  covy = unname(unlist(covid_scenario[sample.int(nrow(covid_scenario), 1), f_00:f_70]));
  covy = rep(covy, each = 2);
  
  # 2. Calculate R0 adjustment needed
  parametersUK1$pop[[1]]$y = covy;
  u_adj = R0s[r] / cm_calc_R0(parametersUK1, 1);
  
  # 3. Pick seeding times
  seed_start = ifelse(london, sample(0:6, length(london), replace = T), sample(0:20, length(london), replace = T));
  
  # 4. Do base model
  
  # 4a. Set parameters
  params = duplicate(parameters);
  for (j in seq_along(params$pop)) {
    params$pop[[j]]$u = params$pop[[j]]$u * u_adj;
    params$pop[[j]]$y = covy;
    params$pop[[j]]$seed_times = rep(seed_start[j] + 0:27, each = 2);
    params$pop[[j]]$dist_seed_ages = cm_age_coefficients(25, 50, 5 * 0:16);
  }
  
  # CALCULATE IMPACT ON R0
  if (analysis == 5) {
    interventions = list(
      `Base`                      = list(),
      `School Closures`           = list(contact = c(1.0, 1.0, 0.0, 1.0,  1.0, 1.0, 0.0, 1.0,  0)),
      `Social Distancing`         = list(contact = c(1.0, 0.5, 1.0, 0.5,  1.0, 0.5, 1.0, 0.5,  0)),
      `Elderly Shielding`         = list(contact = c(1.0, 1.0, 1.0, 1.0,  1.0, 0.25, 1.0, 0.25,  0)),
      `Self-Isolation`            = list(fIs = rep(0.65, 16)),
      `Combination`               = list(contact = c(1.0, 0.5, 0.0, 0.5,  1.0, 0.25, 0.0, 0.25,  0), fIs = rep(0.65, 16)),
      `Intensive, schools open`   = list(contact = c(1, 0.655, 1, 0.59155,  1, 0.25, 1, 0.157375,  0), fIs = rep(0.65, 16)),
      `Intensive, schools closed` = list(contact = c(1, 0.655, 0, 0.59155,  1, 0.25, 0, 0.157375,  0), fIs = rep(0.65, 16)),
      `Lockdown`                  = list(contact = c(1, 0.1, 0.1, 0.1,  1, 0.1, 0.1, 0.1,  0), fIs = rep(0.65, 16))
    );
    
    for (i in seq_along(interventions))
    {
      iR0s = rep(0, length(params$pop));
      iweights = rep(0, length(params$pop));
      for (j in seq_along(params$pop))
      {
        for (k in seq_along(interventions[[i]]))
        {
          params$pop[[j]][[names(interventions[[i]])[k]]] = interventions[[i]][[k]];
        }
        iR0s[j] = cm_calc_R0(params, j);
        iweights[j] = sum(params$pop[[j]]$size);
      }
      
      weighted_R0 = weighted.mean(iR0s, iweights);
      dynamics = rbind(dynamics, data.table(run = r, scenario = names(interventions)[i], R0 = weighted_R0));
    }
    
    next;
  }
  
  # 4b. Set school terms
  iv = cm_iv_build(params)
  cm_iv_set(iv, school_close_b, school_reopen_b, contact = c(1, 1, 0, 1,  1, 1, 0, 1,  1), trace_school = 2);
  params = cm_iv_apply(params, iv);
  
  # 4c. Run model
  run = cm_simulate(params, 1, r);
  run$dynamics[, run := r];
  run$dynamics[, scenario := "Base"];
  run$dynamics[, R0 := R0s[r]];
  save(run);
  totals = add_totals(run, totals);
  dynamics = add_dynamics(run, dynamics, iv);
  peak_t = run$dynamics[compartment == "cases", .(total_cases = sum(value)), by = t][, t[which.max(total_cases)]];
  peak_t_bypop = run$dynamics[compartment == "cases", .(total_cases = sum(value)), by = .(t, population)][, t[which.max(total_cases)], by = population]$V1;
  
  rm(run)
  gc()
  
  # 5. Run interventions
  for (i in seq_along(interventions)) {
    for (duration in option.duration) {
      for (trigger in option.trigger) {
        for (intervention_shift in option.intervention_shift) {
          for (lockdown in option.lockdown) {
            cat(paste0(names(interventions)[i], "...\n"))
            
            # 5a. Make parameters and adjust R0
            params = duplicate(parameters);
            for (j in seq_along(params$pop)) {
              params$pop[[j]]$u = params$pop[[j]]$u * u_adj;
              params$pop[[j]]$y = covy;
              if (!is.na(lockdown)) {
                params$pop[[j]]$observer = observer_lockdown(lockdown);
              }
            }
            
            # 5b. Set interventions
            if (trigger == "national") {
              intervention_start = peak_t - duration / 2 + intervention_shift;
            } else if (trigger == "local") {
              intervention_start = peak_t_bypop - duration / 2 + intervention_shift;
            } else {
              intervention_start = as.numeric(ymd(trigger) - ymd(params$date0));
            }
            
            if (trigger == "local") {
              # Apply interventions to one population at a time.
              for (pi in seq_along(params$pop)) {
                ymd_start = ymd(params$date0) + intervention_start[pi];
                ymd_end = ymd_start + duration - 1;
                iv = cm_iv_build(params)
                cm_iv_set(iv, school_close_i, school_reopen_i, contact = c(1, 1, 0, 1,  1, 1, 0, 1,  1), trace_school = 2);
                cm_iv_set(iv, ymd_start, ymd_end, interventions[[i]]);
                cm_iv_set(iv, ymd_start, ymd_end, trace_intervention = 2);
                params = cm_iv_apply(params, iv, pi);
              }
            } else {
              # Apply interventions to entire population.
              ymd_start = ymd(params$date0) + intervention_start;
              ymd_end = ymd_start + duration - 1;
              iv = cm_iv_build(params)
              cm_iv_set(iv, school_close_i, school_reopen_i, contact = c(1, 1, 0, 1,  1, 1, 0, 1,  1), trace_school = 2);
              cm_iv_set(iv, ymd_start, ymd_end, interventions[[i]]);
              cm_iv_set(iv, ymd_start, ymd_end, trace_intervention = 2);
              params = cm_iv_apply(params, iv);
            }
            
            # 5c. Run model
            run = cm_simulate(params, 1, r);
            
            tag = "";
            if (length(option.duration) > 1)            { tag = paste0(tag, " ", duration + 1, " day"); }
            if (length(option.lockdown) > 1)            { tag = paste0(tag, " ", ifelse(lockdown >= 0, lockdown, "variable"), " lockdown"); }
            if (length(option.trigger) > 1)             { tag = paste0(tag, " ", trigger, " trigger"); }
            if (length(option.intervention_shift) > 1)  { tag = paste0(tag, " ", intervention_shift, " shift"); }
            
            run$dynamics[, run := r];
            run$dynamics[, scenario := paste0(names(interventions)[i], tag)];
            run$dynamics[, R0 := R0s[r]];
            save(run);
            totals = add_totals(run, totals);
            dynamics = add_dynamics(run, dynamics, iv);
            
            rm(run)
            gc()
          }
        }
      }
    }
  }
}
cm_save(totals, paste0(covid_uk_path, analysis, "-totals.qs"));
cm_save(dynamics, paste0(covid_uk_path, analysis, "-dynamics.qs"));
print(Sys.time())