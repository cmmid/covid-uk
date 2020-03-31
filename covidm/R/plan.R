# processes


processes = list(
    list(name = "to_hospital", source = "Is", prob = 0.5, report = NA, delay = cm_dist_gamma(12, 4, 21, 0.25???)),    # 12 day delay from symptom onset to hospital
    list(name = "non-ICU", source = "to-hospital", prob = 0.96, report = "prevalence", delay = cm_dist_gamma(7, 4)),  # of those who present at hospital, 96% are non-ICU
    list(name = "ICU", source = "to-hospital", prob = 0.04, report = "prevalence", delay = cm_dist_gamma(14, 4))      # and 4% are ICU
)


trees = list(
    list(name = "to_hospital", source = "Is", prob = 0.5, report = NA, delay = cm_dist_gamma(12, 4, 21, 0.25???),    # 12 day delay from symptom onset to hospital
        sub = list(
            list(name = "non-ICU", source = "to-hospital", prob = 0.96, report = "prevalence", delay = cm_dist_gamma(7, 4)),  # of those who present at hospital, 96% are non-ICU
            list(name = "ICU", source = "to-hospital", prob = NA, report = "prevalence", delay = cm_dist_gamma(14, 4))      # the rest are ICU
        )
    )
)

# Here's what I would like to see as an interface. At the moment this is all very hypothetical.

cm_parameters = function( some common options here . . . ) { }

# this would return a list() of parameters with some sensible defaults, 
# according to parameters passed to it. we would also have

cm_parameters_england()
cm_parameters_wuhan()

# etc., for the more common use cases.

# actually, a common use case is to compare different interventions.
# so we might want to have parameters actually be a list with elements 1, 2, 3 etc or with named elements,
# so that this can represent a set of scenarios. The number of scenarios would be set in cm_parameters

# A function

cm_make_interventions = function(parameters, school_terms, interventions, ...){}

# would take a parameters object and add interventions to it, which include school terms etc. This should
# warn if there are any existing interventions as they will almost definitely get overwritten.


# the function
cm_simulate = function(parameters, runs = 1, seed = NULL, observer = NULL, all_ages = T, all_regions = T) {}

# would run the simulation with the given parameters, "runs" times, and if a seed is provided
# then it would use that seed to start off. If all_ages or all_regions are F, then these are summed
# up before being returned in the dynamics (to save space...). This could also be done by the user with

cm_remove_ages = function(model_run) {}

# or

cm_remove_regions = function(model_run) {}

# or both.



# cm_simulate returns a list of type model_run, which is of the form

list(
    parameters,
    dynamics,
    seed,
    observer,
    model_version
)

# covidm is a version number in case it is needed.

# in cm_simulate and model_run, observer is either NULL or a function of the form
observer = function(scenario, run, t, dynamics) {}
# where scenario is which set of parameters is being run, run is which run, t is what day it is, dynamics is today's dynamics.
# it should return some kind of observer state as well as something formally similar to temp_contact_changes or temp_dist_changes
# to trigger one of these things. One should be very careful of using both observer changes and scheduled changes as these
# will not pay attention to what the other is doing.

# CHANGES TO PARAMETERS
# - instead of contact_matrix, populations should have a list of matrices, which by default are summed at sim start.
# then scheduled / observer changes can make reference to this list.
# - make parameter names consistent with new terminology.
# - group instead of age_group, named in matrix, named in dynamics
# - named populations in dynamics
# - maybe:
#            P1 -> I1 
#  S -> E -<          > R
#            P2 -> I2
# - what about alternative model structures/populations and their associated parameters?
# within parameters, populations will have to have types, which can in principle have their own parameters.

# there should be various plotting functions operating on model_run


#
    
