# - - - - - - - - - - - - - - - - - - - - - - - 
# UK model: plot outputs
# - - - - - - - - - - - - - - - - - - - - - - - 

library(cowplot)
library(stringr)
library(rlang)

# Set path
# Set this path to the base directory of the repository.
covid_uk_path = "~/Dropbox/COVID-UK"

# covidm options
cm_path = paste0(covid_uk_path, "/covidm/");
if(grepl(Sys.info()["user"], pattern = "^adamkuchars(ki)?$")){cm_path = "~/Documents/GitHub/covidm/"}
source(paste0(cm_path, "/R/covidm.R"))

# Friendly number
friendly = function(x, digits = 2)
{
    ifelse(x == 0, 0,
        ifelse(x > 1000000, paste(signif(x / 1000000, digits), "M"),
            ifelse(x > 1000, paste(signif(x / 1000, digits), "k"),
                signif(x, digits))))
}

# Friendly y axis labels
axis_friendly = function(breaks)
{
    if (max(breaks, na.rm = T) > 1000000) {
        return (ifelse(breaks == 0, 0, paste(breaks / 1000000, "M")))
    } else if (max(breaks, na.rm = T) > 1000) {
        return (ifelse(breaks == 0, 0, paste(breaks / 1000, "k")))
    }
    return (breaks)
}

# Friendly x axis labels for dates
axis_date = function(breaks)
{
    ifelse(month(breaks) == 1 | breaks == breaks[!is.na(breaks)][1], format(breaks, "%Y %b"), format(breaks, "%b"));
}


# amalgamate user compartments into more friendly ones
reflow_dynamics = function(d)
{
    d[compartment == "icu_p", compartment := "beds_icu"];
    d[compartment == "nonicu_p", compartment := "beds_nonicu"];
    d[compartment == "death_o", compartment := "deaths"];
}

reflow_totals = function(t)
{
    t[compartment == "icu_p", compartment := "beds_icu"];
    t[compartment == "nonicu_p", compartment := "beds_nonicu"];
    t[compartment == "death_o", compartment := "deaths"];
}

table_spec = fread(
"compartment, stat, time
cases, total, t
deaths, total, t
cases, peak, week
deaths, peak, week
beds_icu, peak, t
beds_nonicu, peak, t
cases, peak_time, week
trace_lockdown, lockdown_duration, t
subclinical, total, t")

median_ci = function(x, conf = 0.95)
{
    y = quantile(x, c((1 - conf) / 2, 0.5, 1 - (1 - conf) / 2));
    y = as.list(y);
    names(y) = c("lower", "median", "upper");
    return (y)
}

make_table = function(d)
{
    d[, week := t %/% 7]
    results = NULL
    for (spec in seq_len(nrow(table_spec)))
    {
        comp = table_spec[spec, compartment];
        stat = table_spec[spec, stat];
        time = table_spec[spec, time];
        
        if (stat == "total") {
            res = d[region == "United Kingdom" & compartment == comp, .(x = sum(value)), by = .(scenario, run, region)];
            res[, statistic := paste(stat, comp)];
            res = res[, median_ci(x), by = .(scenario, region, statistic)];
            stat_nice = paste("Total", comp);
        } else if (stat == "peak") {
            res = d[region == "United Kingdom" & compartment == comp, .(x = sum(value)), by = c("scenario", "run", time, "region")];
            res = res[, .(x = max(x)), by = .(scenario, run, region)];
            res[, statistic := paste(stat, comp)];
            res = res[, median_ci(x), by = .(scenario, region, statistic)];
            stat_nice = ifelse(time == "t", paste("Peak", comp, "required"),
                paste(comp, "in peak week"));
        } else if (stat == "peak_time") {
            res = d[region == "United Kingdom" & compartment == comp, .(x = sum(value)), by = c("scenario", "run", time, "region")];
            res = res[, .(x = get(time)[which.max(x)]), by = .(scenario, run, region)];
            res[, statistic := paste(stat, comp)];
            res = res[, median_ci(x), by = .(scenario, region, statistic)];
            stat_nice = paste("Time to peak", comp, ifelse(time == "t", "(days)", "(weeks)"));
        } else if (stat == "lockdown_duration") {
            if (d[compartment == comp, .N] == 0) { 
                next;
            }
            res = d[region == "All" & compartment == comp, .(x = mean(value - 1)), by = .(scenario, run, region)];
            res[, statistic := paste(stat, comp)];
            res = res[, median_ci(x), by = .(scenario, region, statistic)];
            stat_nice = paste("Proportion of time spent in", comp);
        } else if (stat == "total_end") {
            res = d[region == "United Kingdom" & compartment == comp & t == max(t), .(x = sum(value)), by = .(scenario, run, region)];
            res[, statistic := paste(stat, comp)];
            res = res[, median_ci(x), by = .(scenario, region, statistic)];
            stat_nice = paste("Number of ", comp, " at simulation end");
        } else {
            stop("Unrecognised stat.");
        }
        stat_nice = str_to_sentence(stat_nice);
        stat_nice = str_replace(stat_nice, "beds_icu", "ICU beds");
        stat_nice = str_replace(stat_nice, "beds_nonicu", "non-ICU beds");
        stat_nice = str_replace(stat_nice, "trace_lockdown", "lockdown");
        stat_nice = str_replace(stat_nice, "S", "susceptibles");
        res[, statistic := stat_nice]
        results = rbind(results, res);
    }
    
    results[, value_str := 
        paste0(friendly(median), " (", friendly(lower), "–", friendly(upper), ")")];
    results[, statistic := factor(statistic, levels = unique(statistic))]
    results[, scenario := factor(scenario, levels = unique(scenario))]
    results
}

save_table = function(tb, filename)
{
    tb2 = dcast(tb, statistic ~ scenario, value.var = "value_str");
    fwrite(tb2, filename)
}

plot_table = function(tb, nrow = NULL)
{
    ggplot(tb) +
        geom_pointrange(aes(x = scenario, y = median, ymin = lower, ymax = upper, colour = scenario), size = 0.25, fatten = 0.2) +
        facet_wrap(~statistic, scales = "free", nrow = nrow) +
        labs(x = NULL, y = NULL) +
        scale_y_continuous(labels = axis_friendly, limits = c(0, NA)) +
        theme(strip.background = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

plot_attackrate = function(t)
{
    ts = t;
    ts[, age_lower := (as.numeric(str_replace(group, "([0-9]+)(-|\\+).*", "\\1")) %/% 10) * 10];
    ts[, age_group := paste0(age_lower, "-", age_lower + 9)];
    ts[age_lower == 70, age_group := "70+"];
    ts[, age_group := factor(age_group, levels = unique(age_group))];
    
    ts[, total := sum(total), by = .(scenario, run, compartment, age_group)];
    ts = ts[compartment %in% c("cases", "deaths"), median_ci(total), by = .(scenario, compartment, age_group)];
    ts[, compartment := paste(str_to_sentence(as.character(compartment), "(thousands)"))];
    ts[, scenario := factor(scenario, levels = unique(scenario))];
    
    ggplot(ts) +
        geom_col(aes(x = age_group, y = median / 1000, fill = scenario)) +
        geom_linerange(aes(x = age_group, ymin = lower / 1000, ymax = upper / 1000), size = 0.25) +
        facet_grid(compartment~scenario, switch = "y", scales = "free") +
        labs(x = NULL, y = NULL) +
        theme(strip.background = element_blank(), strip.placement = "outside",
            axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
}

plot_example = function(d0, t, quant, ymd_start, ymd_truncate = "2050-01-01")
{
    # Copy data and process
    d = duplicate(d0[region %in% c("United Kingdom", "All")])
    d[, scenario := factor(scenario, levels = unique(scenario))];
    d[, compartment := factor(compartment, levels = unique(compartment))];

    # Choose run
    qrun = t[scenario == "Base", sum(total), by = run][, which(rank(V1) == round(1 + (.N - 1) * quant))];
    d = d[run == qrun];
    
    # Merge intervention traces
    trace_school = d[compartment == "trace_school", .(t, trace_school = value - 1, scenario)]
    d = merge(d, trace_school, by = c("t", "scenario"), all.x = T);
    d[, trace_school := trace_school * max(value), by = .(compartment)];

    trace_intervention = d[compartment == "trace_intervention", .(t, trace_intervention = value - 1, scenario)]
    d = merge(d, trace_intervention, by = c("t", "scenario"), all.x = T);
    d[, trace_intervention := trace_intervention * max(value) * 0.75, by = .(compartment)];

    if (d[compartment == "trace_lockdown", .N > 0]) {
        trace_lockdown = d[compartment == "trace_lockdown", .(t, trace_lockdown = value - 1, scenario)]
        d = merge(d, trace_lockdown, by = c("t", "scenario"), all.x = T);
        d[, trace_lockdown := trace_lockdown * max(value) * 0.5, by = .(compartment)];
    } else {
        d[, trace_lockdown := 0];
    }

    d = d[compartment %in% c("cases", "deaths", "beds_icu", "beds_nonicu")];
    
    # Give nice names
    d[compartment == "cases", compartment := "Cases"];
    d[compartment == "deaths", compartment := "Deaths"];
    d[compartment == "beds_icu", compartment := "ICU beds\nrequired"];
    d[compartment == "beds_nonicu", compartment := "Non-ICU beds\nrequired"];
    d[, compartment := factor(compartment, levels = unique(compartment))];
    
    d[, date := ymd(ymd_start) + t];

    # Plot
    ggplot(d[region == "United Kingdom" & date <= ymd(ymd_truncate)]) +
        geom_line(aes(x = date, y = value, colour = scenario), size = 0.25) +
        geom_ribbon(aes(x = date, ymin = 0, ymax = trace_school), fill = "#0000ff", alpha = 0.1) +
        geom_ribbon(aes(x = date, ymin = 0, ymax = trace_intervention), fill = "#ff0000", alpha = 0.1) +
        geom_ribbon(aes(x = date, ymin = 0, ymax = trace_lockdown), fill = "#000000", alpha = 0.2) +
        facet_grid(compartment ~ scenario, switch = "y", scales = "free") +
        scale_y_continuous(labels = axis_friendly, limits = c(0, NA)) +
        scale_x_date(date_breaks = "1 month", labels = axis_date) +
        theme(strip.background = element_blank(), strip.placement = "outside", 
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1), legend.position = "none") +
        labs(x = NULL, y = NULL);
}

plot_epi = function(d0, t, quant, ymd_start, ymd_truncate = "2050-01-01", colours = NULL, maxy = NA, which_region = "United Kingdom", exclude = NULL)
{
    # Copy data and process
    d = duplicate(d0[region %in% c(which_region, "All")])
    d[, scenario := factor(scenario, levels = unique(scenario))];
    d[, compartment := factor(compartment, levels = unique(compartment))];
    
    # Choose runs
    quant = union(quant, 0.5);
    mrun = t[scenario == "Base", sum(total), by = run][, which(rank(V1) == round(1 + (.N - 1) * 0.5))];
    qrun = t[scenario == "Base", sum(total), by = run][, which(rank(V1) %in% round(1 + (.N - 1) * quant))];
    d = d[run %in% qrun];
    
    # Merge intervention traces
    trace_school = d[compartment == "trace_school", .(run, t, trace_school = value - 1, scenario)]
    d = merge(d, trace_school, by = c("run", "t", "scenario"), all.x = T);
    d[, trace_school := trace_school * min(maxy, max(value), na.rm = T), by = .(compartment)];

    trace_intervention = d[compartment == "trace_intervention", .(run, t, trace_intervention = value - 1, scenario)]
    d = merge(d, trace_intervention, by = c("run", "t", "scenario"), all.x = T);
    d[, trace_intervention := trace_intervention *min(maxy, max(value), na.rm = T) * 0.75, by = .(compartment)];

    if (d[compartment == "trace_lockdown", .N > 0]) {
        trace_lockdown = d[run == mrun & compartment == "trace_lockdown", .(t, trace_lockdown = value - 1, scenario)]
        d = merge(d, trace_lockdown, by = c("t", "scenario"), all.x = T);
        d[, trace_lockdown := trace_lockdown * min(maxy, max(value), na.rm = T) * 0.5, by = .(compartment)];
    } else {
        d[, trace_lockdown := 0];
    }

    d = d[compartment %in% c("cases", "deaths", "beds_icu", "beds_nonicu")];

    # Give nice names
    d[compartment == "cases", compartment := "New cases"];
    d[compartment == "deaths", compartment := "Deaths"];
    d[compartment == "beds_icu", compartment := "ICU beds\nrequired"];
    d[compartment == "beds_nonicu", compartment := "Non-ICU beds\nrequired"];
    d[, compartment := factor(compartment, levels = unique(compartment))];

    d[, date := ymd(ymd_start) + t];

    # Plot
    plot = ggplot(d[region == which_region & date <= ymd(ymd_truncate) & !(scenario %in% exclude)]) +
        geom_line(aes(x = date, y = value, colour = scenario, group = run), size = 0.25, alpha = 0.35) +
        geom_ribbon(aes(x = date, ymin = 0, ymax = trace_school, group = run), fill = "#0000ff", alpha = 0.1/length(quant)) +
        geom_ribbon(aes(x = date, ymin = 0, ymax = trace_intervention, group = run), fill = "#ff0000", alpha = 0.2/length(quant)) +
        geom_ribbon(aes(x = date, ymin = 0, ymax = trace_lockdown, group = run), fill = "#000000", alpha = 0.15/length(quant)) +
        geom_line(data = d[region == which_region & date <= ymd(ymd_truncate) & run == mrun & !(scenario %in% exclude)], 
            aes(x = date, y = value, colour = scenario), size = 0.6) +
        facet_grid(compartment ~ scenario, switch = "y", scales = "free") +
        scale_y_continuous(labels = axis_friendly, limits = c(0, NA)) +
        scale_x_date(date_breaks = "1 month", labels = axis_date) +
        theme(strip.background = element_blank(), strip.placement = "outside", 
            axis.text.x = element_text(size = 5, angle = 45, hjust = 1), legend.position = "none") +
        labs(x = NULL, y = NULL);
    
    if (!is.null(colours)) {
        plot = plot + scale_colour_manual(values = colours)
    }
    if (!is.na(maxy)) {
        plot = plot + coord_cartesian(ylim = c(0, maxy));
    }
    return (plot)
}

# set theme
theme_set(theme_cowplot(font_size = 7, line_size = 0.25))

# load data
d1 =   reflow_dynamics(qread(paste0(covid_uk_path, "/1-dynamics.qs")));
t1 =     reflow_totals(qread(paste0(covid_uk_path, "/1-totals.qs")));
d2.1 = reflow_dynamics(qread(paste0(covid_uk_path, "/2.1-dynamics.qs")));
t2.1 =   reflow_totals(qread(paste0(covid_uk_path, "/2.1-totals.qs")));
d2.2 = reflow_dynamics(qread(paste0(covid_uk_path, "/2.2-dynamics.qs")));
t2.2 =   reflow_totals(qread(paste0(covid_uk_path, "/2.2-totals.qs")));
d3 =   reflow_dynamics(qread(paste0(covid_uk_path, "/3-dynamics.qs")));
t3 =     reflow_totals(qread(paste0(covid_uk_path, "/3-totals.qs")));
d4 =   reflow_dynamics(qread(paste0(covid_uk_path, "/4-dynamics.qs")));
t4 =     reflow_totals(qread(paste0(covid_uk_path, "/4-totals.qs")));
d6 =   reflow_dynamics(qread(paste0(covid_uk_path, "/6-dynamics.qs")));
t6 =     reflow_totals(qread(paste0(covid_uk_path, "/6-totals.qs")));
r0s =                  qread(paste0(covid_uk_path, "/5-dynamics.qs"));

# ANALYSIS 1 - 12 WEEK INTERVENTIONS
tb1 = make_table(d1)
pl1 = plot_table(tb1)
save_table(tb1, paste0(covid_uk_path, "/table-12week.csv"));
pl2 = plot_attackrate(t1)
pl3 = plot_epi(d1, t1, (0:10)/10, "2020-01-29", "2020-10-15")
f = plot_grid(pl1, pl2, pl3, ncol = 1, rel_heights = c(6, 6, 10), labels = c("a", "b", "c"), label_size = 9);
ggsave(paste0(covid_uk_path, "/full-1.pdf"), f, width = 20, height = 22, units = "cm", useDingbats = F);

gg_color_hue = function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
cols6 = gg_color_hue(6)
cols5 = gg_color_hue(5)
cols4 = gg_color_hue(4)

pla1 = plot_epi(d1[compartment != "deaths" & compartment != "beds_nonicu" & 
        scenario %in% c("Base", "School Closures", "Social Distancing")], t1, (0:10)/10, "2020-01-29", "2020-10-15", 
        colours = cols6[1:3]);
pla2 = plot_epi(d1[compartment != "deaths" & compartment != "beds_nonicu" & 
        !(scenario %in% c("Base", "School Closures", "Social Distancing"))], t1, (0:10)/10, "2020-01-29", "2020-10-15",
        colours = cols6[4:6]);
tb1[statistic == "Cases in peak week", statistic := "Cases in\npeak week"];
tb1[statistic == "Peak ICU beds required", statistic := "Peak ICU beds\nrequired"];
tb1[statistic == "Peak non-ICU beds required", statistic := "Peak non-ICU beds\nrequired"];
tb1[statistic == "Time to peak cases (weeks)", statistic := "Time to peak\ncases (weeks)"];
plb = plot_table(tb1[statistic != "Deaths in peak week"]) + theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 2, byrow = TRUE)) + labs(colour = NULL)
r0s1 = r0s[!(scenario %like% "Intensive") & !(scenario %like% "Lockdown")];
r0s1[, scenario := factor(scenario, levels = unique(scenario))];
plR = ggplot(r0s1) +
    geom_violin(aes(x = R0, y = scenario, fill = scenario), colour = NA) +
    geom_vline(xintercept = 1, size = 0.25, linetype = "44") +
    theme(legend.position = "none", axis.text.y = element_text(size = 6)) +
    labs(x = expression(R[0]), y = NULL) +
    scale_y_discrete(limits = rev(unique(r0s1$scenario))) +
    xlim(0, NA)
f = plot_grid(pla1, plb, pla2, plR, 
    nrow = 2, ncol = 2, rel_widths = c(3, 2), labels = c("a", "b", "", "c"), label_size = 9, align = "hv", axis = "bottom")
ggsave(paste0(covid_uk_path, "/fig-12week.pdf"), f, width = 20, height = 12, units = "cm", useDingbats = F);

# ANALYSIS 2 - TRIGGERS
d2.1[scenario != "Base", scenario := paste(scenario, "national")]
d2.2[scenario != "Base", scenario := paste(scenario, "local")]
d2 = rbind(d2.1, d2.2[scenario != "Base"])
t2 = rbind(t2.1, t2.2[scenario != "Base"])
d2[scenario == "Combination 0 shift local", scenario := "Local trigger"]
d2[scenario == "Combination 14 shift local", scenario := "Local trigger, +2 weeks"]
d2[scenario == "Combination 28 shift local", scenario := "Local trigger, +4 weeks"]
d2[scenario == "Combination 56 shift local", scenario := "Local trigger, +8 weeks"]
d2[scenario == "Combination 0 shift national", scenario := "National trigger"]
d2[scenario == "Combination 14 shift national", scenario := "National trigger, +2 weeks"]
d2[scenario == "Combination 28 shift national", scenario := "National trigger, +4 weeks"]
d2[scenario == "Combination 56 shift national", scenario := "National trigger, +8 weeks"]
d2[, scenario := factor(scenario, levels = c("Base", "Local trigger", "National trigger",
    "Local trigger, +2 weeks", "National trigger, +2 weeks",
    "Local trigger, +4 weeks", "National trigger, +4 weeks",
    "Local trigger, +8 weeks", "National trigger, +8 weeks"))]
d2 = d2[order(d2$scenario)]
tb2 = make_table(d2)
pl1 = plot_table(tb2)
save_table(tb2, paste0(covid_uk_path, "/table-triggers.csv"));
pl2 = plot_attackrate(t2)
pl3 = plot_epi(d2, t2, (0:10)/10, "2020-01-29")
f = plot_grid(pl1, pl2, pl3, ncol = 1, rel_heights = c(6, 6, 10), labels = c("a", "b", "c"), label_size = 9);
ggsave(paste0(covid_uk_path, "/COVID-UK/full-2.pdf"), f, width = 20, height = 22, units = "cm", useDingbats = F);

pla1 = plot_epi(d2[compartment != "deaths" & compartment != "beds_nonicu" & compartment != "beds_icu" &
        scenario %like% "Local"], t2, (0:10)/10, "2020-01-29", "2020-8-31", exclude = "Base");
pla2 = plot_epi(d2[compartment != "deaths" & compartment != "beds_nonicu" & compartment != "beds_icu" &
        scenario %like% "National"], t2, (0:10)/10, "2020-01-29", "2020-8-31", exclude = "Base");

tb2[statistic == "Cases in peak week", statistic := "Cases in\npeak week"];
tb2[statistic == "Peak ICU beds required", statistic := "Peak ICU beds\nrequired"];
tb2[statistic == "Peak non-ICU beds required", statistic := "Peak non-ICU beds\nrequired"];
tb2[statistic == "Time to peak cases (weeks)", statistic := "Time to peak\ncases (weeks)"];
plb = plot_table(tb2[statistic != "Deaths in peak week"], nrow = 3) + theme(legend.position = "bottom") +
    guides(colour = guide_legend(nrow = 5, byrow = F)) + labs(colour = NULL)

plc1 = plot_epi(d2[compartment != "deaths" & compartment != "beds_nonicu" & compartment != "cases" &
        scenario %like% "Base"], t2, (0:10)/10, "2020-01-29", "2020-6-30", which_region = "Cumbria") + 
    labs(title = "County 1") + theme(strip.text = element_blank(), axis.text.x = element_blank()) + xlim(ymd("2020-03-01"), NA)
plc2 = plot_epi(d2[compartment != "deaths" & compartment != "beds_nonicu" & compartment != "cases" &
        scenario %like% "Base"], t2, (0:10)/10, "2020-01-29", "2020-6-30", which_region = "West Midlands") + 
    labs(title = "County 2") + theme(strip.text = element_blank(), axis.text.x = element_blank()) + xlim(ymd("2020-03-01"), NA)
plc3 = plot_epi(d2[compartment != "deaths" & compartment != "beds_nonicu" & compartment != "cases" &
        scenario %like% "Base"], t2, (0:10)/10, "2020-01-29", "2020-6-30", which_region = "United Kingdom") + 
    labs(title = "UK") + theme(strip.text = element_blank()) + xlim(ymd("2020-03-01"), NA)

pla = plot_grid(pla1, pla2, nrow = 2, align = "hv", axis = "bottom")
plc = plot_grid(plc1, plc2, plc3, nrow = 3, align = "v", axis = "bottom", rel_heights = c(1, 1, 1.2))

f = plot_grid(pla, plb, plc,
    ncol = 3, labels = c("a", "b", "c"), label_size = 9, rel_widths = c(2, 1, .5))
ggsave(paste0(covid_uk_path, "/fig-triggers.pdf"), f, width = 20, height = 8, units = "cm", useDingbats = F);

# ANALYSIS 3 - LOCKDOWN
d3[scenario == "Intensive Interventions NA lockdown", scenario := "Intensive Interventions"];
d3[scenario == "Intensive Interventions 1000 lockdown", scenario := "Lockdown 1000-bed trigger"];
d3[scenario == "Intensive Interventions 2000 lockdown", scenario := "Lockdown 2000-bed trigger"];
d3[scenario == "Intensive Interventions 5000 lockdown", scenario := "Lockdown 5000-bed trigger"];
#d3[compartment == "subclinical"]$value = d3[compartment == "subclinical"]$value + d3[compartment == "cases"]$value

tb3 = make_table(d3)
pl1 = plot_table(tb3)
save_table(tb3, paste0(covid_uk_path, "/table-lockdown.csv"));
pl2 = plot_attackrate(t3)
pl3 = plot_epi(d3[scenario != "Base"], t3, (0:10)/10, "2020-01-29")
f = plot_grid(pl1, pl2, pl3, ncol = 1, rel_heights = c(6, 6, 10), labels = c("a", "b", "c"), label_size = 9);
ggsave(paste0(covid_uk_path, "/full-3.pdf"), f, width = 20, height = 22, units = "cm", useDingbats = F);

tb3[statistic == "Cases in peak week", statistic := "Cases in\npeak week"];
tb3[statistic == "Peak ICU beds required", statistic := "Peak ICU beds\nrequired"];
tb3[statistic == "Peak non-ICU beds required", statistic := "Peak non-ICU beds\nrequired"];
tb3[statistic == "Time to peak cases (weeks)", statistic := "Time to peak\ncases (weeks)"];
tb3[statistic == "Proportion of time spent in lockdown", statistic := "Proportion of time\nspent in lockdown"];
plb = plot_table(tb3[statistic != "Deaths in peak week" & statistic != "Time to peak\ncases (weeks)" & statistic != "Total subclinical" & scenario != "Base"]) + 
    theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 3, byrow = TRUE)) + labs(colour = NULL) +
    scale_colour_manual(values = cols5[2:5])

d3 = d3[scenario != "Base"]
d3 = d3[scenario != "Intensive Interventions" | t < 425]
d3 = d3[scenario == "Intensive Interventions" | t < 575]
pla1 = plot_epi(d3[compartment != "deaths" & compartment != "beds_nonicu" & compartment != "cases" & 
        scenario %in% c("Intensive Interventions", "Lockdown 1000-bed trigger")], t3, (0:10)/10, "2020-01-29", 
        colours = cols5[2:3], maxy = 26000) + geom_hline(data = data.table(scenario = c("Intensive Interventions", "Intensive Interventions", "Lockdown 1000-bed trigger", "Lockdown 1000-bed trigger"), 
            beds = c(4562, 9124, 4562, 9124), lt = c("a", "b", "a", "b")),
            aes(yintercept = beds, linetype = lt), size = 0.25)
pla2 = plot_epi(d3[compartment != "deaths" & compartment != "beds_nonicu" & compartment != "cases" & 
        scenario %in% c("Lockdown 2000-bed trigger", "Lockdown 5000-bed trigger")], t3, (0:10)/10, "2020-01-29",
        colours = cols5[4:5], maxy = 26000) + geom_hline(data = data.table(scenario = c("Lockdown 2000-bed trigger", "Lockdown 2000-bed trigger", "Lockdown 5000-bed trigger", "Lockdown 5000-bed trigger"),
            beds = c(4562, 9124, 4562, 9124), lt = c("a", "b", "a", "b")),
            aes(yintercept = beds, linetype = lt), size = 0.25)

r0s2 = r0s[scenario %like% "Base" | scenario %like% "Intensive" | scenario %like% "Lockdown"]
r0s2[scenario == "Intensive, schools open", scenario := "Intensive,\nschools open"]
r0s2[scenario == "Intensive, schools closed", scenario := "Intensive,\nschools closed"]
r0s2[, scenario := factor(scenario, levels = unique(scenario))]
plR = ggplot(r0s2) +
    geom_violin(aes(x = R0, y = scenario, fill = scenario), colour = NA) +
    geom_vline(xintercept = 1, size = 0.25, linetype = "44") +
    theme(legend.position = "none") +
    labs(x = expression(R[0]), y = NULL) +
    xlim(0, NA) +
    scale_y_discrete(limits = rev(unique(r0s2$scenario))) +
    scale_fill_manual(values = c(cols5[1], cols5[2], cols5[2], "#bbbbbb"))
f = plot_grid(pla1, plb, pla2, plR, 
    nrow = 2, ncol = 2, rel_widths = c(3, 2), labels = c("a", "b", "", "c"), label_size = 9, align = "hv", axis = "bottom")
ggsave(paste0(covid_uk_path, "/fig-lockdown.pdf"), f, width = 20, height = 12, units = "cm", useDingbats = F);

# ANALYSES 4,6 - GRANDPARENTS AND SPORTS/LEISURE
# Grandparents
d4[scenario == "Intensive", scenario := "Background"]
d4[scenario == "Intensive + School", scenario := "School closure"]
d4[scenario == "Intensive + School + G20", scenario := "School closure, 20% care by elderly"]
d4[scenario == "Intensive + School + G50", scenario := "School closure, 50% care by elderly"]
d4[scenario == "Intensive + School + G100", scenario := "School closure, 100% care by elderly"]
tb4 = make_table(d4)
pl1 = plot_table(tb4)
save_table(tb4, paste0(covid_uk_path, "/table-grandparents.csv"));
pl2 = plot_attackrate(t4)
pl3 = plot_epi(d4, t4, (0:10)/10, "2020-01-29", "2020-12-31")
f = plot_grid(pl1, pl2, pl3, ncol = 1, rel_heights = c(6, 6, 10), labels = c("a", "b", "c"), label_size = 9);
ggsave(paste0(covid_uk_path, "/full-4.pdf"), f, width = 20, height = 22, units = "cm", useDingbats = F);

tb4[statistic == "Cases in peak week", statistic := "Cases in\npeak week"];
tb4[statistic == "Peak ICU beds required", statistic := "Peak ICU beds\nrequired"];
tb4[statistic == "Peak non-ICU beds required", statistic := "Peak non-ICU beds\nrequired"];
tb4[statistic == "Time to peak cases (weeks)", statistic := "Time to peak\ncases (weeks)"];
plb = plot_table(tb4[scenario != "Base" & statistic != "Deaths in peak week"]) + theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 3, byrow = TRUE)) + labs(colour = NULL)

# Sports and leisure
d6[scenario == "Background", scenario := "Background"]
d6[scenario == "Background + 0% Sports", scenario := "Spectator sports banned"]
d6[scenario == "Background + 25% Leisure", scenario := "Leisure reduced by 75%"]
tb6 = make_table(d6)
pl1 = plot_table(tb6)
save_table(tb6, paste0(covid_uk_path, "/table-sports.csv"));
pl2 = plot_attackrate(t6)
pl3 = plot_epi(d6, t6, (0:10)/10, "2020-01-29")
f = plot_grid(pl1, pl2, pl3, ncol = 1, rel_heights = c(6, 6, 10), labels = c("a", "b", "c"), label_size = 9);
ggsave(paste0(covid_uk_path, "/full-sports.pdf"), f, width = 20, height = 22, units = "cm", useDingbats = F);

tb6[statistic == "Cases in peak week", statistic := "Cases in\npeak week"];
tb6[statistic == "Peak ICU beds required", statistic := "Peak ICU beds\nrequired"];
tb6[statistic == "Peak non-ICU beds required", statistic := "Peak non-ICU beds\nrequired"];
tb6[statistic == "Time to peak cases (weeks)", statistic := "Time to peak\ncases (weeks)"];
pla = plot_table(tb6[scenario != "Base" & statistic != "Deaths in peak week"]) + theme(legend.position = "bottom") + 
    guides(colour = guide_legend(nrow = 3, byrow = TRUE)) + labs(colour = NULL)
f = plot_grid(pla, plb, nrow = 2, labels = c("a", "b"), label_size = 9, align = "hv", axis = "bottom");
ggsave(paste0(covid_uk_path, "/fig-misc.pdf"), f, width = 9, height = 12, units = "cm", useDingbats = F);


# NINGBO EST.
x = rbeta(100000, 6, 140)/rbeta(100000, 126, 1875)
cm_mean_hdi(x)
