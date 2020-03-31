# covidm_plot.R
# functions to plot model results

# Plot priors and posterior distributions from a model fit
cm_plot_posterior = function(fit)
{
    if (!any("cm.fit" %in% class(fit))) {
        stop("cm_plot_posterior requires a cm.fit object.");
    }
    
    # get prior distributions
    priors = NULL;
    for (pr in 1:length(fit$priors)) {
        ep = cm_evaluate_distribution(fit$priors[[pr]]);
        setDT(ep);
        ep[, p := p / max(p)]
        ep[, parameter := names(fit$priors)[pr]];
        priors = rbind(priors, ep);
    }
    posterior = melt(fit$posterior, id.vars = NULL, measure.vars = 5:ncol(fit$posterior), variable.name = "parameter");
    
    # put in correct order
    priors[, parameter := factor(parameter, levels = names(fit$priors))]
    posterior[, parameter := factor(parameter, levels = names(fit$priors))]
    
    # plot
    ggplot() +
        geom_line(data = priors, aes(x = x, y = p), colour = "grey", linetype = "22", size = 0.25) +
        geom_histogram(data = posterior, aes(value, y = stat(ndensity), fill = parameter), alpha = 0.75, bins = 30) +
        facet_wrap(~parameter, scales = "free") +
        labs(x = NULL, y = NULL) +
        theme(legend.position = "none", strip.background = element_blank())
}

# Plot pairwise posterior distributions from a model fit
cm_plot_pairwise = function(fit)
{
    if (!any("cm.fit" %in% class(fit))) {
        stop("cm_plot_pairwise requires a cm.fit object.");
    }
    
    # get prior distributions
    priors = NULL;
    for (pr in 1:length(fit$priors)) {
        ep = cm_evaluate_distribution(fit$priors[[pr]]);
        setDT(ep);
        ep[, p := p / max(p)]
        ep[, parameter := names(fit$priors)[pr]];
        priors = rbind(priors, ep);
    }
    posterior = melt(fit$posterior, id.vars = NULL, measure.vars = 5:ncol(fit$posterior), variable.name = "parameter");
    
    # put in correct order
    priors[, parameter := factor(parameter, levels = names(fit$priors))];
    posterior[, parameter := factor(parameter, levels = names(fit$priors))];
    
    # assemble plots
    plotlist = list();
    i = 1;
    for (r in 1:length(fit$priors)) {
        for (c in 1:r) {
            rp = names(fit$priors)[r];
            cp = names(fit$priors)[c];
            if (r == c) { # posterior element on its own
                cat("\\");
                plot = ggplot() +
                    geom_line(data = priors[parameter == rp], aes(x = x, y = p), colour = "grey", linetype = "22", size = 0.25) +
                    geom_histogram(data = posterior[parameter == rp], aes(value, y = stat(ndensity)), fill = "#88bbff", bins = 30) +
                    labs(x = NULL, y = NULL, title = rp)
            } else { # pairwise plot
                cat("x");
                data = cbind(posterior[parameter == cp, .(x = value)], posterior[parameter == rp, .(y = value)]);
                plot = ggplot(data) +
                    geom_density2d(aes(x, y), colour = "#ffbb88") +
                    labs(x = NULL, y = NULL)
            }
            plotlist[[i]] = plot;
            i = i + 1;
        }
        for (c in seq_len(length(fit$priors) - r) + r) {
            cat(".");
            plotlist[[i]] = NULL;
            i = i + 1;
        }
        cat("\n");
    }
    
    plot_grid(plotlist = plotlist, nrow = length(fit$priors), ncol = length(fit$priors))
}
