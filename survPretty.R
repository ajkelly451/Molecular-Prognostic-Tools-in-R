require(tidyverse)
require(Rmisc)
require(grDevices)
require(rlist)
require(ggplot2)
require(broom)
require(survival)

## Annotation for the survfit data. The default is 60 months = 5 years.
sf_ann <- function(sf, tm=60) {
  # sf is a survfit object, tm is the time you want to summarize at
  cl <- (sf %>% summary(time=tm))[c('lower','surv','upper')] %>% 
    data.frame() # Summarize survfit object 
  if(dim(sf) > dim(cl)[1]) cl <- rbind(cl, 0)
  cl <- cl %>% apply(., 1, function(x) {
    paste0(round(100*x[2]), '% (95% CI, ', round(100*x[1]), 
           '%-', round(100*x[3]), '%)')}) %>% t() %>% data.frame() # Return text reguarding the upper and lower around surv.
  return(cl)
}

survPretty <- function(form = NULL, data = NULL, pos = c(87.5, .875), cx = 3.5, xlab = "Months", 
                       ylab = "Overall Survival Rate", lbs = list(NA), one.t = FALSE, col = NULL, 
                       sup.lbl = NULL, sup.ann = F, pretty_sf = T) {
  # form must follow Surv(t, s) ~ vars
  # You must have either form or data, and data only works if it is a previous output of the survPretty function. 
  # 
  
  # Run a coxph fun to acquire logrank data for text
  # and Survfit for plotting data
  if (class(lbs) != 'list') lbs <- list(lbs) # Formatting for consistency later on
  if (!is.null(form)) {
    if (!is.null(data)) {
      # Run cox model using your formula if data given
      cox <- summary(coxph(form, data = data))
      sfit <- survfit(form, data = data)
    } else {
      # Run cox model if data isn't given
      cox <- summary(coxph(form))
      sfit <- survfit(form)
    }
  } else {
    # Extract coxph and survfit models if only data is given (rare)
    cox <- data$coxph
    sf <- data$survfit
  }
  sf <-  sfit %>% tidy() # Clean up the survfit output! 
  p <- cox$sctest[3]; hr <- cox$coef[2] # Extract logrank p and extimated HR (hazard ratio)
  if (one.t == TRUE) {
    p <- cox$sctest[3]/2 # If you are only testing for one tailed tests, return one-tailed p-value
    ## Caution! Right now this only works if you know the HR is going the right way, check yourself! 
    print('I see you entered that you want a 1-tailed test. Make sure your HR is going the direction you expected!')
  }
  fit <- sf[, -(6:8)] # Selects the relevant columns
  # Need to do some cleaning for the names of the columns in fit, so ggplot plays nicely. 
  if (!is.na(lbs[[1]][1])) {
    a <- list.rbind(strsplit(fit$strata, ', '))
    for (i in seq_along(lbs)) {
      b <- factor(a[, i])
      levels(b) <- lbs[[i]]
      fit[[paste0('strata', i)]] <- b
    }
  }
  
  # Check to make sure survival times start at 0 for each stratum. 
  # I am unsure of why it happens, but sometimes it is truncated. 
  d <- c(1, which(fit$n.risk == 1) + 1)
  d <- d[-length(d)]
  for (i in sort(d, decreasing = TRUE)) {
    if (fit$time[i] != 0) {
      add <- data.frame(time = 0, n.risk = fit$n.risk[i], n.event = 0, n.censor = 0, 
                 estimate = 1, strata = fit$strata[i])
      if (!is.na(lbs[[1]][1])) { 
        add$strata1 <- fit$strata1[i]
        if (length(lbs) > 1) add$strata2 <- fit$strata2[i]
      }
      if (i == 1) fit <- rbind(add, fit)
      else fit <- rbind(fit[1:(i - 1), ], add, fit[i:nrow(fit), ])
    }
  }
  if (is.na(lbs[[1]][1])) fit$strata1 <- fit$strata
  
  # Plot the results -- general ggplot part. I've used theme_bw because I like it. Should be able 
  # to use whatever you want by adding "+ theme_?()" at the end of your call of this. 
  pl <- ggplot(fit, aes(time, estimate, group = strata)) + theme_bw() + ylim(0, 1) + 
    theme(legend.position = "top") + xlab(xlab) + ylab(ylab)
  pts <- subset(fit, n.censor > 0) # this is our censored data. We insert it into the plot object with below code.  
  if (length(lbs) > 1) { plt <- pl + geom_step(aes(linetype = strata1)) + aes(color = strata2) + 
    scale_color_manual(values = col) + geom_point(data = pts, shape = 3, show.legend = F)} else {
        plt <- pl + geom_step(aes(linetype = strata1)) + geom_point(data = pts, shape = 3, show.legend = F) }
  # You can suppress this annotation, but it prints by default (info about the coxph fit)
  if (sup.ann == F) { plt <- plt + # geom_point(data = pts, shape = 3) + # guides(alpha = FALSE) + 
    annotate("text", pos[1], pos[2], 
             label = paste0('HR=', round(hr, 2), '\np=', formatC(p, digits = 2, format = "g")), cex = cx) }
  
  #Some general tidying of the legends and annotations
  plt$labels$linetype <- ""; plt$labels$colour <- ""
  if (is.na(lbs[[1]][1])) lapply(plt$layers, function(Y) Y$show.legend <- FALSE)
  if (!is.null(sup.lbl)) plt$layers[[sup.lbl]]$show.legend <- FALSE
  if (pretty_sf == T) { sf.pretty <- sf_ann(sfit)
   if (length(lbs) > 1) nm <- unlist(lapply(lbs[[1]], function(x) paste0(as.character(x), lbs[[2]]))) else 
     nm <- c('Low_expr', 'High_expr')
   names(sf.pretty) <- nm
   print(sf.pretty)
  }
  return(plt)
  
}




