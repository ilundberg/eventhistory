
library(tidyverse)
library(survival)

# - Simulation with constant hazard and unmeasured heterogeneity ()
# - See that Exponential does not match Kaplan-Meier well (km_compare)
# - Try Weibull model (survreg)
# - See that matches Kaplan-Meier well (km_compare)
# - See that can reject null of Exponential (summary)
# - See that hazard is not constant (viz_survreg)


simulated <- tibble(id = 1:1e3) |>
  mutate(
    frailty = ifelse(
      id %% 2 == 0, 1, 10
    ),
    t = rexp(n(), rate = frailty),
    event = 1
  )

weibull_model <- survreg(
  Surv(t, event = event) ~ 1,
  data = simulated,
  dist = "weibull"
)

predict_curves <- function(x, newdata = NULL, time_values = NULL) {
  if (is.null(newdata)) {
    if (exists(x$call$data)) {
      newdata <- get(x$call$data)
      message(paste0("Using the object '",x$call$data,"' from your environment as 'newdata'."))
    }
  }
  # Warn if data contains columns with names that will be overwritten
  for (problem_col_name in c("time","shape","scale","survival","density","hazard")) {
    if (any(colnames(newdata) == problem_col_name)) {
      message(paste0("'newdata' contains a column '",problem_col_name,"' that will be overwritten."))
    }
  }

  if (is.null(time_values)) {
    first_event <- x$y[,"time"] |> min()
    last_event <- x$y[,"time"] |> as.numeric() |> max()
    time_values <- seq(first_event, last_event, length.out = 100)
    message("By default, making predictions at 100 evenly spaced values over the range of times observed in the data used to estimate the model.")
  }
  if (x$dist != "weibull") {
    stop("This function currently only handles 'survreg' objects with dist = 'weibull'.")
  } else if (x$dist == "weibull") {
    predicted_scale <- exp(predict(x, newdata = newdata, type = "linear"))
    predicted_shape <- 1 / x$scale
    
    predicted <- newdata |>
      mutate(
        shape = predicted_shape,
        scale = predicted_scale,
        time = map(.x = scale, .f = \(x) time_values),
        survival = map2(
          .x = shape, 
          .y = scale,
          .f = \(x,y) {
            pweibull(time_values, shape = x, scale = y)
          }
        ),
        density = map2(
          .x = shape, 
          .y = scale,
          .f = \(x,y) {
            dweibull(time_values, shape = x, scale = y)
          }
        )
      ) |>
      unnest(cols = c("time","survival","density")) |>
      mutate(hazard = density / survival)
    
    return(predicted)
  }
}

predict_curves(weibull_model, newdata = tibble(id = 1), time_values = seq(.1,25,.1)) |>
  ggplot(aes(x = time, y = hazard)) +
  geom_line()


viz_survreg <- function(x) {
  if (x$dist == "exponential") {
    if (length(unique(x$linear.predictors)) == 1) {
      tibble(
        scale = exp(unique(x$linear.predictors)),
        time = x$y[,"time"]
      ) |>
        mutate(
          survival = pexp(time, rate = 1 / scale, lower.tail = FALSE),
          density = dexp(time, rate = 1 / scale),
          hazard = density / survival
        ) |>
        select(time, survival, hazard) |>
        pivot_longer(cols = -time) |>
        ggplot(
          aes(x = time, y = value)
        )+
        geom_line() +
        facet_wrap(~name, scales = "free_y", ncol = 1)
    } else {
      tibble(
        group = paste(c("Min","Median","Max"),"of X * Beta"),
        scale = c(
          x$linear.predictors |> exp() |> min(),
          x$linear.predictors |> exp() |> median(),
          x$linear.predictors |> exp() |> max()
        )
      ) |>
        mutate(time = map(group, \(g) x$y |> as.numeric())) |>
        unnest(cols = "time") |>
        mutate(
          survival = pexp(time, rate = 1 / scale, lower.tail = FALSE),
          density = dexp(time, rate = 1 / scale),
          hazard = density / survival
        ) |>
        select(group, time, survival, hazard) |>
        pivot_longer(cols = -c("time","group")) |>
        ggplot(
          aes(x = time, y = value, linetype = group)
        )+
        geom_line() +
        facet_wrap(~name, scales = "free_y", ncol = 1)
    }
  } else if (x$dist == "weibull") {
    if (length(unique(x$linear.predictors)) == 1) {
      tibble(
        scale = exp(unique(x$linear.predictors)),
        shape = 1 / x$scale,
        time = x$y[,"time"]
      ) |>
        mutate(
          survival = pweibull(time, shape, scale, lower.tail = FALSE),
          density = dweibull(time, shape, scale),
          hazard = density / survival
        ) |>
        select(time, survival, hazard) |>
        pivot_longer(cols = -time) |>
        ggplot(
          aes(x = time, y = value)
        )+
        geom_line() +
        facet_wrap(~name, scales = "free_y", ncol = 1)
    } else {
      tibble(
        group = paste(c("Min","Median","Max"),"of X * Beta"),
        scale = c(
          x$linear.predictors |> exp() |> min(),
          x$linear.predictors |> exp() |> median(),
          x$linear.predictors |> exp() |> max()
        ),
        shape = 1 / x$scale
      ) |>
        mutate(time = map(group, \(g) x$y[,"time"])) |>
        unnest(cols = "time") |>
        mutate(
          survival = pweibull(time, shape, scale, lower.tail = FALSE),
          density = dweibull(time, shape, scale),
          hazard = density / survival
        ) |>
        select(group, time, survival, hazard) |>
        pivot_longer(cols = -c("time","group")) |>
        ggplot(
          aes(x = time, y = value, linetype = group)
        )+
        geom_line() +
        facet_wrap(~name, scales = "free_y", ncol = 1)
    }
  }
}

weibull_model |>
  viz_survreg()


weibull_model_correct <- survreg(
  Surv(t, event = event) ~ log(frailty),
  data = simulated,
  dist = "weibull"
)

weibull_model_correct |>
  viz_survreg()


expo_model_correct <- survreg(
  Surv(t, event = event) ~ log(frailty),
  data = simulated,
  dist = "exponential"
)

expo_model_correct |>
  viz_survreg()


## Function to compare parametric survival to Kaplan-Meier survival

km_compare <- function(x, subgroup_rows = 1:nrow(x$y)) {
  
  # Fit Kaplan-Meier to compare
  survfit.out <- survfit(x$y[subgroup_rows,] ~ 1)
  km <- tibble(
    time = survfit.out$time,
    surv = survfit.out$surv,
    method = "Kaplan-Meier"
  )
  
  # Now create the same survival points with the parametric model
  time_values <- vector("list", length = length(subgroup_rows))
  for (i in 1:length(time_values)) {
    time_values[[i]] <- km$time
  }
  
  if (x$dist == "exponential") {
    parametric <- tibble(
      rate = 1 / exp(x$linear.predictors)[subgroup_rows],
      time = time_values
    ) |>
      unnest(cols = "time") |>
      mutate(surv = pexp(time, rate = rate, lower.tail = FALSE)) |>
      # Now take population mean survival at each time
      group_by(time) |>
      summarize(surv = mean(surv)) |>
      mutate(method = "Exponential Model")
  } else if (x$dist == "weibull") {
    parametric <- tibble(
      scale = exp(x$linear.predictors)[subgroup_rows],
      shape = 1 / x$scale,
      time = time_values
    ) |>
      unnest(cols = "time") |>
      mutate(surv = pweibull(time, scale = scale, shape = shape, lower.tail = FALSE)) |>      # Now take population mean survival at each time
      group_by(time) |>
      summarize(surv = mean(surv)) |>
      mutate(method = "Weibull Model")
  }
  
  message(
    "Note: Graph compares marginal survival curves: E_X(P(T > t | X)). If the model and Kaplan-Meier are similar for this marginal survival curve, it does not necessarily imply that they are similar within a subgroup on P(T > t | X = x). If a subgroup X = x is well-populated, you can test the model within the subgroup by re-running this comparison on a model estimated within the subgroup."
  )
  
  parametric |>
    bind_rows(km) |>
    ggplot(aes(x = time, y = surv, color = method, linetype = method)) +
    geom_line() +
    labs(
      x = "Time",
      y = "Survival",
      color = "Method",
      linetype = "Method"
    )
}

km_compare(expo_model_correct)

km_compare(weibull_model_correct)


expo_wrong <- survreg(
  Surv(t, event) ~ 1, 
  data = simulated,
  dist = "exponential"
)
km_compare(expo_wrong)
km_compare(expo_wrong, subgroup_rows = which(simulated$frailty == 1))
km_compare(expo_wrong, subgroup_rows = which(simulated$frailty == 10))


weibull_model <- survreg(
  Surv(t, event) ~ 1, 
  data = simulated,
  dist = "weibull"
)
km_compare(weibull_model)
km_compare(weibull_model, subgroup_rows = which(simulated$frailty == 1))
km_compare(weibull_model, subgroup_rows = which(simulated$frailty == 10))

# You can reject the null of a constant hazard,
# even though the DGP was a constant hazard in two subgroups
summary(weibull_model)

1 / weibull_model$scale
weibull_model |>
  viz_survreg()

