## Function to compare parametric survival to Kaplan-Meier survival

km_compare <- function(x, subgroup_rows = 1:nrow(x$y)) {
  
  if (!(x$dist %in% c("exponential","weibull","lognormal"))) {
    stop("km_compare can currently only handle survreg models with dist = exponential, weibull, or lognormal")
  }
  
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
  } else if (x$dist == "lognormal") {
    parametric <- tibble(
      meanlog = x$linear.predictors[subgroup_rows],
      sdlog = x$scale,
      time = time_values
    ) |>
      unnest(cols = "time") |>
      mutate(surv = plnorm(time, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE)) |>      # Now take population mean survival at each time
      group_by(time) |>
      summarize(surv = mean(surv)) |>
      mutate(method = "Lognormal Model")
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
