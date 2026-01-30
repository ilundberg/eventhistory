viz_survreg <- function(x) {
  if (!(x$dist %in% c("exponential","weibull","lognormal"))) {
    stop("viz_survreg can currently only handle survreg models with dist = exponential, weibull, or lognormal")
  }
  
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
  } else if (x$dist == "lognormal") {
    if (length(unique(x$linear.predictors)) == 1) {
      tibble(
        meanlog = unique(x$linear.predictors),
        sdlog = x$scale,
        time = x$y[,"time"]
      ) |>
        mutate(
          survival = plnorm(time, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE),
          density = dlnorm(time, meanlog = meanlog, sdlog = sdlog),
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
        meanlog = c(
          x$linear.predictors |> min(),
          x$linear.predictors |> median(),
          x$linear.predictors |> max()
        ),
        sdlog = x$scale
      ) |>
        mutate(time = map(group, \(g) x$y[,"time"])) |>
        unnest(cols = "time") |>
        mutate(
          survival = plnorm(time, meanlog = meanlog, sdlog = sdlog, lower.tail = FALSE),
          density = dlnorm(time, meanlog = meanlog, sdlog = sdlog),
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
