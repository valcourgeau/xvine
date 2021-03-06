
fit_markov_svines <- function(data, k.markov, ...){
  assertthat::see_if(max(data) <= 1)
  assertthat::see_if(min(data) >= 0)

  n_cores <- parallel::detectCores()
  return(svines::svine(data = data, p = k.markov, cores = n_cores, ...))
}

fit_markov_rvines <- function(data, k.markov, ...){
  assertthat::see_if(max(data) <= 1)
  assertthat::see_if(min(data) >= 0)

  n_cores <- parallel::detectCores()
  d <- ncol(data)
  dt_stack <- build_stack(data, k.markov)
  vine <- rvinecopulib::vine(
    dt_stack, cores = n_cores,
    margins_controls=list(
      xmin=rep(0, d),
      xmax=rep(1, d)
    ),
    ...
  )
  vine$d <- ncol(data)
  return(vine)
}

fit_markov_conditional_rvines <- function(data, k.markov, col_source, u0, ...){
  assertthat::see_if(max(data) <= 1)
  assertthat::see_if(min(data) >= 0)
  d <- dim(data)[2]
  n_cores <- parallel::detectCores()

  # above
  dt_above_stack <- build_above_conditional_stack(
    data = data, k = k.markov, col = col_source, u0 = u0
  )
  vine_above <- rvinecopulib::vine(
    dt_above_stack, cores = 1,
    margins_controls=list(
      xmin=rep(0, d*k.markov), #pmax(1e-5, apply(dt_above_stack, MARGIN = 2, min)),
      xmax=rep(1, d*k.markov) #pmin(1.0, apply(dt_above_stack, MARGIN = 2, max))
    ),
    ...
  )
  vine_above$d <- d

  # below
  dt_below_stack <- build_below_conditional_stack(
    data = data, k = k.markov, col = col_source, u0 = u0
  )
  vine_below <- rvinecopulib::vine(
    dt_below_stack, cores = n_cores-1,
    margins_controls=list(
      xmin=rep(0, d*k.markov), #pmax(0, apply(dt_below_stack, MARGIN = 2, min)),
      xmax=rep(1, d*k.markov) #pmin(1, apply(dt_below_stack, MARGIN = 2, max))
    ),
    ...
  )
  vine_below$d <- d

  return(
    structure(
      list(
        'above'=vine_above,
        'below'=vine_below,
        'd'= dim(data)[2],
        'col'=col,
        'u0'=u0
      ),
      class = 'cond_vine'
    )
  )
}

fit_markov_timewise_rvines <- function(data, k.markov, ...){
  # fits a vine to (x_t, t_{t+s}) for each 1 <= s <= k.markov
  assertthat::see_if(max(data) <= 1)
  assertthat::see_if(min(data) >= 0)
  d <- dim(data)[2]

  n_cores <- parallel::detectCores()
  dt_timewise_stack <- build_timewise_stack(
    data = data, k = k.markov
  )
  vine_timewise <- lapply(
    dt_timewise_stack,
    function(dt_stack){
      v <- rvinecopulib::vine(
        dt_stack, cores = n_cores,
        margins_controls=list(
          xmin=rep(0, d),
          xmax=rep(1, d)
        ), ...
      )
      return(v)
    }
  )

  return(
    structure(
      list(
        'timewise'=vine_timewise,
        'd'=d
      ),
      class = 'timewise_vine'
    )
  )
}

fit_markov_conditional_timewise_rvines <- function(data, k.markov, col_source, u0, ...){
  # fits a vine to (x_t, t_{t+s}) for each 1 <= s <= k.markov
  assertthat::see_if(max(data) <= 1)
  assertthat::see_if(min(data) >= 0)
  d <- dim(data)[2]

  n_cores <- parallel::detectCores()
  dt_above_timewise_stack <- build_above_conditional_timewise_stack(
    data = data, k = k.markov, col = col_source, u0 = u0
  )
  dt_below_timewise_stack <- build_below_conditional_timewise_stack(
    data = data, k = k.markov, col = col_source, u0 = u0
  )

  vine_above_timewise <- lapply(
    dt_above_timewise_stack,
    function(dt_stack){
      v <- rvinecopulib::vine(
        dt_stack, cores = n_cores,
        margins_controls=list(
          xmin=rep(0, d),
          xmax=rep(1, d)
        ), ...
      )
      return(v)
    }
  )

  vine_below_timewise <- lapply(
    dt_below_timewise_stack,
    function(dt_stack){
      v <- rvinecopulib::vine(
        dt_stack, cores = n_cores,
        margins_controls=list(
          xmin=rep(0, d),
          xmax=rep(1, d)
        ), ...
      )
      return(v)
    }
  )

  vine_above_timewise <- structure(
    list(
      'timewise'=vine_above_timewise,
      'd'=d
    ),
    class = 'timewise_vine'
  )

  vine_below_timewise <- structure(
    list(
      'timewise'=vine_below_timewise,
      'd'=d
    ),
    class = 'timewise_vine'
  )

  return(
    structure(
      list(
        'above_timewise'=vine_above_timewise,
        'below_timewise'=vine_below_timewise,
        'd'=d,
        'u0'=u0,
        'col'=col
      ),
      class = 'cond_timewise_vine'
    )
  )
}
