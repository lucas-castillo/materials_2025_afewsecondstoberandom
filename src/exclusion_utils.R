out_of_range <- function(
    sequence, 
    hard_lb, 
    hard_ub, 
    soft_lb, 
    soft_ub, 
    threshold){
  sequence <- sequence[!is.na(sequence)]
  a <- max(sequence) < hard_lb # no items above hard lower bound, OR
  b <- min(sequence) > hard_ub # no items below hard upper bound, OR
  c <- mean( # proportion items outside soft bounds higher than threshold
    sequence < soft_lb | sequence > soft_ub
    ) > threshold
  a | b | c
}
slow_pace_harsh <- function(sequence, expectation, starts){
  expectation <- expectation[1]
  isSparse <- F
  threshold <- expectation * .8 * .5
  for (s in starts){
    if (s < 90){
      t <- starts[starts >= s & starts < (s + 30)]
      if (length(t) < threshold) isSparse = T
    }
  }
  return(isSparse)
}


slow_pace_lenient <- function(sequence, expectation, starts){
  expectation <- expectation[1]
  windows <- c()
  threshold <- 80 * .8
  pace <- 60 / expectation
  paceTick <- pace / .8
  index <- 1
  total <- length(starts)
  while (index < total){
    s <- starts[index]
    max <- index + threshold - 1
    candidateWindow <- starts[index:max]
    nJumps <- length(candidateWindow) - 1
    goodPace <- max((candidateWindow - candidateWindow[1])[2:length(candidateWindow)]) < paceTick * nJumps
    if (max <= total & goodPace){
      while (goodPace & max < total){
        max <- max+1
        candidateWindow <- starts[index:max]
        nJumps <- length(candidateWindow) - 1
        goodPace <- max((candidateWindow - candidateWindow[1])[2:length(candidateWindow)]) < paceTick * nJumps
      }
      if (max != total){max <- max-1}
      min <- index
      windows <- c(windows, min:max)
      if ((max + expectation) < total){index <- max} else{index <- total}
    } else {
      if ((index + threshold - 1) < total){
        index <- index + 1  
      } else{
        index <- total
      }
    }
  }
  r <- rep(F, length(sequence))
  r[windows] <- T
  return(r)
}

slow_pace_p <- function(sequence, expectation, threshold){
  if (length(sequence) < (expectation[1] * threshold)){
    return(F)
  } else{
    return(T)
  }
}