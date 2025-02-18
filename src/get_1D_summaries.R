source("src/compute_all_rand_measures.R")
get_Summaries <- function(v, ub=NULL, lb=NULL, drop.r=T, is.circular=FALSE){
  unlist(lapply(get_Summaries_full(v, ub, lb, drop.r, is.circular), FUN = \(x){mean(x, na.rm=T)}))
}

get_Summaries_full <- function(v, ub=NULL, lb=NULL, drop.r=T, is.circular=FALSE){
  v <- round(v)
  # v <- v[!is.na(v)]
  r <- c(FALSE, v[2:length(v)] == v[1:(length(v)-1)])
  if (sum(!r, na.rm=T) == 1){
    return(list(
      "R"=1,
      "A"=0,
      "TPF"=NA,
      "D"=NA
    ))
  }


  if (is.null(ub)){
    ub <- quantile(v, 1/2 + 1/8, na.rm=T)
  }
  if (is.null(lb)){
    lb <- quantile(v, 1/2 - 1/8, na.rm=T)
  }
  v2 <- v[!r]
  n2 <- length(v2)
  
  r[1] <- NA
  n <- length(v)
  if (!is.circular){
    D <- abs(v[2:n] - v[1:(n-1)])
  } else{
    opts <- length(min(v):max(v))
    D1 <- c(abs(opts + v[2:length(v)] - (v[1:length(v) - 1])) %% opts)
    D2 <- c(abs(v[2:length(v)] - (opts + v[1:length(v) - 1])) %% opts)
    D <- sapply(1:length(D1), \(i){min(D1[i], D2[i])})
  }
  A <- (D < (1 + 1e-8)) & (D > (1 - 1e-8))
  A <- c(NA, A)
  D <- c(NA, D)
  if(drop.r){
    A[r] <- NA
    D[r] <- NA 
  }

  if (!is.circular){
    if (n2 >= 3){
      one <- v2[1:(n2-2)]
      two <- v2[2:(n2-1)]
      three <- v2[3:(n2)]
      
      TP <- (one > two & two < three) | (one < two & two > three)
      if (drop.r){
        for (i in which(r)){
          TP <- append(TP, NA, i - 2)
        }
      } else{
        for (i in which(r)){
          TP <- append(TP, F, i - 2)
        }
      }
      
      TP <- c(NA, NA, TP)  
    } else{
      TP <- NA
    }
    
  } else{
    TP <- NA
  }

  return(list(
    "R"=r,
    "A"=A,
    "TPF"=TP,
    "D"=D
  ))
}

.manage_alternatives <- function(a, s){
  # if (!is.null(a) & length(a) == 1){
  #   stop("a should be a vector of alternatives, not the number of alternatives")
  # }
  if (is.null(a)){
    a <- min(s, na.rm=T):max(s, na.rm=T)
  }
  return(a)
}

redundancy <- function(s, a=NULL){
  s <- s[!is.na(s)]
  a <- .manage_alternatives(a, s)
  
  n <- length(s)
  H_max <- log2(length(a))
  counts <- sapply(a, \(x){sum(s == x)})
  H_single <- log2(n) -
    sum(counts * log2(counts), na.rm = T) / n
  return(1 - H_single/H_max)
}

.response_matrix <- function(s, a, wrap=T){
  M <- matrix(0, length(a), length(a))
  for (i in 2:length(s)){
    last_index <- which(a == s[i-1])
    this_index <- which(a == s[i])
    M[last_index, this_index] <- M[last_index, this_index] + 1
  }
  if (wrap){
    last_index <- which(a == s[i])
    this_index <- which(a == s[1])
    M[last_index, this_index] <- M[last_index, this_index] + 1
  }
  return(M)
}

RNG <- function(s, a=NULL){
  s <- s[!is.na(s)]
  # Evans' RNG. Note T+Neil has a typo here
  a <- .manage_alternatives(a,s)
  M <- .response_matrix(s,a, wrap = T)
  
  row_sums <- rowSums(M)
  top <- sum(M * log10(M), na.rm = T)
  if (top == 0) return(0)
  bottom <- sum(row_sums * log10(row_sums), na.rm = T)
  top / bottom
}


all_measures <- function(s, a=NULL, drop.r=T, is.circular=FALSE){
  data.frame(t(unlist(
    lapply(
      all_measures_full(s, a, drop.r, is.circular), 
      FUN = \(x){apply(x, 2, mean, na.rm=T)}))))
}

all_measures_full <- function(s, a=NULL, drop.r=T, is.circular=FALSE){
  a <- .manage_alternatives(a,s)
  summ <- get_Summaries_full(s, drop.r = drop.r, is.circular=is.circular)
  
  list(  
    data.frame(
      R = summ[[1]],
      A = summ[[2]],
      TPF = summ[[3]],
      D = summ[[4]]
    ),
    data.frame(
      RNG = RNG(s, a),
      SWLZ_ER = randMeasures::SWLZ_ER(s),
      gzip=compute_gzip(s),
      lz76=randMeasures::lz76_complexity(s)
    )
  )
}
