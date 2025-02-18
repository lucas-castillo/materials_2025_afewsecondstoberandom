# Adapted from Angelike, T., & Musch, J. (2024). 
# A comparative evaluation of measures to assess randomness in human-generated 
# sequences. Behavior Research Methods, 56(7), 7831–7848. https://doi.org/10.3758/s13428-024-02456-7
# Code in OSF Repository: https://osf.io/bycth 

#' wrapper function that computes all available randomness indices/measures
#'
#' @param sequences sequences over which the randomness measures will be computed
#' @param min_val minimum value a number can take in the sequences
#' @param max_val maximum value a number can take in the sequences
#' @param options number of distinct values a sequence can consist of (default equals
#' the maximum value if the set of numbers that can be used for generating the
#' sequences starts at 1)
#' @param min_block_size lower bound for the block size used for computation of
#' entropy, complexity, and phi index measures
#' @param max_block_size upper bound for the block size used for computation of
#' entropy, complexity, and phi index measures
#' @param max_block_size_phi upper bound for the block size used for computation of
#' the phi index. This is by default equal to \code{max_block_size} but may be
#' set to a lower value due to its computationally highly time-intensive nature.
compute_all_measures <- function(sequences,
                                 min_val,
                                 max_val,
                                 options = max_val,
                                 min_block_size = 2,
                                 max_block_size = 10,
                                 max_block_size_phi = max_block_size) {

  if (min_block_size == 1 | max_block_size > 12) {
    stop("This wrapper function can only be used for the block sizes 2 to 12. You may want to use other functions in this script instead.")
  }

  ## compute block entropy over different block sizes
  entropy_values <- block_entropy_df(sequences, min_block_size, max_block_size)
  all_indices <- entropy_values

  ## compute traditional measures of randomness
  traditional_indices <- compute_randseqR_indices(sequences, min_val = min_val, max_val = max_val, max_block_size_phi)
  all_indices <- cbind(all_indices, traditional_indices)

  ## compute compression algorithms
  lz_complexity <- compute_compression_algorithms(sequences)
  all_indices <- cbind(all_indices, lz_complexity)

  ## compute measures of algorithmic complexity / BDM
  complexity_indices <- compute_algorithmic_complexity(sequences, options = options, min_block_size, max_block_size)
  all_indices <- cbind(all_indices, complexity_indices)

  return(all_indices)
}

#' function for computing block entropy for different block sizes
#'
#' @param data sequences over which the randomness measures will be computed
#' @param min_block_size lower bound for the block size used for computation of
#' entropy
#' @param max_block_size upper bound for the block size used for computation of
#' entropy
block_entropy_df <- function(data, min_block_size = 1, max_block_size = 1) {
  ## prepare results data frame
  results <- matrix(0, nrow = nrow(data), ncol = max_block_size - min_block_size + 1)
  results <- as.data.frame(results)
  colnames(results) <- paste0("block_entropy_", min_block_size:max_block_size)

  ## compute block entropy for all block sizes
  for (i in min_block_size:max_block_size) {
    current_col_name <- paste0("block_entropy_", i)
    results[, current_col_name] <- apply(data, MARGIN = 1, function(x) randfindR::block_entropy(x, i))
  }

  return(results)
}

#' wrapper function that computes all randseqR measures of randomness
#' for each sequence and extracts them from a list into a data frame
#'
#' @param sequences sequences over which the randomness measures will be computed
#' @param min_val minimum value a number can take in the sequences
#' @param max_val maximum value a number can take in the sequences
#' @param min_block_size lower bound for the block size used for computation of
#' phi index measures
#' @param max_block_size upper bound for the block size used for computation of
#' phi index measures
#'
#' @references
#' Oomens, W., Maes, J. H. R., Hasselman, F., & Egger, J. I. M. (2021).
#' RandseqR: An R Package for Describing Performance on the Random Number
#' Generation Task. Frontiers in Psychology, 12.
#' https://doi.org/10.3389/fpsyg.2021.629012
compute_randseqR_indices <- function(data, min_val, max_val, max_block_size) {

  empty_col <- numeric(length = nrow(data))
  phi_names <- paste0("phi_", 2:max_block_size)

  col_names_without_phi <-
    c(
      "rng",
      "rng_2",
      "coupon",
      "repetition_mean",
      "repetition_median",
      "repetition_mode",
      "null_score",
      "adjacency_asc",
      "adjacency_desc",
      "adjacency_combi",
      "tp_index",
      "runs_index",
      "redundancy_index"
    )

  all_col_names <- c(col_names_without_phi, phi_names)
  randseqR_indices <- matrix(0, nrow = nrow(data), ncol = length(all_col_names))
  randseqR_indices <- as.data.frame(randseqR_indices)
  colnames(randseqR_indices) <- all_col_names

  ## iterate over all rows and compute randomness indices with randseqR
  for (p in 1:nrow(data)) {
    sequence <- as.numeric(data[p, ])
    tryCatch({
      ## compute variants of gap score
      rep_gap_results <- randseqR::repGap(sequence, minScale = min_val, maxScale = max_val)
      randseqR_indices$repetition_mean[p] <- as.numeric(rep_gap_results$RG_mean)
      randseqR_indices$repetition_median[p] <- as.numeric(rep_gap_results$RG_median)
      randseqR_indices$repetition_mode[p] <- as.numeric(rep_gap_results$RG_mode)

      ## compute variants of adjacency index
      adjacency_results <- randseqR::Adjacency(sequence, minScale = min_val, maxScale = max_val)
      randseqR_indices$adjacency_asc[p] <- as.numeric(adjacency_results$Ascending)
      randseqR_indices$adjacency_desc[p] <- as.numeric(adjacency_results$Descending)
      randseqR_indices$adjacency_combi[p] <- as.numeric(adjacency_results$Combined)

      ## compute variants of phi index
      phi_index_results <- randseqR::phiIndex(sequence, minScale = min_val, maxScale = max_val, maxOrder = max_block_size)
      for (i in 1:length(phi_names)) {
        phi_col_name <- phi_names[i]
        randseqR_indices[p, phi_col_name] <- as.numeric(phi_index_results[[i]])
      }

      ## compute remaining measures
      randseqR_indices$rng[p] <- as.numeric(RNG(sequence, minScale = min_val, maxScale = max_val)[[1]])
      randseqR_indices$rng_2[p] <- as.numeric(RNG2(sequence, minScale = min_val, maxScale = max_val)[[1]])
      randseqR_indices$coupon[p] <- as.numeric(Coupon(sequence, minScale = min_val, maxScale = max_val)[[1]])
      randseqR_indices$null_score[p] <- as.numeric(NSQ(sequence, minScale = min_val, maxScale = max_val)[[1]])
      randseqR_indices$tp_index[p] <- as.numeric(TPI(sequence, minScale = min_val, maxScale = max_val)[[1]])
      randseqR_indices$runs_index[p] <- as.numeric(Runs(sequence, minScale = min_val, maxScale = max_val)[[1]])
      randseqR_indices$redundancy_index[p] <- as.numeric(Redundancy(sequence, minScale = min_val, maxScale = max_val)[[1]])
    },
    error = function(e) {
      ## if the function causes unexpected errors, store NA
      randseqR_indices[p, all_col_names] <- NA
    })

  }
  return(randseqR_indices)
}

compute_gzip <- function(v){
  v <- v[!is.na(v)]
  v <- v - min(v) + 1
  codes <- c((10*16+2):(16*25))
  ascii <- sapply(codes, intToUtf8)
  length(memCompress(paste0(ascii[v], collapse = ""), type = "gzip"))
}