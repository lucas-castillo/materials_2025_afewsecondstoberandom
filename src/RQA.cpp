#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//


// [[Rcpp::export]]
NumericVector getDiagonals(NumericMatrix M){
  // Gets diagonals from a square matrix, 
  // returns a vector with the counts of each length
  int N = M.nrow();
  NumericVector counts(N); // in a square Matrix, the longest possible diag is as long as nrows
  for (int col = 0; col < N; col++){ // as many diags as 2N - 1 (N diags * 2 - 1 [not to count the main diag twice])
    int length = 0; // how long is the run of 1s
    int maxRows = N - col; // if we start in col N - 3, can only go three rows down until bump into matrix wall
    int addValue = col == 0 ? 1 : 2; // count diags once in main diag, twice otherwise (as matrix symmetric)
    for (int row = 0; row < maxRows; row++){
      int value = M(row, col+row); // we start in col x, row 0 and go down to row 1, col x + 1, etc.
      if (value == 1){ // if 1, increase length
        length++;
      } else{ // else we update histogram + length as we've finished the diagonal we were counting and starting a new one
        if (length > 0){
          counts(length - 1) +=  addValue; // length - 1 as cpp is 0-indexed  
        }
        length = 0;
      }
      if ( (row == (maxRows - 1)) & (length > 0) ){ // if we've ended the diagonal on a 1, we can't keep going! must update histogram
        counts(length - 1) += addValue; // no need to update length as end of loop
      }
    }
  }
  
  
  return counts;
}

// [[Rcpp::export]]
NumericVector getVerticals(NumericMatrix M){
  // Gets vertical lines from a square matrix, 
  // returns a vector with the counts of each length
  int N = M.nrow();
  NumericVector counts(N); // the longest possible vline is as long as nrows
  for (int col = 0; col < N; col++){
    int length = 0; // how long is the run of 1s
    for (int row = 0; row < N; row++){
      int value = M(row, col);
      if (value == 1){ // if 1, increase length
        length++;
      } else{ // else we update histogram + length as we've finished the line we were counting and starting a new one
        if (length > 0){
          counts(length - 1) += 1; 
        }
        length = 0;
      }
      if ( (row == (N - 1)) & (length > 0) ){ // if we've ended the line on a 1, we can't keep going! must update histogram
        counts(length - 1) += 1;  // no need to update length as end of loop
      }
    }
  }
  
  
  return counts;
}

//[[Rcpp::export]]
NumericMatrix ARMatrix(NumericVector sequence){
  int N = sequence.size();
  NumericMatrix M(N,N);
  
  // Fill matrix
  for (int i = 0; i < N; i++){
    for (int j = 0; j < N; j++){
      if (sequence(i) == sequence(j)){
        M(i,j) = 1;
      }
      
    }
  }
  return(M);
}

double get_entropy(NumericVector diagonalCounts, int lmin){
  double entropy = 0;
  double denominator = 0;
  for (int i = (lmin-1); i < diagonalCounts.size(); i++){
    denominator += diagonalCounts(i);
  }
  for (int i = (lmin-1); i < diagonalCounts.size(); i++){
    double probability = diagonalCounts(i) / denominator;
    if (probability > 0){
      entropy -= probability * log(probability);   
    }
    
  }
  return entropy;
}

// [[Rcpp::export]]
List RQA(NumericVector sequence, int lmin = 2, int vmin = 2) {
  // See Webber Jr. and Zbilut, 1994 for RR and DET; Marwan et al., 2002 for LAM
  // also http://www.recurrence-plot.tk/rqa.php in general

  NumericMatrix M = ARMatrix(sequence);
  int N = M.nrow();
  
  // RR is sum of 1s / N**2 
  int RR_sum= sum(M);
  double RR = (1/ pow(N, 2)) * RR_sum;
  
  // We need counts of each diag length for DET
  NumericVector diagHist = getDiagonals(M);
  double numerator = 0;
  double denominator = 0;
  for (int i = 0; i < diagHist.size(); i++){
    if (i >= (lmin - 1)){ 
      numerator += (i+1) * diagHist(i);
    }
    denominator += (i+1) * diagHist(i);
  }
  
  double DET = numerator / denominator;
  
  // We need counts of each vert length for LAM
  NumericVector vertHist = getVerticals(M);
  numerator = 0;
  denominator = 0;
  for (int i = 0; i < vertHist.size(); i++){
    if (i >= (vmin - 1)){ 
      numerator += (i+1) * vertHist(i);
    }
    denominator += (i+1) * vertHist(i);
  }
  
  double LAM = numerator / denominator;
  
  // longest diagonal [except main]
  int maxDiag = -1;
  for (int i = 0; i < (diagHist.size() - 1); i++){
    if (diagHist(i) > 0){
      maxDiag = i + 1;
    }
  }
  // longest vertical
  int maxVert = -1;
  for (int i = 0; i < vertHist.size(); i++){
    if (vertHist(i) > 0){
      maxVert = i + 1;
    }
  }
  // Trapping Time
  numerator = 0;
  denominator = 0;
  for (int i = 0; i < vertHist.size(); i++){
    if (i >= (vmin - 1)){ 
      numerator += (i+1) * vertHist(i);
      denominator += vertHist(i);
    }
  }
  double trapping_time = numerator / denominator;
  
  return List::create(
    Named("RR") = RR, 
    Named("DET") = DET, 
    Named("LAM") = LAM, 
    Named("maxL") = maxDiag, 
    Named("maxV") = maxVert,
    Named("TT") = trapping_time,
    Named("ENTR") = get_entropy(diagHist, lmin));
}

