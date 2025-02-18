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
NumericVector slice(NumericVector v, int start, int end){
  int size = end - start + 1;
  NumericVector r(size);
  for (int i = 0; i < size; i++){
    r(i) = v(start + i);
  }
  return(r);
}
bool check_equal(NumericVector a, NumericVector b){
  bool equal = true;
  for (int i = 0; i < a.size(); i++){
    if (a(i) != b(i)){
      equal = false;
    }
  }
  return(equal);
}
bool find_match(NumericVector pattern, NumericVector target){
  int p_size = pattern.size();
  int k = 0;
  int ub = k + p_size-1;
  while(ub < target.size()){
    bool check = check_equal(pattern, slice(target, k, k+p_size-1));
    if (check){return(true);}
    k += 1; 
    ub = k + p_size - 1;
  }
  return(false);
}
// [[Rcpp::export]]
NumericVector SWLZ_lengths(NumericVector x) {
  int L = x.size();
  NumericVector lengths(L-1);
  for (int i = 1; i < L; i++){
    NumericVector previous_symbols = slice(x, 0, i-1);
    for (int l=i; l < L; l++){
      NumericVector test = slice(x, i, l);
      bool match = find_match(test, previous_symbols);
      if (!match){
        lengths(i-1) = l-i+1; 
        break;
      }
      if (l == (L-1)){
        // add longest + 1
        lengths(i-1) = l-i+2; 
        
      }
    }
  }

  return(lengths);
}
// [[Rcpp::export]]
double SWLZ_ER(NumericVector sequence){
  NumericVector swlz = SWLZ_lengths(sequence);
  return((std::log2((double)(sequence.size())))/ mean(swlz));
  
  // log2(length(sequence)) / mean(swlz)
}
