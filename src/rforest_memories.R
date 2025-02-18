library(ranger)
library(abcrf)
library(foreach)
library(doParallel)

# Def functions -----------------------------------------------------------
find_leaf <- function(tree, datum){
  currentNodeID <- 0
  terminal <- F
  while (!terminal){
    currentNode <- tree[tree$nodeID == currentNodeID, ]
    terminal <- currentNode$terminal
    if (!terminal){
      if (datum[[currentNode$splitvarName]] > currentNode$splitval){
        currentNodeID <- currentNode$rightChild
      } else{
        currentNodeID <- currentNode$leftChild
      }
    }
  }
  return(c(currentNodeID, as.numeric(currentNode$prediction)))
}
get_tree_memory <- function(tree, training){
  M <- matrix(0, nrow=nrow(tree), ncol=length(levels(tree$prediction)))
  
  for (row in 1:nrow(training)){
    datum <- training[row,]
    RES <- find_leaf(tree, datum)
    true_label <- as.numeric(datum$f_interval)
    M[RES[1], true_label] <- M[RES[1], true_label] + 1
  }
  M
}

get_forest_memory <- function(abcrf_model, training, n_trees=NULL){
  forest <- abcrf_model$model.rf
  # add lda -----
  training.p <- cbind(training, predict(abcrf_model$model.lda, training)$x)
  if (is.null(n_trees)){
    n_trees <- forest$num.trees
  }
  n_labels <- length(levels(forest$predictions))  
  
  M <- foreach (
    treeIndex=1:n_trees, 
    .packages = c("abcrf", "ranger"),
    .export = c("get_tree_memory", "find_leaf")
    ) %dopar% {
    tree <- treeInfo(forest, treeIndex)
    get_tree_memory(tree, training.p)
  }
  M
}
