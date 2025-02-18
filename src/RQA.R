Rcpp::sourceCpp("src/RQA.cpp")

plot_RP <- function(seq=NULL, name="Random", crqaMatrix = NULL){
  if (is.null(crqaMatrix)){
    N <- length(seq)
    M <- ARMatrix(seq)  
  } else{
    M <- matrix(as.numeric(crqaMatrix), nrow = nrow(crqaMatrix))  
    N <- nrow(M)
  }
  
  
  G <- (
    ggplot(data.frame(x = rep(1:N, N), y = rep(1:N, each=N), z = as.vector(M)), aes(x = x, y = y, fill = z)) 
    + geom_raster() 
    + scale_fill_gradientn(colours=c("white","black")) 
    + theme_minimal() 
    + theme(
      legend.position = "none", panel.grid = element_blank(), 
      axis.title = element_blank(), 
      axis.text = element_blank(), 
      plot.title = element_text(face = "bold", color = "red") 
    )
  )
  
  G <- G + labs(title=name)
  return(G)
}
