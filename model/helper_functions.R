# Additional helper functions for BTR project

######
## Some manual functions
# Outer product of multiple vectors
mouter <- function(x1, ...) { 
  r <- x1
  for(vi in list(...)) r <- outer(r, vi)
  return(r)
}

# To reconstruct a tensor
reconstruct_tensor3 <- function(vec_list){
  rank = length(vec_list) # rank
  K = length(vec_list[[1]][,1])
  tensor = array(0, dim = c(K,K,K))
  for(r in 1:rank){
    tensor = tensor + mouter(vec_list[[r]][,1], vec_list[[r]][,2], vec_list[[r]][,3])
  }
  return(tensor)
}
#####