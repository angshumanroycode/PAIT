T_phi<-function(x,y,B=200){
  if(is.numeric(x)==FALSE){
    stop("'x' is not a numeric vector!")
  }
  if(is.numeric(y)==FALSE){
    stop("'y' is not a numeric vector!")
  }
  n=length(x)
  if(n!=length(y)){
    stop("'x' and 'y' are not of equal length!")
  }
  return(T_phiC(x,y,B))
}

T_cor<-function(x,y,B=200){
  if(is.numeric(x)==FALSE){
    stop("'x' is not a numeric vector!")
  }
  if(is.numeric(y)==FALSE){
    stop("'y' is not a numeric vector!")
  }
  n=length(x)
  if(n!=length(y)){
    stop("'x' and 'y' are not of equal length!")
  }
  return(T_corC(x,y,B))
}

T_dcor<-function(x,y,B=200){
  if(is.numeric(x)==FALSE){
    stop("'x' is not a numeric vector!")
  }
  if(is.numeric(y)==FALSE){
    stop("'y' is not a numeric vector!")
  }
  n=length(x)
  if(n!=length(y)){
    stop("'x' and 'y' are not of equal length!")
  }
  return(T_dcorC(x,y,B))
}
