norm = function(x,p=2) {
  return(dist(rbind(x,double(length(x))),method="minkowski",p=p))
}
