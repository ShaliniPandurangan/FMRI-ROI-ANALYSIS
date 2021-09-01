
a = (colSums(summary(diag.pc) != 0))
degree = as.numeric(attr(table(a), "dimnames")$a)
frequency = as.numeric(table(a))


maxval=1000
z=cbind(seq(0, maxval), 0)
z[degree + 1, 2] = frequency




nllik.pa=function(z, alpha){
  
  p = (1:maxval)^(-alpha)
  p=p/sum(p)
  
  -sum(z[-1, 2]*log(p)) 
}


nl=nlminb(start=1, objective=function(x) nllik.pa(z, x), lower=1e-3, upper=100)



p = (1:maxval)^(-nl$par)
p=p/sum(p)

#fitted powerlaw
plot(p, xlim=c(1, 30), ty="l")
points(z[-1, 1], z[-1, 2]/sum(z[-1, 2]), xlim=c(1, 30))

#cumulative powerlaw
plot(cumsum(p), xlim=c(1, 30), ylim=c(0,1), ty="l")
points(z[-1, 1], cumsum(z[-1, 2]/sum(z[-1, 2])), xlim=c(1, 30))



#lets put a cutoff at degree 5
nllik.pa=function(z, alpha, cutoff){
  
  p = (cutoff:maxval)^(-alpha)
  p=p/sum(p)
  
  -sum(z[(cutoff+1):nrow(z), 2]*log(p)) 
}

cutoff=5
nl=nlminb(start=1, objective=function(x) nllik.pa(z, x, cutoff), lower=1e-3, upper=100)



p = (cutoff:maxval)^(-nl$par)
p=p/sum(p)

#fitted powerlaw
plot(cutoff:maxval, p, xlim=c(1, 30), ty="l")
points(cutoff:maxval, z[(cutoff+1):nrow(z), 2]/sum(z[(cutoff+1):nrow(z), 2]), xlim=c(1, 30))

#cumulative powerlaw
plot(cutoff:maxval, cumsum(p), xlim=c(1, 30), ty="l")
points(cutoff:maxval, cumsum(z[(cutoff+1):nrow(z), 2]/sum(z[(cutoff+1):nrow(z), 2])), xlim=c(1, 30))






