table(a)



maxval=100

z=cbind(seq(0, maxval), 0)

z[as.numeric(names(table(a)))+1,2] <-table(a) 

#z[1:NROW(z), 2]=z[, 2]




nllik.pa=function(z, alpha){
  
  p = (1:(nrow(z)-1))^(-alpha)
 # p=cumsum(p)
  p=p/sum(p)
  
  -sum(z[-1, 2]*log(p)) 
  
}

nl=nlminb(start=1, objective=function(x) nllik.pa(z, x), lower=1e-3, upper=4)

p<-(1:maxval)^(-nl$par)

p<-p/sum(p)

hist(a[a!=0],probability = TRUE,breaks = 20) 

lines(1:max(a), p[1:max(a)])

#nl=#nlminb(start=1, objective=function(x) nllik.pa(z, x), lower=1e-3, upper=4)







