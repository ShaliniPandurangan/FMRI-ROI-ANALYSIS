adj<-as.vector(AIC_NETWORK_DIAG_NO_WEIGHTS[1:70,]!=0)

x1<-c(rep(rep(c(0,1),each=35),35),rep(rep(c(1,0),each=35),35))
x2<-as.vector(rbind(cbind(matrix(0,ncol=35,nrow=35),diag(35)),cbind(diag(35),matrix(0,ncol=35,nrow=35))))
x3<- 

summary(glm(adj~x1+x2,family = binomial))
d[,i] <- as.vector(Y[i,,,1]) 

for (i in 1:70) {
  
  d[,i] <- as.vector(Y[i,,,1]) }

s<-s[complete.cases(as.longitudinal(d))]

s# <- as.longitudinal(d[405:9696,])

complete_obs <- complete.cases(s)

z <- s[complete.cases(s), ]
s <- d[complete.cases(d),]
