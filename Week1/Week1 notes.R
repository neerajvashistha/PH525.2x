y<-rnorm(1e6)
x<-cbind(rep(1,1e6),rep(0,1,each=5e5))
beta<-c(1,1)
system.time({sum((y-x %*% beta)^2)})

system.time({sum(sapply(seq_along(y),
                        function(i) y[i] - x[i,1] * beta[1] -x[i,2] *beta[2])^2)})



#-----------------------------------
 
g<-9.8
n<-25
tt<-seq(0,3.4,len=n)

f<-56.67+0*tt-0.5*g*tt^2
y<-f+rnorm(n,sd=1)#error term becoz of obsrvational error.
plot(tt,y,xlab="Time in secs",ylab="Distance in meters")
lines(tt,f,col=2)


# Residual sum of square rss

rss<-function(Beta0,Beta1,Beta2){
  r<-y-(Beta0+Beta1*tt+Beta2*tt^2)
  sum(r^2)
}

Beta2s<-seq(-10,0,len=100)

RSS<-sapply(Beta2s,rss,Beta0=55,Beta1=0)
plot(Beta2s,RSS,type="l")
RSS<-sapply(Beta2s,rss,Beta0=65,Beta1=1)
lines(Beta2s,RSS,type = "l",col=3)


tt2<-tt^2
fit<-lm(y~tt+tt2)
summary(fit)


#-------
#matrix algerba to obt Least square estimates

X<-cbind(rep(1,length(tt)),tt,tt^2)
head(X)

Beta <- matrix(c(55,0,5),3,1)

r <- y - X %*% Beta
RSS <- t(r) %*% r

#or

RSS<-crossprod(r)

#Estimates
betahat <- solve(t(X) %*% X) %*% t(X) %*% y
#or
betahat <-  solve(crossprod(X))%*% crossprod(X,y)


QR <- qr(X)
Q <- qr.Q(QR)
R <- qr.R(QR)

backsolve(R,crossprod(Q,y))
