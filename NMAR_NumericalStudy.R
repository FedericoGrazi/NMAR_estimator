library(np)
library(boot)
options(np.messages = FALSE)
setwd("D:/Uni/KU - Courses/Advanced Non Parametric Statistics")

empirical_loss <- function(gg, p, h, Dl, Dm, l, delta, pn){
  psi1_m <- gen_funct(Dl$x,Dm$x,Dm$y2,Dm$Delta,1,h$bw, gamma = gg, funct = "psi")
  psi2_m <- gen_funct(Dl$x,Dm$x,Dm$y2,Dm$Delta,2,h$bw, gamma = gg, funct = "psi")
  eta1_m <- gen_funct(Dl$x,Dm$x,Dm$y2,Dm$Delta,1,h$bw, gamma = gg, funct = "eta")
  eta2_m <- gen_funct(Dl$x,Dm$x,Dm$y2,Dm$Delta,2,h$bw, gamma = gg, funct = "eta")
  m_moj_m <- eta1_m + (psi1_m / psi2_m ) * (1 - eta2_m)
  
  
  EmpLp <- (1/nrow(Dl))*sum( Dl$Delta*abs(m_moj_m-Dl$y2)^p + (1-Dl$Delta)*(delta/pn)*abs(m_moj_m-Dl$y2)^p )
  return(EmpLp)
}
moj_model <- function(data,perc = 0.7,M = 100, alpha = 0.01, lambda = 0.95, p = 2, subsample = F, verbose = T){
  
  n <- nrow(data)
  pn <- ((log(n)^.25)/(n*lambda)^(1-alpha))^.5
  
  
  m <- floor(perc * n)
  l <- n-m
  
  train <- sample(n,m, replace = F)
  
  Dm <- data[train,]
  Dl <- data[-train,]
  
  delta <- rbinom(l,1,pn )
  # sum(delta)
  # sum((1-Dl$Delta)*delta)
  ripescati <- (1-Dl$Delta)*delta
  # points(Dl$x[ripescati==1],Dl$y2[ripescati==1], col = "cyan", pch = 19)
  if(verbose) cat(paste("Resampled", sum(ripescati), "Non-Respondants \n"))
  
  h_m <- np::npregbw(Dm$y2~Dm$x)
  if(verbose) cat(paste("h_m Obtained ✔ \n"))
  
  optimL2 <- optimize(empirical_loss,interval = c(-M,M), p = p, h = h_m, Dl = Dl, Dm = Dm, delta = delta, pn = pn)
  gamma_hat_L2 <- optimL2$minimum
  if(verbose) cat(paste("Estimated Gamma:", round(gamma_hat_L2,3), "\n"))
  
  
  # Optimal L2 
  psi1_L2m <- gen_funct(Dm$x,Dm$x,Dm$y2,Dm$Delta,1, gamma = gamma_hat_L2, funct = "psi"); if(verbose) cat(paste("Psi1 Function Done ✔ \n"))
  psi2_L2m <- gen_funct(Dm$x,Dm$x,Dm$y2,Dm$Delta,2, gamma = gamma_hat_L2, funct = "psi"); if(verbose) cat(paste("Psi2 Function Done ✔ \n"))
  eta1_L2m <- gen_funct(Dm$x,Dm$x,Dm$y2,Dm$Delta,1, gamma = gamma_hat_L2, funct = "eta"); if(verbose) cat(paste("Eta1 Function Done ✔ \n"))
  eta2_L2m <- gen_funct(Dm$x,Dm$x,Dm$y2,Dm$Delta,2, gamma = gamma_hat_L2, funct = "eta"); if(verbose) cat(paste("Eta2 Function Done ✔ \n"))
  
  # Proposed Estimator
  m_moj_L2m <- eta1_L2m + (psi1_L2m / psi2_L2m ) * (1 - eta2_L2m)
  df_moj_L2m <- data.frame(x = Dm$x, m_moj_L2m)
  out <- df_moj_L2m[order(df_moj_L2m$x),]
  if(verbose) cat(paste("Regression Estimate Done ✔ \n"))
  
  if(subsample) out <- list(model = out, subsample = Dl[ripescati ==1,-3])
  return(out)
}
gen_funct <- function(x,Xi,Yi,D,k,h = NULL, gamma, funct = c("eta", "psi", "nw")){
  
  if(is.null(h)) h <- npregbw(Yi~Xi)$bw
  if(funct == "psi") phi_val <- phi(Yi,gamma) else phi_val <- rep(1,length(Yi))  
  
  
  # m <- sapply(x, \(xx) sum(dnorm((xx-Xi)/h,0,1) * Yi^(2-k) * D * phi_val) / sum(dnorm((xx-Xi)/h,0,1)) )
  m <- sapply(x,
              function(xx, Xi, Yi, h,k,D,phi_val){
                sum(dnorm((xx-Xi)/h,0,1) * Yi^(2-k) * D * phi_val)/sum(dnorm((xx-Xi)/h,0,1))
              },
              h = h,
              Xi = Xi,
              Yi = Yi,
              D = D,
              k = k,
              phi_val = phi_val)
  return(m)
}
phi <- function(Yi, gamma) exp(gamma * y)


N <- 500
sigma <- matrix(0,5,5)
for(i in 1:5){
  for(j in 1:5){
    sigma[i,j] <- 2^{-abs(i-j)}
  }
}
X <- MASS::mvrnorm(N,rep(1,5), sigma)
y <- 2.6 - X[,1] - X[,2]^2 + X[,3]*X[,4] + exp(-X[,5])+ rnorm(N, 0, .5)

plot(density(y))

par(mfrow = c(2,3))
plot(X[,1],y,pch = 19, cex = .5);  lines(sort(X[,1]),fitted(smooth.spline(X[,1],y))[order(X[,1])],col = "red", lwd = 2)
plot(X[,2],y,pch = 19, cex = .5);  lines(sort(X[,2]),fitted(smooth.spline(X[,2],y))[order(X[,2])],col = "red", lwd = 2)
plot(X[,3],y,pch = 19, cex = .5);  lines(sort(X[,3]),fitted(smooth.spline(X[,3],y))[order(X[,3])],col = "red", lwd = 2)
plot(X[,4],y,pch = 19, cex = .5);  lines(sort(X[,4]),fitted(smooth.spline(X[,4],y))[order(X[,4])],col = "red", lwd = 2)
plot(X[,5],y,pch = 19, cex = .5);  lines(sort(X[,5]),fitted(smooth.spline(X[,5],y))[order(X[,5])],col = "red", lwd = 2)


gamma <- -.98
phi_val <- phi(y,gamma)
beta <- c(.6,.8,.25,-.35,-.3,.75)

gx <- cbind(1, X) %*% beta
prob <- (1 + exp(gx) * phi_val)^-1
Delta <-  1-rbinom(N,1,prob)
mean(Delta)


fitted_nw <- fitted_moj <- list(mX1 = NULL,mX2 = NULL,mX3 = NULL,mX4 = NULL,mX5 = NULL)
for(i in 1:5){ 
  fitted_nw[[i]] <- npreg(y[Delta == 1] ~ X[Delta == 1,i])
  dfi <- data.frame(x =fitted_nw[[i]]$eval[[1]], y2 = fitted(fitted_nw[[i]]))
  fitted_nw[[i]] <- dfi[order(dfi$x),]
  fitted_moj[[i]] <- moj_model(data.frame(x = X[,i], y2 = y, Delta = Delta),subsample  = T)
}



par(mfrow = c(2,3))
plot(X[,1],y, pch = ((1-Delta)*3)+1, col = (2-Delta), cex = .4)
lines(fitted_nw[[1]], col = "blue")
lines(fitted_moj[[1]]$model, col = "forestgreen")
plot(X[,2],y, pch = ((1-Delta)*3)+1, col = (2-Delta), cex = .4)
lines(fitted_nw[[2]], col = "blue")
lines(fitted_moj[[2]]$model, col = "forestgreen")
plot(X[,3],y, pch = ((1-Delta)*3)+1, col = (2-Delta), cex = .4)
lines(fitted_nw[[3]], col = "blue")
lines(fitted_moj[[3]]$model, col = "forestgreen")
plot(X[,4],y, pch = ((1-Delta)*3)+1, col = (2-Delta), cex = .4)
lines(fitted_nw[[4]], col = "blue")
lines(fitted_moj[[4]]$model, col = "forestgreen")
plot(X[,5],y, pch = ((1-Delta)*3)+1, col = (2-Delta), cex = .4)
lines(fitted_nw[[5]], col = "blue")
lines(fitted_moj[[5]]$model, col = "forestgreen")

