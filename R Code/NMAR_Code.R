library(np)
library(boot)
library(foreach)
options(np.messages = FALSE)
setwd("D:/Uni/KU - Courses/Advanced Non Parametric Statistics/Project")

# f <- function(x) -5+4*atan(2 * x - 10) - 7*atan(2 * x - 5)+  x^4 / 1e3
# gx <- function(x) (.4 + .1 *x - x^3/100)/10
# phi <- function(y,gamma) exp(gamma * (y^3/500))
# pixy <- function(x,y,gamma) 1/(1+1.8*exp(gx(x)) * phi(y,gamma))
f <- function(x) 2.6+sin(.1 * x + 2)^2 *x * cos(.1 * x +3)^2 + (.1 * x)^3
gx <- function(x) .5+ .25 * x
phi <- function(y,gamma) exp(gamma * y)
pixy <- function(x,y,gamma) 1/(1+exp(gx(x)) * phi(y,gamma))
gen_data <- function(n, noise, gamma ){
  XX <- runif(n,-10,10)
  YY <-  f(XX) + rnorm(n,0,noise)
  YY <- as.numeric(YY)
  pi_vali <- pixy(XX,YY,gamma)
  Deltai <-  1-rbinom(n,1,pi_vali)
  out <- data.frame(x = sort(XX), y2 = YY[order(XX)], Delta = Deltai[order(XX)])
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
empirical_loss <- function(gg, p, h = NULL, Dl=NULL, Dm= NULL, delta = NULL, pn= NULL, param = NULL){
  if(!is.null(param)){
    Dm <- param$Dm
    Dl <- param$Dl
    alpha = param$alpha
    delta = param$delta
    pn = param$pn
  }
  if(is.null(h))  h <-  npregbw(Dm$y2~Dm$x)
  psi1_l <- gen_funct(Dl$x,Dm$x,Dm$y2,Dm$Delta,1,h = h$bw, gamma = gg, funct = "psi")
  psi2_l <- gen_funct(Dl$x,Dm$x,Dm$y2,Dm$Delta,2,h = h$bw, gamma = gg, funct = "psi")
  eta1_l <- gen_funct(Dl$x,Dm$x,Dm$y2,Dm$Delta,1,h = h$bw, gamma = gg, funct = "eta")
  eta2_l <- gen_funct(Dl$x,Dm$x,Dm$y2,Dm$Delta,2,h = h$bw, gamma = gg, funct = "eta")
  l_moj <- eta1_l + (psi1_l / psi2_l ) * (1 - eta2_l)
  
  
  EmpLp <- (1/nrow(Dl))*sum( Dl$Delta*abs(l_moj-Dl$y2)^p + (1-Dl$Delta)*(delta/pn)*abs(l_moj-Dl$y2)^p )
  return(EmpLp)
}
moj_model <- function(data,perc = 0.7,M = 100, alpha = 0.01, lambda = 0.95, p = 2, verbose = T, x= NULL, extras = F,seed =NULL){

  n <- nrow(data)
  pn <- ((log(n)^.25)/(n*lambda)^(1-alpha))^.5
  
  m <- floor(perc * n)
  l <- n-m
  
  train <- sample(n,m, replace = F)
  
  Dm <- data[train,]
  Dl <- data[-train,]

  if(!is.null(seed)) set.seed(seed)
  ripescati <- 0
  while(sum(ripescati)==0){
    ripescati <- 0
    delta <- rbinom(l,1,pn)
    ripescati <- (1-Dl$Delta)*delta
  }
  # points(Dl$x[ripescati==1],Dl$y2[ripescati==1], col = "cyan", pch = 19)
  if(verbose) cat(paste("Resampled", sum(ripescati), "Non-Respondants \n"))

  h_m <- np::npregbw(Dm$y2~Dm$x)
  if(verbose) cat(paste("h_m Obtained ✔ \n"))

  optimL2 <- optimize(empirical_loss,interval = c(-M,M), p = p, h = h_m, Dl = Dl, Dm = Dm, delta = delta, pn = pn,maximum = F)
  gamma_hat_L2 <- optimL2$minimum
  if(verbose) cat(paste("Estimated Gamma:", round(gamma_hat_L2,3), "\n"))
  
  
  # Optimal L2 
  if(is.null(x)){
    psi1_L2m <- gen_funct(Dm$x,Dm$x,Dm$y2,Dm$Delta,1, gamma = gamma_hat_L2, funct = "psi"); if(verbose) cat(paste("Psi1 Function Done ✔ \n"))
    psi2_L2m <- gen_funct(Dm$x,Dm$x,Dm$y2,Dm$Delta,2, gamma = gamma_hat_L2, funct = "psi"); if(verbose) cat(paste("Psi2 Function Done ✔ \n"))
    eta1_L2m <- gen_funct(Dm$x,Dm$x,Dm$y2,Dm$Delta,1, gamma = gamma_hat_L2, funct = "eta"); if(verbose) cat(paste("Eta1 Function Done ✔ \n"))
    eta2_L2m <- gen_funct(Dm$x,Dm$x,Dm$y2,Dm$Delta,2, gamma = gamma_hat_L2, funct = "eta"); if(verbose) cat(paste("Eta2 Function Done ✔ \n"))
    m_moj_L2m <- eta1_L2m + (psi1_L2m / psi2_L2m ) * (1 - eta2_L2m)
    df_moj_L2m <- data.frame(x = sort(Dm$x), m_moj_L2m = m_moj_L2m[order(Dm$x)])
  }else{
    psi1_L2m <- gen_funct(x,Dm$x,Dm$y2,Dm$Delta,1, gamma = gamma_hat_L2, funct = "psi"); if(verbose) cat(paste("Psi1 Function Done ✔ \n"))
    psi2_L2m <- gen_funct(x,Dm$x,Dm$y2,Dm$Delta,2, gamma = gamma_hat_L2, funct = "psi"); if(verbose) cat(paste("Psi2 Function Done ✔ \n"))
    eta1_L2m <- gen_funct(x,Dm$x,Dm$y2,Dm$Delta,1, gamma = gamma_hat_L2, funct = "eta"); if(verbose) cat(paste("Eta1 Function Done ✔ \n"))
    eta2_L2m <- gen_funct(x,Dm$x,Dm$y2,Dm$Delta,2, gamma = gamma_hat_L2, funct = "eta"); if(verbose) cat(paste("Eta2 Function Done ✔ \n"))
    m_moj_L2m <- eta1_L2m + (psi1_L2m / psi2_L2m ) * (1 - eta2_L2m)
    df_moj_L2m <- data.frame(x = sort(x), m_moj_L2m = m_moj_L2m[order(x)])
    
  }
  # Proposed Estimator
  out <- df_moj_L2m
  if(verbose) cat(paste("Regression Estimate Done ✔ \n"))

  if(extras) out <- list(model = out, subsample = Dl[ripescati ==1,-3], param = list(Dm = Dm, Dl = Dl, delta = delta, pn = pn, alpha = alpha, lambda = lambda, indices = train, h = h_m))
  return(out)
}
m0 <- function(eval.x,Xi, Yi, D,gamma, h = NULL){
  
  if(is.null(h)) h <- npregbw(Yi~Xi)$bw
  
  m <- sapply(eval.x,
              function(xx, Xi, Yi, h,k,D,phi_val){
                num_w1i <- (D * dnorm((xx-Xi)/h,0,1) * phi(y = Yi, gamma = gamma))
                den_w1i <- sum(D * dnorm((xx-Xi)/h,0,1) * phi(y = Yi, gamma = gamma))
                sum( (num_w1i/den_w1i) * Yi ) 
              },
              h = h,
              Xi = Xi,
              Yi = Yi,
              D = D)
  return(m)
  
}
kim_gamma <- function(gamma, eval.x, eval.y, Xi, Yi,D, h = NULL, minimize = F){
  y_opt <- m0(eval.x,Xi,Yi,D,gamma = gamma, h = h)
  obj <- sum((eval.y - y_opt))
  if(minimize) obj <- obj^2
  return(obj)
  
}
full_delta <- function(moj_out, n){
  delta_df <- data.frame(indx = c(moj_out$param$indices,(1:n)[-moj_out$param$indices]),
                         deltas = c(rep(0, length(moj_out$param$indices)), moj_out$param$delta))
  delta_vect <- delta_df[order(delta_df$indx),2]
  return(delta_vect)
}
kim_est <- function(Xi,Yi, D, d,h =NULL, opt = F){
  ssi <- (1-D) * d
  kimopt <- optimize(kim_gamma,
                     lower = -10, upper = 10, 
                     eval.x = Xi[ssi==1], 
                     eval.y = Yi[ssi==1],
                     h = np::npregbw(Yi[D == 1] ~ Xi[D == 1])$bw,
                     Xi = Xi,
                     Yi = Yi,
                     D = D,
                     minimize = T)
  
  kim <- m0(eval.x = Xi,Xi,Yi,D,gamma = kimopt$minimum)
  kim_vals <- (D * Yi + (1-D) * kim)
  kimnw <- gen_funct(Xi,Xi,kim_vals, D = 1, k = 1, gamma = 1, funct = "nw")
  out <- data.frame(x = Xi, kim = kimnw)
  if(opt) out <- list(model = out, gamma.opt = kimopt$minimum, h = np::npregbw(kim_vals~Xi)) 
  return(out)
  }


# Define the "unknown" quantities #####

# true regression function
n=1001
set.seed(42)

# Generate data 
gamma <- .98
data <- gen_data(n,1,gamma); x <- data$x; y2 <- data$y2; Delta <- data$Delta
mean(Delta)

data_obs <- subset(data, Delta == 1)
par(mar = c(3,4,3,4))
with(data, plot(x,y2, cex = .6,pch = ((1-Delta)*3)+1, col = (2-Delta), ylab = "y"))
legend("topleft",cex = 0.8, col = c("black","red"), pch = c(1,4), legend = c("Observed", "Missing"), bty = "n")

# Plot Estimators Functions ####
# curve(f, from = -10, to =10, add = T, col = "red", lwd = 2)

h_cc <- np::npregbw(y2~x)
h_obs <- np::npregbw(data_obs$y2~data_obs$x)
# h <- 1

# Nadaraya - Watson
m_nw <- np::npreg(h_obs)
fitted_nw <- data.frame(x = data_obs$x, fitted = fitted(m_nw))
fitted_nw <- fitted_nw[order(fitted_nw$x),]
lines(fitted_nw , col = "royalblue3", lwd = 2)

# NW Completed
# m_nw_cc <- np::npreg(h_cc)
# fitted_nw_cc <- data.frame(x = x, fitted = fitted(m_nw_cc))
# fitted_nw_cc <- fitted_nw_cc[order(fitted_nw_cc$x),]
# lines(fitted_nw_cc , col = "black", lwd = 2)


# Proposed Estimators ####

df_moj_L2m <- moj_model(data = data, extras = T, seed = 42)

lines(df_moj_L2m$model$x, df_moj_L2m$model$m_moj_L2m, col = "springgreen3", lwd = 2)
with(df_moj_L2m$subsample,points(x,y2, col = "darkturquoise", pch = 19, cex = 0.8))

# Add points to NW

data_add <- rbind(data_obs, cbind(df_moj_L2m$subsample,Delta = 1))
h_add <- np::npregbw(data_add$y2~data_add$x)

m_nw_res <- np::npreg(h_add)
fitted_nw_res <- data.frame(x = data_add$x, fitted = fitted(m_nw_res))
fitted_nw_res <- fitted_nw_res[order(fitted_nw_res$x),]
lines(fitted_nw_res , col = "darkturquoise", lwd = 2)


# Kim and Yu (2011) Estimator
delta_vect <- full_delta(df_moj_L2m, nrow(data))
fitted_kim <- kim_est(x, y2, Delta, delta_vect, opt = T)
lines(fitted_kim$model, col = "darkorchid3", lwd = 2)

with(data, plot(x,y2, cex = .4,pch = ((1-Delta)*3)+1, col = scales::alpha(2-Delta,.45), ylab = "y"))
with(df_moj_L2m$subsample,points(x,y2, col = "darkturquoise", pch = 19, cex = 0.8))
lines(fitted_nw_res , col = "darkturquoise", lwd = 2)
lines(fitted_kim$model, col = "darkorchid3", lwd = 2)
lines(df_moj_L2m$model$x, df_moj_L2m$model$m_moj_L2m, col = "springgreen3", lwd = 2)
lines(fitted_nw , col = "royalblue3", lwd = 2)
legend("bottomright",legend = c("Nadaraya-Watson", "Mojirsheibani Est.", 
                  "Nadaraya-Watson + follow up", "Kim and Yu", 
                  "Follow-Up Sub Sample"), col = c("royalblue3", "springgreen3", "darkturquoise", "darkorchid3", "darkturquoise"), 
       pch = c(NA, NA, NA, NA, 19),  
       lty = c(1, 1, 1, 1, NA),      
       lwd = c(2, 2, 2, 2, NA),      
       cex = 0.8, 
       bty = "n")

# Convergence Rate ####
# NW Complete
bw_cc <- np::npregbw(data$y2~data$x)
rate_cc <- (bw_cc$nobs * bw_cc$bw)^-1 + bw_cc$bw^2

# NW with missing data
bw_md <- np::npregbw(data_obs$y2~data_obs$x)
rate_md <- (bw_md$nobs * bw_md$bw)^-1 + bw_md$bw^2

# Proposed Estimator
bw_moj <- df_moj_L2m$param$h
pn <- ((log(bw_moj$nobs)^.25)/(bw_moj$nobs*df_moj_L2m$param$lambda)^(1-df_moj_L2m$param$alpha))^.5
rate_moj <- sqrt(log(n)/(n*bw_moj$bw*pn^2))
rate_moj_fixed <- sqrt(log(n)/(n*bw_moj$bw))

# NW with Kim values
bw_kim <- fitted_kim$h
rate_kim <- (bw_kim$nobs * bw_kim$bw)^-1 + bw_kim$bw^2

rates <- 
  data.frame(
  model = c("Complete NW", "NW with  Missing Data", "NW with Kim Values", "Moj with pn = o(1)", "Moj with pn = c"),
  bw = c(bw_cc$bw, bw_md$bw, bw_kim$bw, bw_moj$bw, bw_moj$bw),
  n = c(bw_cc$nobs, bw_md$nobs, bw_kim$nobs, bw_moj$nobs, bw_moj$nobs),
  rate = c(rate_cc, rate_md, rate_kim, rate_moj, rate_moj_fixed)
)
# Gamma Optimization ####

# Mojirsheibani and Khudaverdyan (2024)

ggg <- seq(-15,15,by = 0.1)
# emp2 <- foreach(gval = ggg, .combine = c) %do% empirical_loss(gval, p=2, param = df_moj_L2m$param)
# emp1 <- foreach(gval = ggg, .combine = c) %do% empirical_loss(gval, p=1, param = df_moj_L2m$param)
load("EmpLossS.rda"); emp1 <- EmpLossS$emp1; emp2 <- EmpLossS$emp2; kimyu <- EmpLossS$kimyu

par(mar = c(5, 4, 4, 4) + 0.3)  
plot(ggg, emp2, type = "l", ylab ="L2 Loss", xlab = "Gamma", main = "Empirical Loss")
par(new = TRUE)
plot(ggg, emp1, type = "l", axes = FALSE, bty = "n", xlab = "", ylab = "", col = "blue")
axis(side=4, at = pretty(range(emp1, na.rm = T)), col="blue",col.axis="blue")
mtext("L1 Loss", side=4, line=3, col ="blue")
abline(v = ggg[which.min(emp2)], lty = "dashed")
abline(v = ggg[which.min(emp1)], lty = "dashed", col = "blue")
abline(v = gamma, lty = "dashed", col = "red")
legend(x = "bottomright", col = c("black","blue", "red"),lty = 2, lwd = 2, bty = "n", legend = c("L2 minimizer", "L1 minimizer", "True Gamma"), cex = 0.8)



par(mar = c(5, 4, 4, 4))  
plot(sort(data$y2), ((1+(phi(y = data$y2, gamma = min(emp2,na.rm = T))))^-1)[order(data$y2)], type = "l", xlab = "y",lwd = 2, ylab = "Propensity", main = "Missing Probability")
# lines(sort(data$y2), ((1+(phi(y = data$y2, gamma = gamma)))^-1)[order(data$y2)], col = "firebrick4", lwd = 2)
points(data$y2, pixy(data$x, data$y2, gamma), pch = ".", col = "red", cex = 1.2)
legend(x = "topright", col = c("black", "red"),lty = c( 1,NA) ,pch = c(NA, 19), lwd = c(2,NA), bty = "n", legend = c("Estimated (γ = -2.7)", "True Probabilities"), cex = 0.8)

# Kim and Yu (2011)  

# kimyu <- foreach(gval = ggg, .combine = c) %do%  kim_gamma(gval,eval.x = x[(1-Delta)*delta_vect ==1],eval.y = y2[(1-Delta)*delta_vect ==1], h = np::npregbw(y2[Delta == 1] ~ x[Delta == 1])$bw, Xi = x, Yi = y2,D = Delta )

par(mar = c(5, 4, 4, 4))  
plot(ggg, kimyu, type = "l", ylab ="Equation value", xlab = "Gamma", main = "Kim and Yu Criterion")
abline(v = ggg[which.min(emp2)], lty = "dashed")
abline(v = ggg[which.min(emp1)], lty = "dashed", col = "blue")
abline(h = 0)
abline(v = fitted_kim$gamma.opt, lty = "dashed", col = "darkorchid3", lwd = 2)
abline(v = gamma, lty = "dashed", col = "red")
legend(x =6.5, y = -10.5, col = c("black","blue", "red","darkorchid3"),lty = 2, lwd = 1.3, bty = "n", legend = c("L2 minimizer", "L1 minimizer", "True Gamma", "Kim and Yu Minimizer"), cex = 0.75)


# Comparison with MSE #####

# DONT RUN, LOAD .rda FILES INSTEAD !!

# N <- c(100,300,750)
# Noise <- c(.5,1.5,3)
# 
# mse <- data.frame(
#   N = rep(N, length(Noise)),
#   Noise = rep(Noise, each =length(N)),
#   NW = 0,
#   sd.NW = 0,
#   NWadd = 0,
#   sd.NWadd = 0,
#   MOJ = 0,
#   sd.MOJ = 0,
#   KIM = 0,
#   sd.KIM = 0
# )
# 
# plot.mse <-
#   data.frame(
#     N = rep(N, length(Noise)),
#     Noise = rep(Noise, each =length(N))
#   )
# 
# niter <- 100
# 
# plot.nw <- cbind(plot.mse,matrix(NA, 9,niter))
# plot.moj <- cbind(plot.mse,matrix(NA, 9,niter))
# plot.nwadd <- cbind(plot.mse,matrix(NA, 9,niter))
# plot.kim <- cbind(plot.mse,matrix(NA, 9,niter))
# 
# for(i in 1:nrow(mse)){
#   cur_nw <- cur_moj <- cur_nw_add <- cur_kim <-  NULL
#     nn <- mse$N[i]
#     noise <- mse$Noise[i]
#   for(j in 1:niter){
#     start.time <- Sys.time()
#     # Data Generation at Each Step
#     datai <- gen_data(n = nn, noise = noise, gamma = 2)
# 
#     # Moj Estimator
#     m_moj <- moj_model(datai,x = datai$x, verbose = F, extras = T)
#     addi <- m_moj$subsample
#     m_moj_L2i <- m_moj$model$m_moj_L2m
# 
#     # NW Estimator
#     data_obsi <- subset(datai, datai$Delta == 1)
#     h_obsi <- np::npregbw(data_obsi$y2~data_obsi$x)
#     nwi <- gen_funct(datai$x, data_obsi$x, data_obsi$y2, D = 1, h = h_obsi$bw, k = 1, gamma = 0, funct = "nw")
# 
#     # NW add Estimator
#     data_addi <- rbind(data_obsi, cbind(addi, Delta = 1))
#     h_addi <- np::npregbw(data_addi$y2~data_addi$x)
#     nwaddi <- gen_funct(datai$x, data_addi$x, data_addi$y2, D = 1, h = h_addi$bw, k = 1, gamma = 0, funct = "nw")
# 
#     # Kim Estimator
#     deltai <- full_delta(m_moj, nrow(datai))
# 
#     kimnw <- kim_est(datai$x, datai$y2, D = datai$Delta, d = deltai)$kim
# 
#     # MSE
#     cur_nw <- c(cur_nw,nn^-1 * sum(abs(nwi-datai$y2)^2))
#     cur_nw_add <- c(cur_nw_add,nn^-1 * sum(abs(nwaddi-datai$y2)^2))
#     cur_moj <- c(cur_moj,n^-1 * sum(abs(m_moj_L2i-datai$y2)^2))
#     cur_kim <- c(cur_kim,n^-1 * sum(abs(kimnw-datai$y2)^2))
# 
#     end.time <- Sys.time()
#     rm(datai, m_moj, m_moj_L2i, data_obsi, h_obsi, nwi, nwaddi, data_addi, h_addi, deltai, kimnw)
#     cat(paste("Time taken for cycle", i,", rep",j,":", round(end.time - start.time),"s \n"))
#   }
#   plot.nw[i,3:(niter+2)] <- cur_nw
#   plot.moj[i,3:(niter+2)] <- cur_moj
#   plot.nwadd[i,3:(niter+2)] <- cur_nw_add
#   plot.kim[i,3:(niter+2)] <- cur_kim
# 
# 
#   mse$NW[i] <- mean(cur_nw,na.rm = T)
#   mse$MOJ[i] <- mean(cur_moj,na.rm = T)
#   mse$NWadd[i] <- mean(cur_nw_add,na.rm = T)
#   mse$KIM[i] <- mean(cur_kim,na.rm = T)
#   mse$sd.NW[i] <- sd(cur_nw,na.rm = T)
#   mse$sd.MOJ[i] <- sd(cur_moj,na.rm = T)
#   mse$sd.NWadd[i] <- sd(cur_nw_add,na.rm = T)
#   mse$sd.KIM[i] <- sd(cur_kim,na.rm = T)
#   beepr::beep()
# 
# }
# plotlist <- list(plot.nw,plot.moj,plot.nwadd, plot.kim)

load("rda files/PlotList.rda"); load("rda files/MSEtable.rda")
plot.moj <- plotlist[[2]]; plot.nw <- plotlist[[1]]; plot.nwadd <- plotlist[[3]]; plot.kim <- plotlist[[4]]

minmatX <- sapply(plotlist, function(x){
  mat <- x [,3:(niter+2)]
  foreach(j = 1:9, .combine = c) %do% {
    dens <- density(as.numeric(mat[j,]), na.rm =T )$x
    pmax(0,min(dens, na.rm  =T))
  }
})
maxmatX <- sapply(plotlist, function(x){
  mat <- x [,3:(niter+2)]
  foreach(j = 1:9, .combine = c) %do% {
    dens <- density(as.numeric(mat[j,]), na.rm =T )$x
    max(dens, na.rm  =T)
  }
})
maxmatY <- sapply(plotlist, function(x){
  mat <- x [,3:(niter+2)]
  foreach(j = 1:9, .combine = c) %do% {
    dens <- density(as.numeric(mat[j,]), na.rm =T )$y
    max(dens, na.rm  =T)
  }
})

xlims <- foreach(j = 1:9, .combine = rbind) %do%  c(min(minmatX[j,]), max(maxmatX[j,]))
ylims <- foreach(j = 1:9, .combine = rbind) %do%  c(0, max(maxmatY[j,]))
par(mfrow = c(3,3), mar = c(2,4,4,2))
for(i in 1:9){
  plot(density(as.numeric(t(plot.moj[i,3:(niter+2)])),na.rm = T), xlab = "",main = paste("Total Sample Size: ",mse$N[i], "\n Noise Variance:", mse$Noise[i]), ylab = "density", type = "l",
       xlim = xlims[i,], ylim = ylims[i,])
  lines(density(as.numeric(t(plot.nw[i,3:(niter+2)])),na.rm = T), col = "blue")
  lines(density(as.numeric(t(plot.nwadd[i,3:(niter+2)])),na.rm = T), col = "cyan")
  lines(density(as.numeric(t(plot.kim[i,3:(niter+2)])),na.rm = T), col = "darkorchid3")
  legend("topright", col = c("black","blue", "cyan", "darkorchid3"), lty = 1, legend = c("MOJ","NW", "NW+", "KIM"), bty = "n")
}

CI <- data.frame(
  N = mse$N,
  Noise = mse$Noise
)
for(i in c(3,5,7,9)){
  CI <- cbind(CI, CI =paste("(", round(mse[,i]-1.96*mse[,i+1],2), ", ",round(mse[,i]+1.96*mse[,i+1],2), ")", sep = ""))
}
colnames(CI) <- c("Size","Noise","NW","NW+","MOJ", "KIM")

CI
# Bootstrap #####

bootstrap_moj <- function(data, indices){
  d <- data[indices, ]
  model_boot <- moj_model(d, verbose = F)
  return(model_boot$m_moj_L2m) 
}
bootstrap_res <- function(data, indices){
  n <- nrow(data)
  pn <- ((log(n)^.25)/(n*0.95)^(1-.01))^.5
  m <- floor(.7 * n)
  train <- sample(n,m, replace = F)
  Dl <- data[-train,]
  ripescati <- 0
  while(sum(ripescati)==0){
    ripescati <- 0
    delta <- rbinom(n-m,1,pn)
    ripescati <- (1-Dl$Delta)*delta
  }
  out <- Dl[ripescati == 1,2]
  possible.out <- matrix(0,nrow = 30, ncol =1)
  possible.out[1:length(out),1] <- out
  return(possible.out) 
}
bootstrap_nw <- function(data, indices){
  d <- data[indices, ]
  h_bt <- np::npregbw(d$y2~d$x)
  m_nw_bootstrap <- np::npreg(h_bt)
  fitted_nw <- data.frame(x = d$x, fitted = fitted(m_nw_bootstrap))
  out <- fitted_nw[order(fitted_nw$x),2]
  return(out)
}
bootstrap_nwadd <- function(data, indices){
  d <- data[indices, ]
  n <- nrow(d)
  # SubSample
    m <- floor(0.7*n)
    pn <- ((log(n)^.25)/(n*0.95)^(1-0.01))^.5
    train <- sample(n,m,replace = F)
    Dl <- d[-train,]
    ripescati <- 0
    while(sum(ripescati)==0){
      ripescati <- 0
      delta <- rbinom(n-m,1,pn)
      ripescati <- (1-Dl$Delta)*delta
    }
    subsample <- Dl[ripescati==1,1:2]
  # NW
  dobs <- subset(d, d$Delta ==1)
  dadd <- rbind(dobs,cbind(subsample, Delta = 1))
  h_add <- np::npregbw(dadd$y2~dadd$x)
  nwadd <- gen_funct(d$x, dadd$x, dadd$y2, D = 1, h = h_add$bw, k = 1, gamma = 0, funct = "nw")
  fitted_nwadd <-  data.frame(x = d$x, fitted = nwadd)
  
  out <- fitted_nwadd[order(fitted_nwadd$x),2]
  return(out)
}
bootstrap_kim <- function(data, indices){
  d <- data[indices, ]
  d <- d[order(d$x),]
  n <- nrow(d)
  pn <- ((log(n)^.25)/(n*0.95)^(1-0.01))^.5
  m <- floor(0.7*n)
  train <- sample(n,m,replace = F)
  Dl <- d[-train,]
  ripescati <- 0
  while(sum(ripescati)==0){
    ripescati <- 0
    delta <- rbinom(n-m,1,pn)
    ripescati <- (1-Dl$Delta)*delta
  }
  usp <- list(param = list(indices = train, delta = delta))
  delta_vec <- full_delta(usp, n)
  out <- kim_est(d$x, d$y2, D = d$Delta, d = delta_vec)$kim 
  return(out)
}
alpha_perc <- 0.05  

## CI for NW

# DONT RUN FUNCTIONS, LOAD .rda INSTEAD !

# boot_nw <- boot(data = data_obs, statistic = bootstrap_nw, R = 1000); save(boot_nw, file = "BootCInw.rda")
load("rda files/BootCInw.rda")

ave_nw <- apply(boot_nw$t, 2, function(x) quantile(x, probs =.5, na.rm = T))
lower_nw <- apply(boot_nw$t, 2, function(x) quantile(x, probs = alpha_perc / 2, na.rm = T))
upper_nw <- apply(boot_nw$t, 2, function(x) quantile(x, probs = 1 - alpha_perc / 2, na.rm = T))

## CI for moj

# boot_moj <- boot(data = data, statistic = bootstrap_moj, R = 1000); save(boot_moj, file = "BootCImoj.rda")
load("rda files/BootCImoj.rda")

ave_moj <- apply(boot_moj$t, 2, function(x) quantile(x, probs = .5, na.rm = T))
lower_moj <- apply(boot_moj$t, 2, function(x) quantile(x, probs = alpha_perc / 2, na.rm = T))
upper_moj <- apply(boot_moj$t, 2, function(x) quantile(x, probs = 1 - alpha_perc / 2, na.rm = T))

## CI for NW+

# boot_nwadd <- boot(data = data, statistic = bootstrap_nwadd, R = 500); save(boot_nwadd, file = "BootCInwad.rda")
load("rda files/bootCInwad.rda")

ave_nwadd <- apply(boot_nwadd$t, 2, function(x) quantile(x, probs = .5, na.rm = T))
lower_nwadd <- apply(boot_nwadd$t, 2, function(x) quantile(x, probs = alpha_perc / 2, na.rm = T))
upper_nwadd <- apply(boot_nwadd$t, 2, function(x) quantile(x, probs = 1 - alpha_perc / 2, na.rm = T))

## CI for Kim+

# boot_kim <- boot(data = data, statistic = bootstrap_kim, R = 1000); save(boot_kim, file = "BootCIkim.rda")
load("rda files/bootCIkim.rda")

ave_kim <- apply(boot_kim$t, 2, function(x) quantile(x, probs = .5, na.rm = T))
lower_kim <- apply(boot_kim$t, 2, function(x) quantile(x, probs = alpha_perc / 2, na.rm = T))
upper_kim <- apply(boot_kim$t, 2, function(x) quantile(x, probs = 1 - alpha_perc / 2, na.rm = T))


## Boostrap for Resampled

boot_res <- boot(data = data, statistic = bootstrap_res, R = 2000)
boot_fres <- as.vector(t(boot_res$t))[as.vector(t(boot_res$t))!=0]
boot_matres <- data.frame(
  x = data$x[match(boot_fres, data$y2)],
  y2 = boot_fres
)
boot_contres <- MASS::kde2d(boot_matres$x, boot_matres$y2)


## CI Plot  ####
# par(mfrow = c(2,2), mar = c(3,2.7,3,2))
par(mfrow = c(1,1), mar = c(3,2.7,3,2))
with(data, plot(x,y2, cex = .15,pch = 19, col = scales::alpha((2-Delta),.5)))
lines(fitted_nw , col = "royalblue3", lwd = 2)
lines(fitted_nw$x,fitted(smooth.spline(fitted_nw$x,lower_nw,spar=.1)), col = "blue", lty = "dashed", lwd = 2)
lines(fitted_nw$x,fitted(smooth.spline(fitted_nw$x,upper_nw, spar = .1)), col = "blue", lty = "dashed", lwd = 2)
lines(df_moj_L2m$model$x, fitted(smooth.spline(df_moj_L2m$model$x,lower_moj)), col = scales::alpha("forestgreen",.6), lty = "dashed", lwd = .5)
lines(df_moj_L2m$model$x, fitted(smooth.spline(df_moj_L2m$model$x,upper_moj)), col = scales::alpha("forestgreen",.6), lty = "dashed", lwd = .5)


with(data, plot(x,y2, cex = .15,pch = 19, col = scales::alpha((2-Delta),.5)))
lines(df_moj_L2m$model, col = "springgreen3", lwd = 2)
lines(df_moj_L2m$model$x, fitted(smooth.spline(df_moj_L2m$model$x,lower_moj)), col = "forestgreen", lty = "dashed", lwd = 2)
lines(df_moj_L2m$model$x, fitted(smooth.spline(df_moj_L2m$model$x,upper_moj)), col = "forestgreen", lty = "dashed", lwd = 2)
lines(fitted_nw$x,fitted(smooth.spline(fitted_nw$x,lower_nw,spar=.1)), col = scales::alpha("blue",0.6), lty = "dashed", lwd = .5)
lines(fitted_nw$x,fitted(smooth.spline(fitted_nw$x,upper_nw, spar = .1)), col = scales::alpha("blue",0.6), lty = "dashed", lwd = .5)

with(data, plot(x,y2, cex = .15,pch = 19, col = scales::alpha((2-Delta),.5)))
lines(fitted_nw_res, col = "darkslategray3", lwd = 2)
lines(sort(data$x), fitted(smooth.spline(sort(data$x),lower_nwadd)), col = "darkslategray3", lty = "dashed", lwd = 2)
lines(sort(data$x), fitted(smooth.spline(sort(data$x),upper_nwadd)), col = "darkslategray3", lty = "dashed", lwd = 2)
lines(fitted_nw$x,fitted(smooth.spline(fitted_nw$x,lower_nw,spar=.1)), col = scales::alpha("blue",0.6), lty = "dashed", lwd = .5)
lines(fitted_nw$x,fitted(smooth.spline(fitted_nw$x,upper_nw, spar = .1)), col = scales::alpha("blue",0.6), lty = "dashed", lwd = .5)
lines(df_moj_L2m$model$x, fitted(smooth.spline(df_moj_L2m$model$x,lower_moj)), col = scales::alpha("forestgreen",.6), lty = "dashed", lwd = .5)
lines(df_moj_L2m$model$x, fitted(smooth.spline(df_moj_L2m$model$x,upper_moj)), col = scales::alpha("forestgreen",.6), lty = "dashed", lwd = .5)

with(data, plot(x,y2, cex = .15,pch = 19, col = scales::alpha((2-Delta),.5)))
lines(fitted_kim$model, col = "darkorchid3", lwd = 2)
lines(data$x, fitted(smooth.spline(data$x,lower_kim)), col = "darkorchid3", lty = "dashed", lwd = 2)
lines(data$x, fitted(smooth.spline(data$x,upper_kim)), col = "darkorchid3", lty = "dashed", lwd = 2)
lines(fitted_nw$x,fitted(smooth.spline(fitted_nw$x,lower_nw,spar=.1)), col = scales::alpha("blue",0.6), lty = "dashed", lwd = .5)
lines(fitted_nw$x,fitted(smooth.spline(fitted_nw$x,upper_nw, spar = .1)), col = scales::alpha("blue",0.6), lty = "dashed", lwd = .5)
lines(df_moj_L2m$model$x, fitted(smooth.spline(df_moj_L2m$model$x,lower_moj)), col = scales::alpha("forestgreen",.6), lty = "dashed", lwd = .5)
lines(df_moj_L2m$model$x, fitted(smooth.spline(df_moj_L2m$model$x,upper_moj)), col = scales::alpha("forestgreen",.6), lty = "dashed", lwd = .5)

## Resampled Plot
layout(matrix(c(1,1,2), nrow = 1, ncol = 3))
contour(boot_contres, col = "navyblue", labcex = .7, lwd = 1.3, xlim = c(-10,10),ylim= c(-3.5,7), main = "Resampled Dneisty")
# image(boot_contres,col = scales::alpha(hcl.colors(50),.5))
with(data, points(x,y2, cex = .4,pch = 19, col = (2-Delta)))
hist(apply(t(boot_res$t),2, function(x) sum(x!=0)), breaks = 17, xlab = "Resampled", main = "Size of Follow-up Sample", col = "whitesmoke", ylab = "")
     

