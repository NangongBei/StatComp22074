## ---- results='asis'----------------------------------------------------------
Speed <- cars$speed
Dist <- cars$dist
print(xtable::xtable(head(cars)),type = "html")

## -----------------------------------------------------------------------------
summary(cars)

## -----------------------------------------------------------------------------
lm.cars <- lm(Dist~Speed)
summary(lm.cars)$coef

## -----------------------------------------------------------------------------
#par(mfrow=c(2,2))
plot(lm.cars)

## ----eval=TRUE----------------------------------------------------------------
rdistPareto = function(n = 1e3,a = 2,b = 2,picture = FALSE){
  u <- runif(n) # Generate uniformly distributed random numbers
  x <- b/(u)^(1/a) # calculate F^(-1)(U')
  if(picture){
    # Draw a histogram
    hist(x, prob = TRUE, breaks = 'scott',main = expression(f(x)==a*b^a/x^(a+1),x>=b))
    # Sampling and plotting probability density curves
    sampling <- seq(b,max(x),.1)
    lines(sampling,a*b^a/sampling^(a+1))
  }
  return(x)
}
x <- rdistPareto(1e3,2,2,TRUE)
rm(list=ls())

## ----eval=TRUE----------------------------------------------------------------
rdistBETA = function(n = 1e3,alpha = 3,beta = 2,picture = FALSE){
  sample_total <- m <- 0 # draw the mth sample and Output sampling total
  y <- numeric(n) # Record the generated random numbers
  while(m < n)
  {
    u <- runif(1) # Generate random numbers from U(0,1)
    sample_total <- sample_total + 1
    x <- runif(1) # Generate random numbers from g(x)
    if(u <= x^(alpha-1)*(1-x)^(beta-1)) # Determine whether u <= rho(x)
    {
      # accept x and record it
      m <- m + 1
      y[m] <- x
    }
  }
  if(picture){
    # Draw a histogram
    hist(y, prob = TRUE, breaks = 'scott',main = expression(f(x)==x^(alpha-1)*(1-x)^(beta-1)/B(alpha,beta)))
    # Sampling and plotting probability density curves
    sampling <- seq(0,1,.01)
    lines(sampling,sampling^(alpha-1)*(1-sampling)^(beta-1)/beta(alpha,beta))
  }
  return(y)
}
y <- rdistBETA(1e3,3,2,TRUE)
rm(list=ls())

## ----eval=TRUE----------------------------------------------------------------
rdistGammaExp = function(n = 1e3,r = 4,beta = 2,picture = FALSE){
  lambda_sampling <- rgamma(n,r,beta) # Generate gamma distributed random numbers
  y <- rexp(lambda_sampling) # the length of x and lambda is n = 1000
  if(picture){
    # Draw a histogram to show the random sample
    hist(y, prob = TRUE, breaks = 'scott',main = expression(y %~% Exp(lambda %~% Gamma(r,beta))))
    # Sampling and plotting probability density curves
    sampling <- seq(0,10,.01)
    lines(sampling,r*beta^r/(beta+sampling)^(r+1))
  }
  return(y)
}
x <- rdistGammaExp(1e3,4,2,TRUE)
rm(list=ls())

## ----eval=TRUE----------------------------------------------------------------
quick_sort<-function(x){
  num<-length(x)
  if(num==0||num==1){return(x)
  }else{
    a<-x[1]
    y<-x[-1]
    lower<-y[y<a]
    upper<-y[y>=a]
    return(c(quick_sort(lower),a,quick_sort(upper)))}
}

## ----eval=TRUE----------------------------------------------------------------
HW2_Ex1 <- function(n=100){
  set.seed(1)
  sample_n <- c(1,2,4,6,8) *1e4 # The array size
  test_time_sample <- numeric(length(sample_n)) # record spending time
  # calculate time
  for(i in 1:length(sample_n)){
    test_number<-sample(1:sample_n[i])
    test_time_sample[i] <- system.time(for(j in 1:n) quick_sort(test_number))[1]/n
  }
  cat('The average time:',test_time_sample)
  t <- sample_n*log(sample_n)
  t_test_time_lm <- lm(t~test_time_sample) # linear regression
  summary(t_test_time_lm)
  # graph the scatter plot and regression line
  plot(test_time_sample,t,main="Scatter plots and regression line")
  abline(t_test_time_lm,col='red')
  legend("topleft",legend="Regression line",col="red",lty=1)
}

## ----eval=TRUE----------------------------------------------------------------
HW2_Ex1(100)
rm(list = ls())

## ----eval=TRUE----------------------------------------------------------------
HW2_Q_5_7 <- function(n=1e3){
  set.seed(1)
  u <- runif(2*n) # Take the first n random numbers for Uniform distribution
  u1 <- u[1:n] # Take the first n random numbers
  theta_MC <- mean(exp(u)) # simple MC estimator
  theta <- mean((exp(u1)+exp(1-u1))/2) # antithetic variate approach
  cat('Simple MC estimator:',theta_MC,'\n')
  cat('Antithetic Variate estimator:',theta,'\n')
  var_theta <- var((exp(u1)+exp(1-u1))) / (4*n) # variance of antithetic variate approach
  var_theta_MC <- var(exp(u)) / (2*n) # variance of simple MC estimator
  cat('The percent reduction in variance:',1 - var_theta / var_theta_MC)
}

## ----eval=TRUE----------------------------------------------------------------
HW2_Q_5_7(1e3)
rm(list = ls())

## -----------------------------------------------------------------------------
x <- seq(1,10,0.01) #sampling
y <- x^2 * exp(- x^2/2) / sqrt(2*pi) # function of g(x)
plot(x,y,col = 'black',xlim = c(1,10),type = 'l',ylab = 'y',main = 'Graph of g(x)')

## -----------------------------------------------------------------------------
set.seed(1)
n <- 1e4 # random numbers size
x1 <- -log(1-runif(n)) + 1 # X1~F1
x2 <- -2 * log(1-runif(n)) + 1 # X2~F2

## -----------------------------------------------------------------------------
theta_hat1 <- x1^2 * exp(- x1^2/2) / sqrt(2*pi) / exp(- (x1 - 1))
theta_hat2 <- x2^2 * exp(- x2^2/2) / sqrt(2*pi) / (exp(- (x2 - 1)/2) / 2)
var_theta1 <- var(theta_hat1) / n # variance of f1 important sampling
var_theta2 <- var(theta_hat2) / n # variance of f2 important sampling
print(c(mean(theta_hat1),mean(theta_hat2)))
print(c(var_theta1,var_theta2))

## -----------------------------------------------------------------------------
x <- seq(1,10,0.01)
y <- x^2 * exp(- x^2/2) / sqrt(2*pi) # g(x)
f1_ing <- exp(- (x - 1)) # f1(x)
f2_ing <- exp(- (x - 1)/2) / 2 # f2(x)
plot(x,f1_ing,col = 'red',xlim = c(1,10),type = 'l',ylab = 'y',main = 'Graph of three functions')
lines(x,f2_ing,col = 'blue')
lines(x,y,col = 'black')
legend('topright',c('y','f1','f2'),col = c('black','red','blue'),lty=1)

## -----------------------------------------------------------------------------
plot(x,y / f2_ing,col = 'blue',xlim = c(1,10),type = 'l',ylab = 'ratio',main = 'Graph of g(x)/f(x)')
lines(x,y / f1_ing, col = 'red')
legend('topright',c('g(x)/f1(x)','g(x)/f2(x)'),col = c('red','blue'),lty=1)

## -----------------------------------------------------------------------------
set.seed(1)
k <- 5 # numbers of layers
n <- 10000# total sampling
m <- n / k # sampling number in each layer
theta_hat <- theta_var <- 0
for(j in 1:k)
{
  x <- -log(exp(-(j-1)/5)-(exp(-(j-1)/5)-exp(-j/5))*runif(m)) # generate xj~Fj
  theta <- (exp(-(j-1)/5)-exp(-j/5)) / (1 + x^2) # calculate g(Xj)/fj(Xj)
  theta_hat <- theta_hat + mean(theta) # estimate theta
  theta_var <- theta_var + var(theta)/m # estimate varience of theta
}
theta_hat
theta_var
rm(list=ls())

## -----------------------------------------------------------------------------
# generate and save data
set.seed(1)
sampling <- function(n = 1e2, m = 1e3, mu = 0, sigma = 1){
  y <- rnorm(n*m,mean = mu, sd = sigma) # generate normal random sample
  x <- exp(y) # generate log normal sample from normal random sample
  return(x) # Clean up the memory
}

## -----------------------------------------------------------------------------
# analyzing data
analyzing <- function(n = 1e2, m = 1e3,alpha=0.05,mu = mu){
  x_all <- sampling(n,m,mu,sigma=1)
  result <- numeric(m)
  # calculate confidence interval
  for(i in 1:m){
    x <- x_all[((i-1)*n+1):(i*n)]
    lower <- mean(log(x))-qt(1-alpha/2,length(x)-1)*sqrt(var(log(x)))/sqrt(length(x))
    upper <- mean(log(x))+qt(1-alpha/2,length(x)-1)*sqrt(var(log(x)))/sqrt(length(x))
    result[i] <- (lower < mu) * (upper > mu) # judge whether mu in confidence interval
  }
  result
}

## -----------------------------------------------------------------------------
# result and summary
showing <- function(n = 1e2, m = 1e3, mu = 0, alpha = 0.05){
  sampling(n=n,m=m,mu = mu)
  result <- analyzing(n=n,m=m,alpha = alpha, mu = mu)
  print(mean(result))
}

## -----------------------------------------------------------------------------
# test
set.seed(1)
showing(n = 1e2, m = 1e3, mu = 0, alpha = 0.05)
showing(n = 1e2, m = 1e3, mu = 10, alpha = 0.05)
showing(n = 1e2, m = 1e3, mu = 100, alpha = 0.05)
rm(list=ls())

## -----------------------------------------------------------------------------
# generate and save data
data_generate <- function(nx=10,ny=10,sigma1=1,sigma2=2,m=1e3){
  x <- rnorm(nx*m, 0, sigma1)
  y <- rnorm(ny*m, 0, sigma2)
  return(list(x=x,y=y))
}

## -----------------------------------------------------------------------------
# count5test statistic
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
power5test <- function(nx=10, ny=10,m=1e3,data){
  # Import data
  x_all <- data$x
  y_all <- data$y
  result <- numeric(m)
  for(i in 1:m){
    x <- x_all[((i-1)*nx+1):(i*nx)]
    y <- y_all[((i-1)*ny+1):(i*ny)]
    result[i] <- count5test(x, y) # count five test
  }
  print(mean(result))
}

## -----------------------------------------------------------------------------
# count5test statistic
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
  # return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}
power5test <- function(nx=10, ny=10,m=1e3,data){
  # Import data
  x_all <- data$x
  y_all <- data$y
  result <- numeric(m)
  for(i in 1:m){
    x <- x_all[((i-1)*nx+1):(i*nx)]
    y <- y_all[((i-1)*ny+1):(i*ny)]
    result[i] <- count5test(x, y) # count five test
  }
  print(mean(result))
}

## -----------------------------------------------------------------------------
# F-test
powerFtest <- function(nx=10, ny=10,m=1e3,alpha = 0.055,data){
  # Import data
  x_all <- data$x
  y_all <- data$y
  result <- numeric(m)
  for(i in 1:m){
    x <- x_all[((i-1)*nx+1):(i*nx)]
    y <- y_all[((i-1)*ny+1):(i*ny)]
    # Determine whether in the reject field
    result[i] <- 
      1 - as.integer((var(x)/var(y) < qf(1-alpha/2,nx-1,ny-1))*
                       (var(x)/var(y) > qf(alpha/2,nx-1,ny-1)))
  }
  print(mean(result))
}

## -----------------------------------------------------------------------------
# result and summary
power <- function(nx=10, ny=10,m=1e3,alpha = 0.055,sigma1=1,sigma2=10){
  data <- data_generate(nx=nx, ny=ny,sigma1=sigma1,sigma2=sigma2,m=m)
  print('Power of Count Five Test:')
  power5test(nx=nx, ny=ny,m=m,data=data)
  print('Power of F Test:')
  powerFtest(nx=nx, ny=ny,m=m,alpha = alpha,data=data)
}

## -----------------------------------------------------------------------------
# test
set.seed(1)
power(nx=20, ny=20,sigma1=1,sigma2=1.5) #small size
power(nx=80, ny=80,sigma1=1,sigma2=1.5) #median size
power(nx=200, ny=200,sigma1=1,sigma2=1.5) #large size
rm(list=ls())

## -----------------------------------------------------------------------------
rm(list = ls()) # Clear memory
library(boot);
# get data
data <- aircondit$hours
set.seed(1)
B <- 1e4
lambda_star <- numeric(B)
# generate data
generate_data <- function(data = data,B = B){
  for(b in 1:B){
    datastar <- sample(data,replace=TRUE)
    lambda_star[b] <- 1 / mean(datastar)
  }
  return(lambda_star) #save data
}
# calculate data
result_calculate <- function(data = data,B=B){
  lambda_star <- generate_data(data,B)
  bias <- mean(lambda_star) - 1/mean(data) # calculate bias
  se_boot <- sd(lambda_star) # calculate bias standard error
  list(bias = bias, se_boot = se_boot)
}
# show the result
show_result <- function(data = data,B = B){
  result_calculate(data,B)
}

## -----------------------------------------------------------------------------
#show
show_result(data,B)
rm(list=ls()) # Clear memory

## -----------------------------------------------------------------------------
library(boot)
set.seed(1)
# get data
data <- aircondit$hours
# Estimate function
boot_theta <- function(datastar,i) mean(datastar[i])
B <- 1e4
# generate data
obj <- boot(data = data, statistic = boot_theta,R = B)
# analyse data
CI <- boot.ci(obj,type=c("norm","basic","perc","bca"))
CI_norm <- CI$norm[2:3]
CI_basic <- CI$basic[4:5]
CI_perc <- CI$percent[4:5]
CI_bca <- CI$bca[4:5]
#show the result
cat('theta_hat=',obj$t0)
cat('CI_norm = (', CI_norm, ')\n',
    'CI_basic = (', CI_basic, ')\n',
    'CI_perc = (', CI_perc, ')\n',
    'CI_bca = (', CI_bca, ')\n')
rm(list=ls()) # clear memory

## ----eval=TRUE----------------------------------------------------------------
library(boot)
set.seed(1)
boot_mean <- function(datastar,i) mean(datastar[i])
B <- 1e3;n <- 20
CI_norm <- CI_basic <- CI_perc  <- matrix(NA,B,2)
for(i in 1:B){
  # generate data
  data <- rnorm(n)
  obj <- boot(data = data, statistic = boot_mean,R = B)
  # analyse data
  CI <- boot.ci(obj,type=c("norm","basic","perc"))
  CI_norm[i,] <- CI$norm[2:3]
  CI_basic[i,] <- CI$basic[4:5]
  CI_perc[i,] <- CI$percent[4:5]
}
coverage_norm <- 1 - mean(CI_norm[,1] > 0)-mean(CI_norm[,2] < 0);
coverage_basic <- 1 -  mean(CI_basic[,1] > 0)-mean(CI_basic[,2] < 0);
coverage_perc <- 1 - mean(CI_perc[,1] > 0)-mean(CI_perc[,2] < 0);
left_miss_norm <- mean(CI_norm[,1] > 0); right_miss_norm <- mean(CI_norm[,2] < 0); 
left_miss_basic <- mean(CI_basic[,1] > 0); right_miss_basic <- mean(CI_basic[,2] < 0); 
left_miss_perc <- mean(CI_perc[,1] > 0); right_miss_perc <- mean(CI_perc[,2] < 0); 
# show result
cat('left_miss_norm = ',left_miss_norm, 'right_miss_norm = ',right_miss_norm,'coverage_norm = ',coverage_norm)
cat('left_miss_basic = ',left_miss_basic, 'right_miss_basic = ',right_miss_basic,'coverage_basic = ',coverage_basic)
cat('left_miss_perc = ',left_miss_perc, 'right_miss_perc = ',right_miss_perc,'coverage_perc = ',coverage_perc)
pander::pander(data.frame("method"=c("norm","basic","percentile"),"coverage probabilities"=c(coverage_norm,coverage_basic,coverage_perc),"left miss proportion"=c(left_miss_norm,left_miss_basic,left_miss_perc),"right miss proportion"=c(right_miss_norm,right_miss_basic,right_miss_perc)))
rm(list = ls()) # Clear memory

## -----------------------------------------------------------------------------
rm(list=ls()) # clear memory
library(bootstrap)
# read data
get_data <- function(){
  scor
}
# calculate jackknife estimate of theta
theta.jacknife <- function(scor){
  n <- length(scor[,1])
  theta_jack <- numeric(n)
  for(i in 1:n)
  {
    lambda_jack <- eigen(var(scor[-i,]))$values
    theta_jack[i] <- max(lambda_jack) / sum(lambda_jack)
  }
  theta_jack
}
# show the result
Result <- function(scor,func = theta.jacknife){
  n <- length(scor[,1])
  lambda_hat <- eigen(var(scor))$values
  theta_hat <- max(lambda_hat) / sum(lambda_hat) # calculate theta_hat
  theta_jack <- func(scor)
  bias <- (n - 1)*(mean(theta_jack) - theta_hat)
  se <- sqrt((n-1)*mean((theta_jack - mean(theta_jack))^2))
  cat('bias =',bias,'\n')
  cat('standard error =',se)
}

## -----------------------------------------------------------------------------
data <- get_data()
Result(data,theta.jacknife)
rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
library(DAAG)
#get data
get_data <- function(){
  ironslag
}
leave_two_out_estimate <- function(ironslag){
  n <- length(ironslag[,1])
  e1 <- e2 <- e3 <- e4 <- 0
  for(i in 1:(n-1)){
    for(k in (i+1):n){
      j <- c(i,k)
      y1 <- ironslag[-j,]$magnetic
      x1 <- ironslag[-j,]$chemical
      y2 <- ironslag[j,]$magnetic
      x2 <- ironslag[j,]$chemical
      # Model 1
      J1 <- lm(y1~x1)
      yhat1 <- J1$coef[1] + J1$coef[2] * x2
      e1 <- e1 + mean((y2 - yhat1)^2)
      # Model 2
      J2 <- lm(y1~x1+I(x1^2))
      yhat2 <- J2$coef[1] + J2$coef[2] * x2 + J2$coef[3] * x2 ^ 2
      e2 <- e2 + mean((y2 - yhat2)^2)
      # Model 3
      J3 <- lm(log(y1)~x1)
      yhat3 <- J3$coef[1] + J3$coef[2] * x2
      e3 <- e3 + mean((y2 - exp(yhat3))^2)
      # Model 4
      J4 <- lm(log(y1)~log(x1))
      yhat4 <- J4$coef[1] + J4$coef[2] * log(x2)
      e4 <- e4 + mean((y2 - exp(yhat4))^2)
    }
  }
  return(c(e1,e2,e3,e4) * 2/ ((n-1)*n))
}
# show result
Result <- function(){
  ironslag <- get_data()
  result <- leave_two_out_estimate(ironslag)
  print(result)
}

## -----------------------------------------------------------------------------
Result()
rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
set.seed(1)
# get data
get_data <- function(n = 50, dependence = 1){
  x <- rnorm(n)
  if(dependence == 1) y <- 5*x+rnorm(n)
  else if(dependence == 0) y <- rnorm(n)
  z <- cbind(x,y)
}
# simulate
Spearman_permutation <- function(B = 1e3, z){
  x <- z[,1];y <- z[,2];
  n <- length(x); K <- length(x)
  spearman0 <- cor.test(x,y,method = 'spearman')$estimate
  spearman <- numeric(B)
  for(i in 1:B){
    k <- sample(1:K,n,replace = FALSE)
    x1 <- z[,1];y1 <- z[k,2]
    spearman[i] <- cor.test(x1,y1,method = 'spearman')$estimate
  }
  p_value <- mean(abs(c(spearman0,spearman))>=abs(spearman0))
  cat('Spearman_permutation p-value =',p_value,'\n')
  cat('Spearman p-value =',cor.test(x,y,method = 'spearman')$p.value,'\n')
}
# result
Result <- function(n = 50, dependence = 1,B = 1e3){
  z <- get_data(n,dependence)
  Spearman_permutation(B,z)
}

## -----------------------------------------------------------------------------
Result(n = 50, dependence = 1,B = 1e3)
Result(n = 50, dependence = 0,B = 1e3)
rm(list=ls()) # clear memory

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(2022)
#generate data
generate_data <- function(Xt,sigma){
  Y <- rnorm(1,Xt,sigma)
  fy <- exp(-abs(Y))/2
  fx <- exp(-abs(Xt))/2
  alpha <- min(1,fy/fx) #Acceptance Probability
  U <- runif(1)
  if(U<=alpha) {list(X=Y,accept=1)}
  else list(X=Xt,accept=0)
}
# one chain to converge
one_chain <- function(X0,sigma){
  N <- 1
  m <- length(X0)
  X <- phi <- Y <- acceptance_prob <- numeric(N*m)
  dim(X) <- dim(phi) <- c(N,m)
  X[N,] <- phi[N,] <- X0
  R <- 1e3
  while(R>=1.2){
    for(k in 1:m){
      result <- generate_data(X[N,k],sigma)
      Y[k] <- result$X
      acceptance_prob[k] <- acceptance_prob[k] + result$accept
  }
    X <- rbind(X,Y)
    phi <- rbind(phi,colMeans(X))
    N <- N+1
    B <- N/(m-1)*sum((colMeans(phi)-mean(phi))^2)
    W <- mean((t(phi) - colMeans(phi))^2)
    var.hat <- N/(N - 1)*W + B/N
    R <- var.hat/W
  }
  plot(phi[,1],type = 'l',ylab = bquote(phi),ylim = c(-max(X0),max(X0)),col = 2,xlab = bquote(sigma==.(sigma)))
  for(k in 2:m){
    lines(phi[,k],col=k+1)
  }
  legend('topright',legend = X0,col=1:m+1,lty=1,title = 'X0')
  list(X=rowMeans(X),N=N,acceptance_prob = mean(acceptance_prob/(N-1)),R=R,phi=phi)
}
# show result
show_result <- function(X0,sigma){
  N <- R <- accept <- numeric(length(sigma))
  for(k in 1:length(sigma)){
    result_X <- one_chain(X0,sigma[k])
    X <- result_X$X
    accept[k] <- result_X$acceptance_prob
    N[k] <- result_X$N
    R[k] <- result_X$R
    plot(X,type = 'l',ylab = 'x',xlab = bquote(sigma==.(sigma[k])),main = bquote(Acceptance_rate ==.(accept[k])))
  }
  list(accept=accept,N=N,R=R)
}
# result
X0 <- c(-10,-5,0,5,10);sigma <- c(0.5,1,5,10)
result_X_sigma <- show_result(X0,sigma)
pander::pander(data.frame('sigma'=sigma,'Acceptance_Probability'=result_X_sigma$accept,'N'=result_X_sigma$N,'R'=result_X_sigma$R))
rm(list=ls())

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(1)
# generate data
gibbs_sample <- function(X0,Y0,M = 1e3){
  N <- 1
  m <- length(X0)
  if(m != length(Y0)){
    print('Wrong: length(X0) is not equal to length(Y0)!')
    return(0)
  }
  X <- phi_X <- phi_Y <- Y <- numeric(N*m)
  dim(X) <- dim(Y) <- dim(phi_X) <- dim(phi_Y) <- c(N,m)
  X[1,] <- phi_X[1,] <- X0; Y[1,] <- phi_Y[1,] <- Y0
  RX <- RY <- 1e3
  while(RX>=1.2|RY>=1.2){
    # generate data
    NX <- rnorm(m,0.9*Y[N,],sqrt(0.19))
    X <- rbind(X,NX)
    NY <- rnorm(m,0.9*X[N+1,],sqrt(0.19))
    Y <- rbind(Y,NY)
    N <- N+1
    # calculate R_x and R_y
    phi_X <- rbind(phi_X,colMeans(X))
    phi_Y <- rbind(phi_Y,colMeans(Y))
    BX <- N/(m-1)*sum((colMeans(phi_X)-mean(phi_X))^2)
    WX <- mean((t(phi_X) - colMeans(phi_X))^2)
    var.hat.X <- N/(N - 1)*WX + BX/N
    RX <- var.hat.X/WX
    BY <- N/(m-1)*sum((colMeans(phi_Y)-mean(phi_Y))^2)
    WY <- mean((t(phi_Y) - colMeans(phi_Y))^2)
    var.hat.Y <- N/(N - 1)*WY + BY/N
    RY <- var.hat.Y/WY
  }
  for(i in 1:M){
    # generate data
    NX <- rnorm(m,0.9*Y[N-1+i,],sqrt(0.19))
    X <- rbind(X,NX)
    NY <- rnorm(m,0.9*X[N+i,],sqrt(0.19))
    Y <- rbind(Y,NY)
    phi_X <- rbind(phi_X,colMeans(X))
    phi_Y <- rbind(phi_Y,colMeans(Y))
  }
  plot(phi_X[,1],type = 'l',ylab = bquote(phi_X),ylim = c(-max(X0),max(X0)),col = 2,xlab = 'plot for X')
  for(k in 2:m){
    lines(phi_X[,k],col=k+1)
  }
  legend('topright',legend = X0,col=(1:m)+1,lty=1,title = 'X0')
  plot(phi_Y[,1],type = 'l',ylab = bquote(phi_Y),ylim = c(-max(Y0),max(Y0)),col = 2,xlab = 'plot for Y')
  for(k in 2:m){
    lines(phi_Y[,k],col=k+1)
  }
  legend('topright',legend = Y0,col=(1:m)+1,lty=1,title = 'Y0')
  list(X=rowMeans(X),Y=rowMeans(Y),N=N,RX=RX,phi_X=phi_X,RY=RY,phi_Y=phi_Y)
}
library('lmtest')
# show result
gibbs_result <- function(X0,Y0,M=1e3){
  result <- gibbs_sample(X0,Y0,M)
  X <- result$X[(result$N+1):(result$N+M)]
  Y <- result$Y[(result$N+1):(result$N+M)]
  plot(X,type = 'l',ylab = 'X',main = 'The chain of X after discarding')
  plot(Y,type = 'l',ylab = 'Y',main = 'The chain of Y after discarding')
  cat('N =',result$N,'\t','RX =',result$RX,'\t','RY =',result$RY,'\n')
  lm.YX <- lm(Y~X)
  cat('Fitted linear model:','\n')
  print(summary(lm.YX))
  # test normality
  cat('\n','Test normality:','\n')
  print(shapiro.test(lm.YX$residuals))
  # Breusch-Pagan Test
  cat('\n','Breusch-Pagan Test:','\n')
  print(bptest(lm.YX))
}
X0 <- c(-5,-1,0,1,5);Y0 <- c(-5,-1,0,1,5)
gibbs_result(X0,Y0,M=1e3)
rm(list=ls())

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(2022)
# generate data
generate_data <- function(N=1e3,alpha,beta,gamma_x=1){
  am <- ay <- 0
  X <- rnorm(N)
  M <- am + alpha*X + rnorm(N)
  Y <- ay + beta*M + gamma_x*X + rnorm(N)
  data_XYM <- data.frame(X=X,M=M,Y=Y)
  return(data_XYM)
}

## -----------------------------------------------------------------------------
library(boot)
library(mediation)
# permutation get p-value
# Condition 1
Permutation_pvalue1 <- function(data_XYM){
  Tn_med <- function(data_XYM,ix){
    data <- data.frame(X=data_XYM$X[ix],M=data_XYM$M,Y=data_XYM$Y)
    model.m <- lm(M ~ X, data = data)
    model.y <- lm(Y ~ M + X, data = data)
    out.1 <- mediate(model.m, model.y, sims = 100, treat = "X",mediator = "M")
    return(out.1$d.avg.p)
  }
  boot.med <- boot(data_XYM,Tn_med,R=1e2,sim = "permutation") 
  ps <- c(boot.med$t0,boot.med$t)
  p.value <- mean(ps<=ps[1])
  return(p.value)
}
# Condition 2
Permutation_pvalue2 <- function(data_XYM){
  Tn_med <- function(data_XYM,ix){
    data <- data.frame(X=data_XYM$X,M=data_XYM$M,Y=data_XYM$Y[ix])
    model.m <- lm(M ~ X, data = data)
    model.y <- lm(Y ~ M + X, data = data)
    out.1 <- mediate(model.m, model.y, sims = 100, treat = "X",mediator = "M")
    return(out.1$d.avg.p)
  }
  boot.med <- boot(data_XYM,Tn_med,R=1e2,sim = "permutation") 
  ps <- c(boot.med$t0,boot.med$t)
  p.value <- mean(ps<=ps[1])
  return(p.value)
}
# Condition 3
Permutation_pvalue3 <- function(data_XYM){
  Tn_med <- function(data_XYM,ix){
    data <- data.frame(X=data_XYM$X,M=data_XYM$M[ix],Y=data_XYM$Y)
    model.m <- lm(M ~ X, data = data)
    model.y <- lm(Y ~ M + X, data = data)
    out.1 <- mediate(model.m, model.y, sims = 100, treat = "X",mediator = "M")
    return(out.1$d.avg.p)
  }
  boot.med <- boot(data_XYM,Tn_med,R=1e2,sim = "permutation") 
  ps <- c(boot.med$t0,boot.med$t)
  p.value <- mean(ps<=ps[1])
  return(p.value)
}

## ----eval=TRUE----------------------------------------------------------------
N <- 1e2
gamma_x <- 1
alpha <- c(0,0,1)
beta <- c(0,1,0)
p_value1 <- numeric(length(alpha))
p_value2 <- numeric(length(alpha))
p_value3 <- numeric(length(alpha))
for(i in 1:length(alpha)){
  data_XYM <- generate_data(N,alpha[i],beta[i],gamma_x)
  p_value1[i] <- Permutation_pvalue1(data_XYM)
  p_value2[i] <- Permutation_pvalue2(data_XYM)
  p_value3[i] <- Permutation_pvalue3(data_XYM)
}
pander::pander(data.frame("alpha"=alpha,"beta"=beta,"Condition1"=p_value1,"Condition2"=p_value2,"Condition3"=p_value3))
rm(list=ls())

## -----------------------------------------------------------------------------
rm(list=ls())
set.seed(1)
# solve function
Find_alpha <- function(N=1e6,b1=0,b2=1,b3=-1,f0=0.1){
  x1 <- rpois(N,1);x2 <- rexp(N);
  x3 <- sample(0:1,N,replace = TRUE)
  G <- function(alpha){
    P <- 1 / (1+exp(-alpha-b1*x1-b2*x2-b3*x3))
    return(mean(P) - f0)
  }
  solution <- uniroot(G,c(-50,0))
  # return alpha
  return(solution$root)
}

## -----------------------------------------------------------------------------
# input N,b1,b2,b3,f0
N <- 1e6
b1 <- 0;b2 <- 1;b3 <- -1
f0 <- c(0.1,0.01,0.001,0.0001)
n <- length(f0)
alpha <- numeric(n)
# output alpha
for(i in 1:n){
  alpha[i] <- Find_alpha(N,b1,b2,b3,f0[i])
}
cbind(alpha,f0)

## -----------------------------------------------------------------------------
# scatter plot
plot(alpha,log(f0),xlim = c(-15,0),main = 'Scatter plot of f0 vs. alpha')
rm(list=ls())

## -----------------------------------------------------------------------------
# data
u <- c(11,8,27,13,16,0,23,10,24,2)
v <- u+1
# iteration time
N <- 20
lambda <- numeric(N)
lambda[1] <- 1
# EM algorithm
for(i in 2:N){
  x <- (u*exp(-lambda[i-1]*u)-v*exp(-lambda[i-1]*v))/(exp(-lambda[i-1]*u)-exp(-lambda[i-1]*v))+1/lambda[i-1]
  lambda[i] <- length(u)/sum(x)
}
plot(lambda,ylab = 'lambda',type = 'l',main = 'EM algorithm solution results')
lambda[N]

score_function <- function(lambda){
  sum((u*exp(-lambda*u)-v*exp(-lambda*v))/(exp(-lambda*u)-exp(-lambda*v)))
}
# solve directly
solution <- uniroot(score_function,c(0,10))
solution$root
score_function(lambda[N])
solution$f.root

## -----------------------------------------------------------------------------
rm(list=ls())
a <- list(1,2,3,4,5)
unlist(a);typeof(unlist(a))
as.vector(a);typeof(as.vector(a))

## -----------------------------------------------------------------------------
1 == "1"
-1 < FALSE
"one" < 2

## -----------------------------------------------------------------------------
1 == as.numeric("1")
-1 < as.numeric(FALSE)
"one" < as.character(2)
rm(list=ls())

## -----------------------------------------------------------------------------
x <- c(1,2,3,4,5)
dim(x)

## -----------------------------------------------------------------------------
x <- matrix(c(1,2,3,4),2)
is.matrix(x)
is.array(x)
rm(list=ls())

## -----------------------------------------------------------------------------
x <- data.frame(X=c(1,2,3,7),Y=c('a','b','c',8))
attributes(x)

## -----------------------------------------------------------------------------
x
as.matrix(x)

## -----------------------------------------------------------------------------
y <- data.frame(X=c(),Y=c())
y
rm(list=ls())

## -----------------------------------------------------------------------------
scale01 <- function(x) {
rng <- range(x, na.rm = TRUE)
(x - rng[1]) / (rng[2] - rng[1])
}

## -----------------------------------------------------------------------------
set.seed(1)
# generate data
N <- 4*10
x <- rnorm(N);dim(x) <- c(10,4)
data <- data.frame(x)
data
cx <- sapply(data,scale01)
cx

## -----------------------------------------------------------------------------
#generate data
y <- rep(c('a','c'),c(3,7))
data2 <- data.frame(x,y=y)
data2
tx <- data.frame(lapply(data2,function(x){if (is.numeric(x)) scale01(x) else x}))
tx
rm(list=ls())

## -----------------------------------------------------------------------------
set.seed(1)
N <- 4*10
x <- rnorm(N);dim(x) <- c(10,4)
y <- rep(c('a','c'),c(3,7))
data <- data.frame(x) # numeric
data2 <- data.frame(x,y=y) # mix

## -----------------------------------------------------------------------------
data
cx <- vapply(data,sd,numeric(1))
cx

## -----------------------------------------------------------------------------
data2
tx <- vapply(data2[,vapply(data2, is.numeric, logical(1))],sd,numeric(1))
tx
rm(list=ls())

## -----------------------------------------------------------------------------
RandomR <- function(x0=0,y0=0,N=1e3){
  x <- y <- numeric(N)
  x[1] <- x0;y[1] <- y0
  for(i in 2:N){
    x[i] <- rnorm(1,0.9*y[i-1],sqrt(0.19))
    y[i] <- rnorm(1,0.9*x[i],sqrt(0.19))
  }
  return(list(x=x,y=y))
}

## -----------------------------------------------------------------------------
set.seed(1)
library(Rcpp)
sourceCpp('../src/RandomC.cpp')
N <- 1e3;x0 <- 0;y0 <- 0;
# Rcpp function
ResultC <- RandomC(x0,y0,N)
# R function
ResultR <- RandomR(x0,y0,N)

## -----------------------------------------------------------------------------
qqplot(ResultC[,1],ResultR$x,xlab='X(Rcpp)',ylab='X(R)',main='Q-Q plot for X')
qqplot(ResultC[,2],ResultR$y,xlab='Y(Rcpp)',ylab='Y(R)',main='Q-Q plot for Y')

## -----------------------------------------------------------------------------
library(microbenchmark)
ts <- microbenchmark(ResultC = RandomC(x0,y0,N),ResultR = RandomR(x0,y0,N))
summary(ts)[,c(1,3,5,6)]
rm(list=ls())

