####Tai stability analysis
###Tai(Entry,Enviroment,Rep,y,"Name our variable")


#Example 1
dat <- read.csv("F-W.csv",head=T)
head(dat)
Tai(dat$Entry,dat$Env,dat$Rep,dat$y,"y")
##Example 2

##Tai(dat$entry,dat$env,dat$rep,dat$y,"Y")
##dat <- read.csv("E-R.csv",head=T)
##Tai(dat$Entry,dat$Env,dat$Rep,dat$y,"Grain kl ha")

Tai <- function( geno, env, rep,y,trait) {
  
  # Everything as factor
  
  geno <- as.factor(geno)
  env <- as.factor(env)
  rep <- as.factor(rep)
  y <- as.numeric(y)
  
  # Check data

  
  geno.num <- nlevels(geno)
  env.num <- nlevels(env)
  rep.num <- nlevels(rep)
  

  # Compute interaction effects matrix
  
  int.mean <- tapply(y, list(geno,env), mean, na.rm = T)
  
  overall.mean <- mean(int.mean)
  env.mean <- apply(int.mean, 2, mean)
  geno.mean <- apply(int.mean, 1, mean)
  int.eff <- int.mean + overall.mean
  for (i in 1:env.num) int.eff[,i] <- int.eff[,i] - geno.mean
  for (i in 1:geno.num) int.eff[i,] <- int.eff[i,] - env.mean
  
  # ANOVA
  
  model <- aov(y ~ geno + env +rep %in% env + geno:env)  
  at <- summary(model)
  
  
  # Compute Tai values alpha and lambda
  
  slgl <- int.eff
  for (i in 1:geno.num) slgl[i,] <- slgl[i,]*(env.mean - overall.mean)/(env.num-1)
  alpha <- apply(slgl, 1, sum)/(at[[1]][2,3]-at[[1]][3,3])*geno.num*rep.num
  
  s2gl <- int.eff
  for (i in 1:geno.num) s2gl[i,] <- s2gl[i,]^2/(env.num-1)
  lambda <- (apply(s2gl, 1, sum) - alpha*apply(slgl, 1, sum))/(geno.num-1)/
    at[[1]][5,3]*geno.num*rep.num
  
  # plot lambda limits
  conf <- 0.95
  
  lmax <- max(c(lambda, qf(1-(1-conf)/2, env.num-2, env.num*(geno.num-1)*(rep.num-1))))*1.1
  
  # Prediction interval for alpha
  
  lx <- seq(0,lmax,lmax/100)
  ta <- qt(1-(1-conf)/2, env.num-2)
  pi.alpha <- ta * ((lx * (geno.num-1) * at[[1]][5,3] * at[[1]][2,3]) /
                      (at[[1]][2,3] - at[[1]][3,3]) /
                      ((env.num-2) * at[[1]][2,3] - (ta^2 + env.num - 2)
                       * at[[1]][3,3]))^.5
  
  # plot alpha limits
  
  amax <- max(c(abs(alpha), pi.alpha))
  
  # Tai plot
  
 
  
  plot(1, type = "n", xlim = c(-0.05*lmax, lmax), ylim = c(-amax, amax),
       main = paste("Tai stability analysis for ", trait, sep = ""), xlab = expression(lambda), ylab = expression(alpha))
  points(lambda, alpha, col = "black", lwd = 2, pch = 16)
  text(lambda, alpha, labels = names(alpha), col = "blue", pos = 1, offset = .3)
  points(lx, pi.alpha, type = "l", lty = 1, col = "black")
  points(lx, -pi.alpha, type = "l", lty = 1, col = "black")
  abline(v = qf((1-conf)/2, env.num-2, env.num*(geno.num)*(rep.num-1)),
         lty = 1, col = "red")
  abline(v = qf(1-(1-conf)/2, env.num-2, env.num*(geno.num)*(rep.num-1)),
         lty = 1, col = "red")
  
  # Output
  
  coords <- cbind(alpha, lambda)
  rs <- as.data.frame(coords)
  res_frame<-cbind(rs,geno.mean) 
  names(res_frame) <- c("Alpha", "lamba", "Genotype Mean")
  result <- list(Anova=at,Stability=res_frame)
  return(result)
  
}
