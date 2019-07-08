# Regression Analysis of Stability 
# (Finlay and Wilkinson's Model)
###
#Finlay and Wilkinson defined that b = 1.0 indicates
#average ','stability" , b < 1.0 means " below average stability" and b > 
#1,0 suggests "above average stability". 
#fw(Entry,Enviroment,Rep,y,"Name our variable")

dat <- read.csv("F-W.csv",header = T)

dat <- dplyr::select(dat,Entry,Loc,Rep,y)
names(dat) <- c("Entry","Env","Rep","y")
head(dat)



f<- fw(dat$Entry,dat$Env,dat$Rep,dat$y,"Yield")

f

fw<- function(Entry,Env,Rep,y,trait){


Entry <-as.factor (Entry)
Loc <- as.factor(Env)
rep <- as.factor(Rep)
y <- as.numeric(y)
geno.num <- nlevels(Entry)
env.num <- nlevels(Loc)
rep.num <- nlevels(rep)
Entry_names<- levels(Entry)


#####################anova

m1 <- aov(y~Entry+Loc+rep:Loc+Entry:Loc)
an1 <- anova(m1)


################# Regression anovo

int.mean <- tapply(y,list(Entry,Loc), mean)

overall.mean <- mean(int.mean, na.rm=T)
env.mean <- apply(int.mean, 2, mean, na.rm=T)
geno.mean <- apply(int.mean, 1, mean, na.rm=T)


rgss<- apply(int.mean, 2, function(x){
  x/mean(x)
})


RGss<- (apply(rgss, 1, sum)^2)
RGss <- (rep.num)^2*(sum(RGss))
Rgms <- RGss/(geno.num-1)
Drss <- an1$`Sum Sq`[4]-RGss
Drms <- Drss/((geno.num-1)*(env.num-1)-(rep.num-1))
SS <- c(RGss,Drss )
MS <- c(Rgms,Drms)
f <- MS/SS
d_f <- c((geno.num-1),((geno.num-1)*(env.num-1)-(rep.num-1)))
pr <-  c(pf(f[1],d_f[1], an1$Df[5]),pf(f[2],d_f[2], an1$Df[5]))

DF <- data.frame(`Df`=d_f,`Sum sq`=SS,`Ms sq`=MS,F=f,`Pr(>F)`= pr,check.names=FALSE)
rownames(DF) <- c("Regression","Dev.Regr") 
class(DF) <- c("anova","data.frame")


### Regression-stability for environments


b <- NULL
se <- NULL

for (i in 1:geno.num){
  modelo <- lm(int.mean[i,] ~ env.mean)
  b[i] <- coef(modelo)[2]
  se[i] <- summary.lm(modelo)$coefficients[2,2]
}

stability <- cbind(geno.mean ,b, se)
stability <- as.data.frame(stability )
names(stability) <- c("Genotype","bi","SE")


plot(1, type = 'n', xlim = range(c(stability[,1])) ,ylim = range(stability[,2]), xlab = trait,
     ylab ="bi",main=paste("Finlay and Wilkinson's Model for ", trait, sep = ""))

text(stability[,1], stability[,2], labels = Entry_names, adj = c(0.5, 0.5), col = "blue",pos = 3)
points(stability[,1], stability[,2], col = "black", lwd = 2, pch =16)
abline(h = 1, v=overall.mean ,lty = 1)

result <- list(ANOVA=an1,FIW=DF,stability=stability)

return(result)

}
