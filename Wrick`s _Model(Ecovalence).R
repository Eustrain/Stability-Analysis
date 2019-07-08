
###Wricke's Model: Ecovalence
##As percentage of ecovalence (Wi) is 
##inversely associated with phenotypic stability, 
##a low percentage of Wi would  indicate high stability of performance.


dat <- read.csv("dat.csv",head=T)
head(dat)
dat <- dplyr::select(dat,Entry,Loc,Rep,y)

Wi(dat$geno,dat$env,dat$rep,dat$y)


Wi <- function(Entry,Loc,Rep,y){

Entry <- as.factor(Entry)
Loc <- as.factor(Loc)
re <- as.factor(Rep)
y <- as.numeric(y)

ixg <-tapply(y,list(Entry,Loc),mean) 
overall.mean <- mean(ixg)
env.mean <- apply(ixg, 2, mean)
geno.mean <- apply(ixg, 1, mean)

int.eff <- ixg+overall.mean-geno.mean
int.eff <- t(t(int.eff) - env.mean)
s <- apply(int.eff , 1,function(x){
  x^2
})

WI<- apply(s, 2, sum)
WI_POR<- WI/sum(WI)*100

print(cat(" low percentage of Wi  indicate high stability of performance" ))
df <- data.frame("Genotype"=geno.mean,"Wi"=WI,"Wi(%)"=WI_POR,check.names = FALSE)
return(df)
}



