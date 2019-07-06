##########Eberhart and Russell Model of stability


###Example 1
dat <- read.csv("E-R.csv",head=T)
names(dat) <- c("","Entry","Env","Rep","y")
head(dat)


####Example 2
#dat<- read.csv("dat.csv",head=T)
#E_R(dat$env,dat$rep,dat$geno,dat$y)


    


ige<- E_R(dat$Env,dat$Rep,dat$Entry,dat$y)
ige$Stability


E_R <- function(Env,Rep,Entry,y){

################# Begin

Rep <- as.factor(Rep)
Entry <- as.factor(Entry)
Loc <- as.factor(Env)
y <- as.numeric(y)

###########################
Entry_names<- levels(Entry)
r <- length(levels(Rep))
s <- length(levels(Loc))
t <- length(levels(Entry))

##########################
m1<- aov(y ~ Entry+ Loc + Loc: Entry)
anv1 <- anova(m1)
S2<- anv1$`Mean Sq`[4]


#################################################################
mydf = data.frame(aggregate (y ~ Entry + Loc,FUN=mean))
mydf$Loc <- as.factor(mydf$Loc)
mydf$Entry <- as.factor(mydf$Entry)
m2 <- aov(y~Entry+Loc+Loc:Entry,data = mydf)
anv2 <- anova(m2)
###########
Lxg <- tapply(y, list(Entry,Loc), mean)
Lmean <- apply(Lxg,2,mean)
Lsum <- apply(Lxg,2,sum)
Emean <- apply(Lxg,1,mean)
over.mean <- mean(Lxg)

#############################
I <- Lmean-over.mean
I2 <- sum(I^2)
YijI <- Lxg%*%I
bi <- YijI/I2 
Y2ij <-apply(Lxg^2, 1,sum) 
Y2 <- apply(Lxg, 1,sum)^2
Y2 <- Y2/s
VarS<- Y2ij-Y2
bi_I <- bi*YijI 
dij2 <-VarS -bi_I

####################
s2 <- anv1$`Mean Sq`[4]
sdi <- (dij2/(s-2))-(s2/r)
#####################################
SS.varieties <-anv2$`Sum Sq`[1]
SS.Env_varietiesxenv<- sum(Y2ij)-sum(apply(Lxg, 1,sum)^2)/s
SS.enviromen_linear <- anv2$`Sum Sq`[2]
SS.ExL_Lineal<- anv2$`Sum Sq`[2]
SS.varietiesxenviroment_linar <- sum(bi_I)-SS.ExL_Lineal
SS.pool.desviation <- sum(dij2)
ss.sd <- dij2 
SS.poolerr <-  anv1 $"Sum Sq"[4] / r
SS <-c(SS.varieties,SS.Env_varietiesxenv,SS.enviromen_linear,SS.varietiesxenviroment_linar
       ,SS.pool.desviation,ss.sd,SS.poolerr)
######################################
df_varieties<- t-1
df_Env_varietiesxenv<- t*(s-1)
df_enviromen_linear <- 1
df_varietiesxenviroment_linar <- t-1
df_pool.desviation <- t*(s-2)
df_sd <- rep(s-2,length(dij2))
df_pooerr <- s*t*(r-1)
df <- c(df_varieties,df_Env_varietiesxenv,df_enviromen_linear,df_varietiesxenviroment_linar
        ,df_pool.desviation,df_sd,df_pooerr)
##################################


ms_varieties<- SS.varieties /df_varieties
ms_Env_varietiesxenv<-SS.Env_varietiesxenv/df_Env_varietiesxenv
ms_enviromen_linear <- SS.enviromen_linear/df_enviromen_linear
ms_varietiesxenviroment_linar <- SS.varietiesxenviroment_linar/df_varietiesxenviroment_linar
ms_pool.desviation <- SS.pool.desviation/df_pool.desviation 
ms_sd <-ss.sd/df_sd
ms_pooerr<- SS.poolerr/df_pooerr
  
ms <- c(ms_varieties,ms_Env_varietiesxenv,ms_enviromen_linear, ms_varietiesxenviroment_linar,
        ms_pool.desviation,ms_sd,ms_pooerr )

#################fc
fvarieties <- ms_varieties/ms_varietiesxenviroment_linar 
fvarietiesxenviroment_linar <- ms_varietiesxenviroment_linar /ms_pool.desviation

f <-ms_sd/ ms_pooerr

ft <- c(fvarieties,NA,NA,fvarietiesxenviroment_linar,NA,f,NA)

PLINES <-  1- pf(f,df[7], df_pooerr)
pval = c(1- pf(fvarieties, df_varieties, df_pool.desviation ), NA, NA,1- pf(fvarietiesxenviroment_linar, df[4], df_pool.desviation),NA,PLINES, NA)


anovadf <- data.frame(df, `Sum Sq`=SS, `Mean Sq`=ms, `F value`= ft,`Pr(>F)`= pval,check.names=FALSE)
rownames(anovadf) <- c("Genotypes","Env + (Gen x Env)", "Env (linear)", " Gen x Env(linear)", "Pooled deviation",levels(Entry), "Pooled error" ) 
class(anovadf) <- c("anova","data.frame")

###########stability

stb <- data.frame("Genotype"=Entry_names,"Means"=Emean,"bi"=bi,"Sdi"=sdi)
stb

cat(" Stable genotype have bij = 1 and sdij = 0", "\n\n"  )


plot(1, type = 'n', xlim = range(c(stb[,3])) ,ylim = range(c(stb[,3], stb[,4])), xlab = "bi",
     ylab ="Sdi")

text(stb[,3], stb[,4], labels = Entry_names, adj = c(0.5, 0.5), col = "blue")
abline(h = 0, v=1,lty = 1)

result <-list(ANOVA=anv1,IGE=anovadf,Stability=stb)
return(result)
}
