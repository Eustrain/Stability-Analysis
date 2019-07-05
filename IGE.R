######## MODEL AMMI AND SREG 

###method=AMMI
###method=2SREG
###model=1 RCBD Design
####model=2 Lattice
####biplot= 0 and 1 plot x=yield, y=PC1 (1 style differente that 0)
####biplot=2 and 3  plot x=CP1 and y=CP2 (2 style different that 3)

#IGE(Block,Enviroment,Repetition,Entry,y,method="AMMI",model=1 or 2,biplot= 0, 1,2, or 3)

#we have not blocks :
##IGE(NULL,Enviroment,Repetition,Entry,y,method="SREG",model=1 or 2,biplot= 0, 1,2, or 3)


library(dplyr)
library(agricolae)

fm <- read.csv("Dat.csv",head=T)

###Example 1 without Blocks (BRCBD) model=1
data("plrv")


Ex1<- IGE(NULL,plrv$Locality,plrv$Rep, plrv$Genotype,plrv$Yield,method = "AMMI",model = 2,biplot = 2)


###Example 2 with Blocks (Lattice), model=2  and SREG

Ex2<- IGE(Block = fm$BLOCK,fm$env,fm$rep,fm$entry,fm$y,method="SREG",model = 2,biplot =1)

tab1 <- cbind(Ex2$Means_Entry,Ex2$PC_Entry[,2])
cor(tab1$MEDIA,tab1$`Ex2$PC_Entry[, 2]`)

####The correlation between CP1 and Yield is high, then we can
### interpret PC1 like yield and PC2 like the effect of IGE
### a genotype with high value in CP1 and low value in CP2 is ideal

IGE <- function(Block=NULL,Env,Rep,Entry,y,method = c("AMMI", "SREG"),model=c(1,2),biplot = c(0,1,2,3)){

  
  
if (!is.null(Block)){
Block <- as.factor(Block)
}

Loc <- as.factor(Env)
Rep <- as.factor(Rep)
Entry <- as.factor(Entry)
y <- as.numeric(y)
######################################################
Entry_names<- levels(Entry)
Loc_names <- levels(Loc)
Entry.num <- length(Entry_names)
Loc.num <- length(Loc_names)
nrep <- length(levels(Rep))
#######################################################
######### Statistics#################################
overall.mean <- mean(y)
Loc.mean <- tapply(y, Loc, mean)
Entry.mean <- tapply(y, Entry, mean)
########################################################
if(model==1){
m1 <- aov(y~Entry+Rep:Loc+Entry:Loc+Loc)##################RCBD Design
anv1 <- anova(m1)
}
else if( model==2){
m1<- aov(y~Entry+Rep:Loc+Block:Rep:Loc+Entry:Loc+Loc)######Lattice
anv1 <- anova(m1)
}
##################AMMI
if(method=="AMMI"){
geno_env <- model.tables(m1, type = "effects", cterms = "Entry:Loc")
z<- geno_env$tables$`Entry:Loc`
}
##################SREG

else if( method =="SREG"){
A <- tapply(y, list(Entry, Loc), mean, na.rm = TRUE)
env.mean <- apply(A, 2, mean)
z <- t(t(A) - env.mean)
}
############################################################

############Descomposition of values singulars

DVS <- svd(z)
D <- diag(DVS$d)
U <- DVS$u %*%sqrt(D)
V <- DVS$v %*%sqrt(D)


############## PCA##############

PC <- min(Loc.num, Entry.num) - 1
PC.num <- paste0("PC", c(1:PC))
PC.dsv <- DVS$d[1:PC]^2
PC.count <- PC.dsv / sum(PC.dsv) * 100
PC.acum <- cumsum(PC.count)
Ecolnumb <- c(1:PC)
Ecolnames <- paste("PC", Ecolnumb, sep = "")
###############ANOVA PC`S##########
residual.SS <- deviance(m1)
residual.DF <- m1$df.residual
residual.MS <- residual.SS/residual.DF 

####################Sum Square PC`S##########
int.SS <- (t(as.vector(z)) %*% as.vector(z))*nrep

PC.SS <- DVS$d[1:PC]^2*nrep
PC.SS <- round(PC.SS,4)
PC.DF <- Entry.num + Loc.num - 1 - 2*Ecolnumb
PC.residual.SS <- round(int.SS - sum(PC.SS),5)
PC.residual.DF <- ((Entry.num - 1)*(Loc.num - 1)) - sum(PC.DF)

PC.SS[PC + 1] <- PC.residual.SS
PC.DF[PC + 1] <- PC.residual.DF
#########################Means square####

MS <- PC.SS/PC.DF
MS <- round(MS,4)
F <- MS/residual.MS
F <- round(F,2)
probab <- pf(F, PC.DF, residual.DF, lower.tail = FALSE)
probab <- round(probab,4)

#####################Table anova################
percSS <- PC.SS/int.SS
rowlab <- c(Ecolnames, "Residuals")
######################Result ANova###################
anv_gl <- c(anv1$Df[1],anv1$Df[2],anv1$Df[4])
anv_ss <- c(anv1$`Sum Sq`[1],anv1$`Sum Sq`[2],anv1$`Sum Sq`[4] )
anv_ms <- c(anv1$`Mean Sq`[1],anv1$`Mean Sq`[2],iga_ms <- anv1$`Mean Sq`[4])
anv_f <- c(anv1$`Pr(>F)`[1],anv1$`Pr(>F)`[2],anv1$`Pr(>F)`[4])
effec <-c("GEN","ENV","ENV*GEN") 
anv.count <- anv_ss / sum(anv_ss) * 100
anv.acom <- cumsum(anv.count)
ftab <- c(anv1$`F value`[1],anv1$`F value`[2],anv1$`F value`[4])
tab_anv <- data.frame(Effect=effec,SS=anv_ss,PORCENT=anv.count,ACOMULATED=anv.acom
                      ,DF=anv_gl,MS=anv_ms,F=ftab,Prob=anv_f)



#######################Result PC`S####################
PC.count <- round(PC.count,5)
PC.acum <- round(PC.acum,5)
Porcent_pc<- c(PC.count,"Residuals")
Acum_pc <- c(PC.acum,"Residuals")
PC_aov<- data.frame(Effect = rowlab, SS = PC.SS,PORCENT=Porcent_pc,ACOMULATED=Acum_pc,
                         DF = PC.DF, MS = MS, F = F, Pr.F = probab)



###############Plot AMMI model AMMI1
if(biplot== 0){
  
  plot(1, type = 'n', xlim = range(c(Loc.mean, Entry.mean)), ylim = range(c(V[,1], U[,1])), xlab = "Yield",
       ylab =paste("PC1",round(PC.count[1],2),"%",sep=" "))
  
  text(Entry.mean, U[,1], labels = Entry_names, adj = c(0.5, 0.5), col = "blue")
  text(Loc.mean , V[,1], labels = Loc_names, adj = c(0.5, 0.5), col = "red")
  abline(h = 0, v = overall.mean, lty = 1)
  for (i in 1:Loc.num){  #    enviroment
    lines(c(overall.mean,Loc.mean [i]), c(0, V[i, 1]),col = "black", lty = 1,lwd=2) }
}

else if(biplot ==1){
  plot(1, type = 'n', xlim = range(c(Loc.mean , Entry.mean)), ylim = range(c(V[,1], U[,1])), xlab = "Yield",
       ylab =paste("PC1",round(PC.count[1],2),"%",sep=" "))
  text(Entry.mean, U[,1], labels = Entry_names, adj = c(0.5, 0.5), col = "blue",pos = 3)
  text(Loc.mean, V[,1], labels = Loc_names, adj = c(0.5, 0.5), col = "red",pos = 3)
  points(Loc.mean, V[,1], col = "red", lwd = 2, pch = 17)
  points(Entry.mean, U[,1], col = "blue", lwd = 2, pch = 15)
  abline(h = 0, v = overall.mean, lty = 1)
  for (i in 1:Loc.num){  #    enviroment
    lines(c(overall.mean,Loc.mean[i]), c(0, V[i, 1]),col = "black", lty = 1,lwd=2) }
}


####################### Plot Model AMMI2
else if(biplot == 2) {

plot(1, type = 'n', asp = 1, xlim = range(c(U[, 1], V[, 1])), ylim = range(c(U[, 2], V[, 2])),
     xlab = paste("PC1",round(PC.count[1],2),"%",sep=" "), ylab = paste("PC2",round(PC.count[2],2),"%",sep =" "))
text(U[, 1], U[, 2], labels = Entry_names, adj = c(0.5, 0.5), col = "blue", pos = 3)

text(V[, 1], V[, 2], labels = Loc_names, adj = c(0.5, 0.5), col = "red", pos = 3)
abline(h = 0, v = 0, col = "black", lty = 1)
#####rows
for (i in 1:Entry.num) {  ###genotype
  arrows(0,0,c(0, U[i, 1]), c(0, U[i, 2]), col = "green", lty = 1,lwd=2)}
                                                                   
for (i in 1:Loc.num) {  ###enviroment
  arrows(0,0,c(0, V[i, 1]), c(0, V[i, 2]), col = "green", lty = 1,lwd=2) }
                                                                       
} 

else if(biplot==3){
  
  plot(1, type = 'n', asp = 1, xlim = range(c(U[, 1], V[, 1])), ylim = range(c(U[, 2], V[, 2])),
xlab = paste("PC1",round(PC.count[1],2),"%",sep=" "), ylab = paste("PC2",round(PC.count[2],2),"%",sep =" "))
points(U[, 1], U[, 2], col = "blue", lwd = 2, pch = 15)
text(U[, 1], U[, 2], labels = Entry_names, adj = c(0.5, 0.5), col = "blue", pos = 3)
points(V[, 1], V[, 2], col = "red", lwd = 2, pch = 17)
text(V[, 1], V[, 2], labels = Loc_names, adj = c(0.5, 0.5), col = "red", pos = 3)
abline(h = 0, v = 0, col = "black", lty = 1)
  

  
############# aggregate lines
for (i in 1:Entry.num) {  #  genotype
  lines(c(0, U[i,1]), c(0, U[i, 2]), col = "black", lty = 1,lwd=2)  }
                                                       

for (i in 1:Loc.num){  #    enviroment
  lines(c(0, V[i,1]), c(0, V[i, 2]),col = "black", lty = 1,lwd=2) }
}


###################Statistics

Means<- as.data.frame(Entry.mean)
Loc_means <- as.data.frame(Loc.mean)
##########ASV

if( method=="AMMI"){
P.C_GEN <- round(U[,1:2],5)
P <-PC.SS[1]/PC.SS[2]
ASV <- round(sqrt(P*(P.C_GEN[,1])^2+(P.C_GEN[,2])^2),4)
rk_asv<-rank(ASV)
tab_out <-data.frame("Entry"=Entry_names,"MEDIA"=Means$Entry.mean,
                     "ASV"=ASV,"Rank ASV"=rk_asv)

}


else if(method=="SREG"){

tab_out <-data.frame("Entry"=Entry_names,"MEDIA"=Means$Entry.mean)
}

Entry_PC <- U[,1:2]
Loc_PC <- V[,1:2]

PC_Entry <- data.frame("Entry"=Entry_names,"CP1"=Entry_PC[,1],"PC2"=Entry_PC[,2])
PC_ENV <- data.frame("Env"=Loc_names,"CP1"=Loc_PC[,1],"PC2"=Loc_PC[,2])
##########################
    
result <-list(ANOVA=anv1,ANOVA_AMMI=tab_anv,ANOVA_PC=PC_aov,Means_Entry=tab_out,
              Means_Loc=Loc_means,PC_Entry=PC_Entry,PC_ENV=PC_ENV)
return(result) 
}



