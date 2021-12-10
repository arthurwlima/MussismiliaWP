##################### Hypergometric test ###############################
#### Comparar a ocorrencia de cada familia viral nos subsets ###########
#### de transcritos exclusivos (H/WP), comparado com o eperado #########
#### pela distribuicao global. Corrigido por Bonferroni???     #########

dX <- read.table("../VirusFamily_Trat.csv", header=T, sep="\t")

dX$Fam <- as.character(dX$Fam)
dHypFam <- data.frame(cbind(dX$Fam, phyper(dX$H_exc, dX$TOTAL, sum(dX$TOTAL)-dX$TOTAL, sum(dX$H_exc)), phyper(dX$WP_exc, dX$TOTAL, sum(dX$TOTAL)-dX$TOTAL, sum(dX$WP_exc))))


> dHypFam[ dHypFam$Healthy_p > 1- 0.05,]
             Family Healthy_p      WP_p
1   Phycodnaviridae 0.9592761 0.6432431
9  Marseilleviridae 0.9707572 0.7703248
13    Herpesviridae 0.9998060 0.2251565
> dHypFam[ dHypFam$Healthy_p < 0.05,]
      Family  Healthy_p      WP_p
7 Myoviridae 0.01057043 0.9995494
> 
> 
> dHypFam[ dHypFam$WP_p > 1- 0.05,]
         Family  Healthy_p      WP_p
7    Myoviridae 0.01057043 0.9995494
11 Siphoviridae 0.62011628 0.9997382
> dHypFam[ dHypFam$WP_p < 0.05,]
        Family  Healthy_p       WP_p
2  Mimiviridae 0.09735253 0.01892303
6 Retroviridae 0.70335544 0.02351148


#################### Bonferroni Correction ###########################
> 0.05/15
[1] 0.003333333


dHypFam[ dHypFam$Healthy_p > 1- 0.05/15,]
          Family Healthy_p      WP_p
13 Herpesviridae  0.999806 0.2251565

> dHypFam[ dHypFam$Healthy_p <  0.05/15,]
[1] Family    Healthy_p WP_p     
<0 rows> (or 0-length row.names)


dHypFam[ dHypFam$WP_p > 1- 0.05/15,]
         Family  Healthy_p      WP_p
7    Myoviridae 0.01057043 0.9995494
11 Siphoviridae 0.62011628 0.9997382

> dHypFam[ dHypFam$WP_p < 0.05/15,]
[1] Family    Healthy_p WP_p     
<0 rows> (or 0-length row.names)



#########################  Barplot ####################################


> rownames(dRev)[12] <- "unc. dsDNA viruses"
> dRev$dumb <- 0

png("ViralFamilies.png", 500, 600)
par(mar=c(5,10,2,2), cex=1)
barplot(as.matrix(t(dRev[,c(11, 12, 13, 5)])), horiz=T, beside=T, col=c('lightblue', 'orangered4', 'white', 'gray'), las=1, width=c(0.9, 0.9, 0.4, 1.1), space=c(0,1.2), xlim=c(0,.3))
legend('bottomright', c('Total', 'WP', 'H'), pch=15, col=c('gray', 'orangered4', 'lightblue'), inset=0.01, cex=1.1)
### marcando as familias over/under representadas, Bonferroni
text(t(dRev[c(3),11])+0.01, mp[1,c(3)] -0.5, labels= "*", cex=1.3)
text(t(dRev[c(9,5),12])+0.01, mp[2,c(9,5)] -0.5, labels= "*", cex=1.3)
dev.off()





############# HIGHER TAXONOMY #########################################
dFam <- data.frame( tapply(dTax$MBH1, dTax$V20, sum), tapply(dTax$MBH3, dTax$V20, sum), tapply(dTax$MBH4, dTax$V20, sum),
tapply(dTax$WP02, dTax$V20, sum), tapply(dTax$WP04, dTax$V20, sum), tapply(dTax$WP05, dTax$V20, sum),
tapply(dTax$WP06, dTax$V20, sum), tapply(dTax$WP07, dTax$V20, sum), tapply(dTax$WP09, dTax$V20, sum))

data.frame( tapply(dTax$MBH1, dTax$V21, sum), tapply(dTax$MBH3, dTax$V21, sum), tapply(dTax$MBH4, dTax$V21, sum),
tapply(dTax$WP02, dTax$V21, sum), tapply(dTax$WP04, dTax$V21, sum), tapply(dTax$WP05, dTax$V21, sum),
tapply(dTax$WP06, dTax$V21, sum), tapply(dTax$WP07, dTax$V21, sum), tapply(dTax$WP09, dTax$V21, sum))

