##############################################
#
# I did this when I was young and needed the money ;)
# Seriously, although I think the results are valid,
# I would run and report my GLMMs differently as of
# today. - Roland, 2015-11-09
#
##############################################

rm(list=ls(all=TRUE))

library(fmsb)
library(boot)
library(DAAG)
library(lmtest)

##############################################
#   ---- Data input/reordering/transf. ----

# Load data.
arti <- read.table("allQsMay2013.csv", header=TRUE, sep="\t", quote="\"")

# Create / Reorder factors.

arti$Kurz <- factor(arti$Kurz, levels=c("0", "1"))
arti$Gen <- factor(arti$Gen, levels=c("Fem", "Mask", "Neut"))
arti$Kas <- factor(arti$Kas, levels=c("Akk", "Nom", "Dat"))
arti$Adj <- factor(arti$Adj, levels=c("0", "1"))
arti$Linksk <- factor(arti$Linksk, levels=c("And", "Praep", "Int", "So"))
arti$Linksv <- factor(arti$Linksv, levels=c("0", "1"))
arti$Rechtsv <- factor(arti$Rechtsv, levels=c("0", "1"))

##############################################
#   ------- Calc and report model  --------

arti.glm <- glm(Kurz~+Gen*Kas+Adj+Linksk, data=arti, family=binomial(link=logit))
arti.glm.sum <- summary(arti.glm)
print(arti.glm.sum)

cat("\n --- Odds ratios ---\n")
print(exp(coef(arti.glm)))
print(exp(confint(arti.glm)))

cat("\n --- GenKas-Coeffs ---\n")
cat("Mask Nom: ", as.numeric(arti.glm$coefficients["(Intercept)"]+arti.glm$coefficients["GenMask"]+arti.glm$coefficients["KasNom"]+arti.glm$coefficients["GenMask:KasNom"]),"\n")
cat("Mask Akk: ", as.numeric(arti.glm$coefficients["(Intercept)"]+arti.glm$coefficients["GenMask"]),"\n")
cat("Mask Dat: ", as.numeric(arti.glm$coefficients["(Intercept)"]+arti.glm$coefficients["GenMask"]+arti.glm$coefficients["KasDat"]+arti.glm$coefficients["GenMask:KasDat"]),"\n")
cat("Neut Nom: ", as.numeric(arti.glm$coefficients["(Intercept)"]+arti.glm$coefficients["GenNeut"]+arti.glm$coefficients["KasNom"]+arti.glm$coefficients["GenNeut:KasNom"]),"\n")
cat("Neut Akk: ", as.numeric(arti.glm$coefficients["(Intercept)"]+arti.glm$coefficients["GenNeut"]),"\n")
cat("Neut Dat: ", as.numeric(arti.glm$coefficients["(Intercept)"]+arti.glm$coefficients["GenNeut"]+arti.glm$coefficients["KasDat"]+arti.glm$coefficients["GenNeut:KasDat"]),"\n")
cat("Fem Nom: ", as.numeric(arti.glm$coefficients["(Intercept)"]+arti.glm$coefficients["KasNom"]),"\n")
cat("Fem Akk: ", as.numeric(arti.glm$coefficients["(Intercept)"]),"\n")
cat("Fem Dat: ", as.numeric(arti.glm$coefficients["(Intercept)"]+arti.glm$coefficients["KasDat"]),"\n")

cat("\n --- Nagelkerke R2 ---\n")
arti.r2 <- NagelkerkeR2(arti.glm)
print(arti.r2$R2)

#cat("\n --- Likelihood ratio test ---\n")
#print(lrtest(arti.glm0, arti.glm))

print(waldtest(arti.glm))

cat("\n --- 10-fold cross-validation ---\n")
cv.binary(arti.glm)
#cost <- function(r,pi=0) mean(abs(r-pi)>0.5)
arti.cv <- cv.glm(arti, arti.glm, K = 10)
print(arti.cv$delta)

# Plot coeffs
arti.glm.coefs <- arti.glm.sum$coefficients
arti.glm.coefs.intercept <- arti.glm.coefs["(Intercept)",]
arti.glm.coefs.nintercept <- arti.glm.coefs[which(rownames(arti.glm.coefs)!="(Intercept)"),]
sortorder <- sort.list(arti.glm.coefs.nintercept[,"Estimate"])
sigind <- as.numeric(which(arti.glm.coefs.nintercept[,"Pr(>|z|)"]<0.1))
arti.glm.finalcoefs <- arti.glm.coefs[c(1,intersect(sortorder, sigind)+1),]
cplot.vals <- as.numeric(arti.glm.finalcoefs[,"Estimate"]) # use exp to make ORs
cplot.ger.vals <- sub("\\.",",",as.character(round(cplot.vals,3)),perl=T)
cplot.names <- as.character(rownames(arti.glm.finalcoefs))
cplot.density <- c(20,rep(0,length(cplot.vals-1)))
cplot.ypos <- ifelse(cplot.vals<=0,cplot.vals-0.2,cplot.vals+0.1)
cplot.ylim <- c(min(cplot.vals)*1.6,max(cplot.vals)*1.4)
cplot.offsets <- rep(c(0.2,0),length(cplot.names)/2)
cplot.offsets <- cplot.offsets[1:length(cplot.names)]
cplot.nameypos <- cplot.offsets+(cplot.ylim[1]+0.1)


svg("coefplot.svg")
cplot.plot <- barplot(cplot.vals, names.arg="", density=cplot.density, ylim=cplot.ylim)
text(cplot.plot, cplot.ypos, labels=cplot.ger.vals, cex=0.9, family="serif")
text(cplot.plot, cplot.nameypos, cplot.names, cex=0.7, family="serif")
dev.off()