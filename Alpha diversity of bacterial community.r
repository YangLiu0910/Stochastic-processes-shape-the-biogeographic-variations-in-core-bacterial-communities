###ANOVA and  mutiple comparisons for  alpha diversity 
library(MASS)
#ANOVA
fit <- aov(Shannon ~ Compartment*Sites,data=alpha_diversity)
summary(fit)
# Multiple comparisons
library(multcomp)
tuk <- glht(fit, linfct = mcp(group = "Tukey"))
plot(cld(tuk, level = 0.05), col = "lightgrey")

#Boxplots for alpha diveristy
library(ggplot2)
library(magrittr)
library(ggpubr)
df=alpha_diversity
df$Compartment<-factor(df$Compartment,levels=c("BS","RS","R","S","L"))
df$Sites<-factor(df$Sites,levels=c("d","n","k","h","s"))
df$Split<-factor(df$Split,levels=c("Total","Core"))

p <- ggboxplot(df, x = "Compartment", y = "ACE", fill = "Compartment", 
       facet.by="Split",size=0.1,alpha=0.8)
p









