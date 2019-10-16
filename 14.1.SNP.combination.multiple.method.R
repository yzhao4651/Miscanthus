##import the data file with total SNPs detected 
filename <- read.csv("Result.8.3/flo.8.3.csv")
str(filename)
###change all of traits to the four temporal scales 
levels(filename$Trait.name)[levels(filename$Trait.name)=="HD_1"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HW_1"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="GHW_1"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HM_1"] <- "M"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HD_50"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HW_50"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="GHW_50"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="HM_50"] <- "M"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FD_1"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FW_1"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="GFW_1"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FM_1"] <- "M"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FD_50"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FW_50"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="GFW_50"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="FM_50"] <- "M"
levels(filename$Trait.name)[levels(filename$Trait.name)=="fprind"] <- "D"
levels(filename$Trait.name)[levels(filename$Trait.name)=="fprinW"] <- "W"
levels(filename$Trait.name)[levels(filename$Trait.name)=="fprinGW"] <- "GW"
levels(filename$Trait.name)[levels(filename$Trait.name)=="fprinM"] <- "M"
library(tidyr)
library(dplyr)
##select the column need for study
t <- as.data.frame(ftable(filename[,c(15,18,19)]))
t.T.1 <- droplevels(t)
str(t.T.1)
###unit the SNP. method and trait.nmes 
t.T.1.u <-  unite_(t.T.1, "IndMT", c("Ind.SNP","Method", "Trait.name"),sep=";")
##change the column name Freq to  "SNPs.Total.detected"
colnames(t.T.1.u)[colnames(t.T.1.u)=="Freq"] <- "SNPs.Total.detected"
###import the SNP hit with flowering genes 
hitflowering <- read.csv("flo.8.4/hitflowering.csv")
##select the SNPs hit with flowering but in the totole detected SNPs. 
hitflo <- filename[filename$SNP %in% hitflowering$SNP,]
hitflo <- droplevels(hitflo)
tflo <- as.data.frame(ftable(hitflo[,c(15,18,19)]))
nlevels(hitflo$SNP)
###unite the Method and Trait.name
###unit the SNP. method and trait.nmes in order to use to merge with the total detected SNPs
t.T.2.u <-  unite_(tflo, "IndMT", c("Ind.SNP", "Method", "Trait.name"),sep=";")
## change the column name Freq to "SNPs.hit.with.flowering.genes"
colnames(t.T.2.u)[colnames(t.T.2.u)=="Freq"] <- "SNPs.hit.with.flowering.genes"
##merge two files t.T.1.u, and t.T.2.u
#t.1 <- merge(x=t.T.1.u, y=t.T.2.u, by ="IndMT", by.x=TRUE)## this one does not work 
t.2 <- plyr::join_all(list(t.T.1.u,t.T.2.u), by="IndMT")
t.2[is.na(t.2)] <- 0
##separate MT to method and Trait.name after merge

t.13 <- separate(t.2, "IndMT", c("Ind.SNP","Method", "Trait.name"),sep=";")

##transform from wide form to vertical form. 
library(reshape2)
library(reshape)
t.14 <- melt(t.13, id.vars = c("Ind.SNP","Method","Trait.name"))
str(t.14)
###F116+36088im
levels(as.factor(t.14$Method))
t.15 <- t.14[t.14$Method=="MLM+SUPER" | t.14$Method=="CMLM+SUPER",]
t.16 <- t.15[t.15$Ind.SNP=="F116+36088im",]
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.detected"))
t.16$Ind.SNP <- factor(t.16$Ind.SNP,levels = c("F116+36088im", "F106+4185", "F116+3077", "F122+2272"))
t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
MLMSUPER36088im <- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) + facet_grid(~Method+Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=0), axis.text.y=element_text(size=12))+
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="F116+36088im")

###106+4185, F116+3077, and F122+2772
t.14 <- melt(t.13, id.vars = c("Ind.SNP","Method","Trait.name"))
levels(as.factor(t.14$Method))
t.15 <- t.14[t.14$Method=="MLM+SUPER" | t.14$Method=="CMLM+SUPER",]
t.16 <- t.15[!(t.15$Ind.SNP=="F116+36088im"),]
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.detected"))
t.16$Ind.SNP <- factor(t.16$Ind.SNP,levels = c("F116+36088im", "F106+4185", "F116+3077", "F122+2272"))
t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
MLMSUPER.CMLMSUPER<- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) + facet_grid(~Method+Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=0), axis.text.y=element_text(size=12))+
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="MLM+SUPER, CMLM+SUPER")

###F122+2272
t.15 <- t.14[!(t.14$Method=="CMLM+SUPER" | t.14$Method=="MLM+SUPER"),]
t.16 <- t.15[t.15$Ind.SNP=="F122+2272",]
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.detected"))
t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
F2272.Ind.SNP <- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) + facet_grid(~Method)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=0), axis.text.y=element_text(size=12))+
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="F122+2272")

###for all of multiple methods
t.14 <- melt(t.13, id.vars = c("Ind.SNP","Method","Trait.name"))
str(t.14)
t.16 <- t.14[!(t.14$Method=="CMLM+SUPER"|t.14$Method=="MLM+SUPER"),]
#levels(as.factor(t.15$Method))
#t.16 <- t.15[t.15$Ind.SNP=="F122+2272",]
t.16$variable<- factor(t.16$variable, levels= c("SNPs.hit.with.flowering.genes", "SNPs.Total.detected"))
write.csv(t.16, file="data/SNPs.M.csv")
t.16$Method <- factor(t.16$Method,levels = c("FASTmrMLM", "mrMLM", "pLARmEB", "ISIS EM-BLASSO","FarmCPU","FASTmrEMMA","pKWmEB"))
t.16[t.16 == 0] <- NA
t.16$Ind.SNP <- factor(t.16$Ind.SNP,levels = c("F116+36088im", "F106+4185", "F116+3077", "F122+2272"))
t.16$Trait.name<- factor(t.16$Trait.name,levels = c("D", "W", "GW", "M"))
allmultiple <- ggplot(data=t.16, aes(x=Trait.name, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=t.16$value),size=4, position = position_stack(vjust = 0.5)) + facet_grid(Method~Ind.SNP)+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12), axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=0), axis.text.y=element_text(size=12))+
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="all")


### this is used for the outpur together 
library(ggplot2)
#install.packages("ggpubr")
library(ggpubr)
#theme_set(theme_pubr())
#install.packages("gridExtra")
library(gridExtra)

#fig.all <- ggarrange(F36088im.Ind.SNP, F4185.Ind.SNP,F3077.Ind.SNP,F2272.Ind.SNP,labels = c("A", "B","C","D","E","F"), ncol = 1, nrow = 6)

#figure106 <- annotate_figure(figure106,top = text_grob("Combination of F106+4185", color = "red", face = "bold", size = 14),left = text_grob("The number of QTNs", color = "red", rot = 90))
#ggexport(figure106, filename = "flo.8.4/4185.Multiple.pdf")

                 
