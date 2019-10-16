###import the data 

setwd("/Users/yonglizhao/Documents/R-corde for miscanthus project/Miscanthus")
trait <- read.csv("Result.8.3/Total.traits.overlap.csv")

trait$SNPs.hit.with.flowering.genes.[trait$SNPs.hit.with.flowering.genes.== 0] <- NA
str(trait)
#install.packages("ggrepel")
#install.packages("reshape2")
#install.packages("reshape")
#library(ggrepel)
library(reshape2)
library(reshape)
trait.1 <- melt(trait, id.vars = c("Type1","Scale.of.measurement"))
str(trait.1)

##this one is for PC1 ith four temporal measure scales 
trait.2 <- trait.1[trait.1$Type1=="DAY"|trait.1$Type1=="Week"|trait.1$Type1=="Group Week"|trait.1$Type1=="Month",]
trait.2 <- droplevels(trait.2)
levels(trait.2$variable)
levels(trait.2$Scale.of.measurement)
###

trait.2$Scale.of.measurement<- factor(trait.2$Scale.of.measurement,levels = c("D", "PC1.D","D+PC1.D","W","PC1.W","W+PC1.W","GW","PC1.GW","GW+PC1.GW","M","PC1.M","M+PC1.M"))
trait.2$variable<- factor(trait.2$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
OverPCA <- ggplot(data=trait.2, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=trait.2$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="PC1 and four temporal measure scales ")


###this one for Heading_1
###checking the levels
levels(trait.1$Type1)
trait.h <- trait.1[trait.1$Type1=="Heading_1",]
trait.h <- droplevels(trait.h)
levels(trait.h$variable)
levels(trait.h$Scale.of.measurement)
###
trait.h$Scale.of.measurement<- factor(trait.h$Scale.of.measurement,levels = c("D", "W","GW","M","D+W","D+GW","D+M","W+GW","GW+M","W+M","D+W+GW","D+GW+M","D+W+M","W+GW+M","D+W+GW+M"))
trait.h$variable<- factor(trait.h$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
h <- ggplot(data=trait.h, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=trait.h$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
       labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="First Heading Time")


###this one for Heading_50
###checking the levels
levels(trait.1$Type1)
trait.h.5 <- trait.1[trait.1$Type1=="Heading_50",]
trait.h.5 <- droplevels(trait.h.5)
levels(trait.h.5$variable)
levels(trait.h.5$Scale.of.measurement)
###
trait.h.5$Scale.of.measurement<- factor(trait.h.5$Scale.of.measurement,levels = c("D", "W","GW","M","D+W","D+GW","D+M","W+GW","GW+M","W+M","D+W+GW","D+GW+M","D+W+M","W+GW+M","D+W+GW+M"))
trait.h.5$variable<- factor(trait.h.5$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
h.5 <- ggplot(data=trait.h.5, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=trait.h.5$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
        labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="Half Heading Time")

###this one for Flowering_1
###checking the levels
levels(trait.1$Type1)
trait.f <- trait.1[trait.1$Type1=="Flowering_1",]
trait.f <- droplevels(trait.f)
levels(trait.f$variable)
levels(trait.f$Scale.of.measurement)
###
trait.f$Scale.of.measurement<- factor(trait.f$Scale.of.measurement,levels = c("D", "W","GW","M","D+W","D+GW","D+M","W+GW","GW+M","W+M","D+W+GW","D+GW+M","D+W+M","W+GW+M","D+W+GW+M"))
trait.f$variable<- factor(trait.f$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
f <- ggplot(data=trait.f, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=trait.f$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="First Flowering Time")

###this one for Flowering_50
###checking the levels
levels(trait.1$Type1)
trait.f.5 <- trait.1[trait.1$Type1=="Flowering_50",]
trait.f.5 <- droplevels(trait.f.5)
levels(trait.f.5$variable)
levels(trait.f.5$Scale.of.measurement)
###
trait.f.5$Scale.of.measurement<- factor(trait.f.5$Scale.of.measurement,levels = c("D", "W","GW","M","D+W","D+GW","D+M","W+GW","GW+M","W+M","D+W+GW","D+GW+M","D+W+M","W+GW+M","D+W+GW+M"))
trait.f.5$variable<- factor(trait.f.5$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
f.5 <- ggplot(data=trait.f.5, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=trait.f.5$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
        labs(x="D=Day,W=Week,GW=Group Week,M=Month", y="The number of QTNs",title="Half Flowering Time")

###this one for PC1
###checking the levels
levels(trait.1$Type1)
traitPC <- trait.1[trait.1$Type1=="PC1",]
traitPC <- droplevels(traitPC)
levels(traitPC$variable)
levels(traitPC$Scale.of.measurement)
###
traitPC$Scale.of.measurement<- factor(traitPC$Scale.of.measurement,levels = c("D", "W","GW","M","D+W","D+GW","D+M","W+GW","GW+M","W+M","D+W+GW","D+GW+M","D+W+M","W+GW+M","D+W+GW+M"))
traitPC$variable<- factor(traitPC$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
PC <- ggplot(data=traitPC, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=traitPC$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
       labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="PC1 of Flowering Time")

###this one for Total.F
###checking the levels
levels(trait.1$Type1)
traitTF <- trait.1[trait.1$Type1=="Total.F",]
traitTF <- droplevels(traitTF)
levels(traitTF$variable)
levels(traitTF$Scale.of.measurement)
###
traitTF$Scale.of.measurement<- factor(traitTF$Scale.of.measurement,levels = c("D", "W","GW","M","D+W","D+GW","D+M","W+GW","GW+M","W+M","D+W+GW","D+GW+M","D+W+M","W+GW+M","D+W+GW+M"))
traitTF$variable<- factor(traitTF$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
TF <- ggplot(data=traitTF, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=traitTF$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
       labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="Total of Flowering Time")


###this one for Total.F
###checking the levels
levels(trait.1$Type1)
traitTF <- trait.1[trait.1$Type1=="Total.F",]
traitTF <- droplevels(traitTF)
levels(traitTF$variable)
levels(traitTF$Scale.of.measurement)
###
traitTF$Scale.of.measurement<- factor(traitTF$Scale.of.measurement,levels = c("D", "W","GW","M","D+W","D+GW","D+M","W+GW","GW+M","W+M","D+W+GW","D+GW+M","D+W+M","W+GW+M","D+W+GW+M"))
traitTF$variable<- factor(traitTF$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
TF <- ggplot(data=traitTF, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="D:Day  W:Week  GW:Group Week  M:Month",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=traitTF$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="D=Day,W=Week,GW=Group Week,M=Month",y="The number of QTNs",title="Four temporal measure scales")



###PS
###checking the levels
levels(trait.1$Type1)
traitPS <- trait.1[trait.1$Type1=="PS",]
traitPS <- droplevels(traitPS)
levels(traitPS$variable)
levels(traitPS$Scale.of.measurement)
###
traitPS$Scale.of.measurement<- factor(traitPS$Scale.of.measurement,levels = c("PS_N","PS_Y","PS_Y and PS_N"))
traitPS$variable<- factor(traitPS$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
PS <- ggplot(data=traitPS, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="PS_N:No Population Structure\nPS_Y:Population Structure\nPS_Y and PS_N:Overlap",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=traitPS$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12), axis.text.y=element_text(size=12))+
  labs(x="PS_N=Population Structure not in Model,PS_Y=Population Structure in Model,PS_Y and PS_N=Overlap",y="The number of QTNs",title="Population structure with or without")

###Ind.SNP
###checking the levels
levels(trait.1$Type1)
traitInd.SNP <- trait.1[trait.1$Type1=="Ind.SNP",]
traitInd.SNP <- droplevels(traitInd.SNP)
levels(traitInd.SNP$variable)
levels(traitInd.SNP$Scale.of.measurement)

levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im"] <- "116im"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185"] <- "106"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+3077"] <- "116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F122+2272"] <- "122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F116+3077"] <- "106+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F122+2272"] <- "106+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+3077;F122+2272"] <- "116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F116+3077;F122+2272"] <- "106+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F106+4185"] <- "116im+106"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F116+3077"] <- "116im+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F122+2272"] <- "116im+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F106+4185;F116+3077"] <- "116im+106+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F116+3077;F122+2272"] <- "116im+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)== "F116+36088im;F106+4185;F106+4185;F116+3077;F122+2272"] <- "116im+106+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)== "F116+36088im+all.non.missing"] <- "116im+Non.Missing"
levels(traitInd.SNP$Scale.of.measurement)
traitInd.SNP <- traitInd.SNP[!(traitInd.SNP$Scale.of.measurement=="116im+106"|
                                 traitInd.SNP$Scale.of.measurement=="116im+116"|traitInd.SNP$Scale.of.measurement=="116im+122"|traitInd.SNP$Scale.of.measurement=="116im+106+116"
                                |traitInd.SNP$Scale.of.measurement=="116im+116+122"|traitInd.SNP$Scale.of.measurement=="116im+106+116+122"),]

###
traitInd.SNP$Scale.of.measurement<- factor(traitInd.SNP$Scale.of.measurement,levels = c("116im","106","116","122", "106+116","106+122",
                                                                                        "116+122","106+116+122", "116im+106","116im+Non.Missing"))

traitInd.SNP$variable<- factor(traitInd.SNP$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
Ind.SNP <- ggplot(data=traitInd.SNP, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="116im:F116+36088im\n106:F106+4185\n116:F116+3077\n122:F122+2272",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=traitInd.SNP$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="116im=F116+36088im,106=F106+4185,116=F116+3077,122=F122+2272", y="The number of QTNs",title="Four Combinations of accessions and SNPs")




###Ind.SNP
###checking the levels
levels(trait.1$Type1)
traitInd.SNP <- trait.1[trait.1$Type1=="Ind.SNP",]
traitInd.SNP <- droplevels(traitInd.SNP)
levels(traitInd.SNP$variable)
levels(traitInd.SNP$Scale.of.measurement)

levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im"] <- "116im"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185"] <- "106"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+3077"] <- "116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F122+2272"] <- "122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F116+3077"] <- "106+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F122+2272"] <- "106+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+3077;F122+2272"] <- "116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F116+3077;F122+2272"] <- "106+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F106+4185"] <- "116im+106"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F116+3077"] <- "116im+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F122+2272"] <- "116im+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F106+4185;F116+3077"] <- "116im+106+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F116+3077;F122+2272"] <- "116im+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)== "F116+36088im;F106+4185;F106+4185;F116+3077;F122+2272"] <- "116im+106+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)== "F116+36088im+all.non.missing"] <- "116im+Non.Missing"
levels(traitInd.SNP$Scale.of.measurement)
traitInd.SNP <- traitInd.SNP[!(traitInd.SNP$Scale.of.measurement=="116im+106"|
                                 traitInd.SNP$Scale.of.measurement=="116im+116"|traitInd.SNP$Scale.of.measurement=="116im+122"|traitInd.SNP$Scale.of.measurement=="116im+106+116"
                               |traitInd.SNP$Scale.of.measurement=="116im+116+122"|traitInd.SNP$Scale.of.measurement=="116im+106+116+122"),]

###
traitInd.SNP$Scale.of.measurement<- factor(traitInd.SNP$Scale.of.measurement,levels = c("116im","106","116","122", "106+116","106+122",
                                                                                        "116+122","106+116+122", "116im+106","116im+Non.Missing"))

traitInd.SNP$variable<- factor(traitInd.SNP$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
Ind.SNP <- ggplot(data=traitInd.SNP, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="116im:F116+36088im 106:F106+4185 116:F116+3077 122:F122+2272",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=traitInd.SNP$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position = "bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="116im=F116+36088im,106=F106+4185,116=F116+3077,122=F122+2272", y="The number of QTNs",title="Four Combinations of accessions and SNPs")


###Ind.SNP
###checking the levels
levels(trait.1$Type1)
traitInd.SNP <- trait.1[trait.1$Type1=="Ind.SNP",]
traitInd.SNP <- droplevels(traitInd.SNP)
levels(traitInd.SNP$variable)
levels(traitInd.SNP$Scale.of.measurement)

levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im"] <- "116im"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185"] <- "106"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+3077"] <- "116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F122+2272"] <- "122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F116+3077"] <- "106+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F122+2272"] <- "106+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+3077;F122+2272"] <- "116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F106+4185;F116+3077;F122+2272"] <- "106+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F106+4185"] <- "116im+106"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F116+3077"] <- "116im+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F122+2272"] <- "116im+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F106+4185;F116+3077"] <- "116im+106+116"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)=="F116+36088im;F116+3077;F122+2272"] <- "116im+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)== "F116+36088im;F106+4185;F106+4185;F116+3077;F122+2272"] <- "116im+106+116+122"
levels(traitInd.SNP$Scale.of.measurement)[levels(traitInd.SNP$Scale.of.measurement)== "F116+36088im+all.non.missing"] <- "116im+Non.Missing"
levels(traitInd.SNP$Scale.of.measurement)
traitInd.SNP <- traitInd.SNP[!(traitInd.SNP$Scale.of.measurement=="116im"|traitInd.SNP$Scale.of.measurement=="116im+106"|
                                 traitInd.SNP$Scale.of.measurement=="116im+116"|traitInd.SNP$Scale.of.measurement=="116im+122"|traitInd.SNP$Scale.of.measurement=="116im+106+116"
                               |traitInd.SNP$Scale.of.measurement=="116im+116+122"|traitInd.SNP$Scale.of.measurement=="116im+106+116+122"|traitInd.SNP$Scale.of.measurement=="116im+Non.Missing"),]

###
traitInd.SNP$Scale.of.measurement<- factor(traitInd.SNP$Scale.of.measurement,levels = c("106","116","122", "106+116","106+122",
                                                                                        "116+122","106+116+122"))

traitInd.SNP$variable<- factor(traitInd.SNP$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
Ind.SNP.nonm <- ggplot(data=traitInd.SNP, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="106:F106+4185\n116:F116+3077\n122:F122+2272",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=traitInd.SNP$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(x="116im=F116+36088im,106=F106+4185,116=F116+3077,122=F122+2272", y="The number of QTNs",title="Non Missing Combinations of accessions and SNPs")





###Ind.SNP
###checking the levels
levels(trait.1$Type1)
traitInd.SNP.T <- trait.1[trait.1$Type1=="Ind.SNP.T",]
traitInd.SNP.T <- droplevels(traitInd.SNP.T)
levels(traitInd.SNP.T$variable)
levels(traitInd.SNP.T$Scale.of.measurement)

###
traitInd.SNP.T$Scale.of.measurement<- factor(traitInd.SNP.T$Scale.of.measurement,levels = c("Imputed SNP","Non-missing SNPs","overlaps"))

traitInd.SNP.T$variable<- factor(traitInd.SNP.T$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))
Ind.SNP.T <- ggplot(data=traitInd.SNP.T, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.5)+
  scale_fill_discrete(name="SNPs",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_text(aes(label=traitInd.SNP.T$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12), axis.text.y=element_text(size=12))+
  labs(y="The number of QTNs",title="Imputed SNP, Non Misssing SNPs, and overlap")


###checking the levels
levels(trait.1$Type1)
trait.M.T <- trait.1[trait.1$Type1== "Method",]
trait.M.T <- droplevels(trait.M.T)
levels(trait.M.T$variable)
levels(trait.M.T$Scale.of.measurement)
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="M.FarmCPU"] <- "F"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="M.FarmCPU;M.mrMLM"] <- "F+M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="M.mrMLM"] <- "M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER"] <- "C.S"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER;S.MLM+SUPER"] <- "C.S+M.S"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER;S.MLM+SUPER;M.FarmCPU"] <- "C.S+M.S+F"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER;S.MLM+SUPER;M.FarmCPU;M.mrMLM"] <- "C.S+M.S+F+M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER;S.MLM+SUPER;M.mrMLM"] <- "C.S+M.S+M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.MLM+SUPER"] <- "M.S"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.MLM+SUPER;M.FarmCPU"] <- "M.S+F"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.MLM+SUPER;M.FarmCPU;M.mrMLM"] <- "M.S+F+M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.MLM+SUPER;M.mrMLM"] <- "M.S+M"
levels(trait.M.T$Scale.of.measurement)
trait.M.T$Scale.of.measurement<- factor(trait.M.T$Scale.of.measurement,levels = c("C.S","M.S","F","M","C.S+M.S","M.S+F","C.S+M.S+F","F+M","M.S+M","C.S+M.S+F+M","C.S+M.S+M","M.S+F+M"))

trait.M.T$variable<- factor(trait.M.T$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))


M.T <- ggplot(data=trait.M.T, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="C.S:CMLM+SUPER\nM.S:MLM+SUPER\nF:FarmCPU\nM:Methods in mrMLM.Package",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_label_repel(hjust = 0.5,aes(label=trait.M.T$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="The number of QTNs",title="Single-locus method,Multiple-locus method, and Overlap")



###checking the levels
levels(trait.1$Type1)
trait.M.T <- trait.1[trait.1$Type1== "Method",]
trait.M.T <- droplevels(trait.M.T)
levels(trait.M.T$variable)
levels(trait.M.T$Scale.of.measurement)
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="M.FarmCPU"] <- "F"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="M.FarmCPU;M.mrMLM"] <- "F+M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="M.mrMLM"] <- "M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER"] <- "C.S"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER;S.MLM+SUPER"] <- "C.S+M.S"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER;S.MLM+SUPER;M.FarmCPU"] <- "C.S+M.S+F"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER;S.MLM+SUPER;M.FarmCPU;M.mrMLM"] <- "C.S+M.S+F+M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.CMLM+SUPER;S.MLM+SUPER;M.mrMLM"] <- "C.S+M.S+M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.MLM+SUPER"] <- "M.S"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.MLM+SUPER;M.FarmCPU"] <- "M.S+F"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.MLM+SUPER;M.FarmCPU;M.mrMLM"] <- "M.S+F+M"
levels(trait.M.T$Scale.of.measurement)[levels(trait.M.T$Scale.of.measurement)=="S.MLM+SUPER;M.mrMLM"] <- "M.S+M"
levels(trait.M.T$Scale.of.measurement)
trait.M.T$Scale.of.measurement<- factor(trait.M.T$Scale.of.measurement,levels = c("C.S","M.S","F","M","C.S+M.S","M.S+F","C.S+M.S+F","F+M","M.S+M","C.S+M.S+F+M","C.S+M.S+M","M.S+F+M"))

trait.M.T$variable<- factor(trait.M.T$variable, levels= c("SNPs.hit.with.flowering.genes.", "SNPs.Total.detected."))

M.T <- ggplot(data=trait.M.T, aes(x=Scale.of.measurement, y=value, fill=variable)) +
  geom_bar(stat="identity",width=0.8)+
  scale_fill_discrete(name="C.S:CMLM+SUPER, M.S:MLM+SUPER\nF:FarmCPU, M:Methods in mrMLM.Package",labels=c("SNPs.Hit.With.Flowering.Genes","SNPs.Total.Detected"))+   
  geom_label_repel(hjust = 0.5,aes(label=trait.M.T$value),size=4,position = position_stack(vjust = 0.5))+
  theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5), plot.subtitle=element_text(size=12),axis.title.x=element_blank(),
        axis.text=element_text(size = 12),
        axis.title=element_text(size = 12),
        axis.title.y=element_text(size=12),
        plot.caption=element_text(size=12),
        strip.text = element_text(face="bold", size=12),
        legend.title = element_text(size=12, face="bold"),
        legend.text = element_text(size=12, face="bold"),
        legend.position="bottom",
        axis.text.x=element_text(face="bold",size=12,angle=45), axis.text.y=element_text(size=12))+
  labs(y="The number of QTNs",title="Single-locus methods,Multiple-locus methods, and Overlap")


##
#install.packages("grid")
#library(grid)
#my_text <- "C.S=CMLM+SUPER\nM.S=MLM+SUPER\nF=FarmCPU\nM=Methods from mrMLM.Package"
#my_grob = grid.text(my_text, x=0.7,  y=0.8, gp=gpar(fontsize=12, fontface="bold"))
#M.T + annotation_custom(my_grob)

figureoverlap1 <- ggarrange(h,h.5,labels = c("A", "B"), ncol = 1, nrow = 2, common.legend =TRUE, legend="bottom" )

figureoverlap2 <- ggarrange(f,f.5, labels = c("C", "D"), ncol = 1, nrow = 2,common.legend = TRUE, legend="bottom" )


figureoverlap7 <- ggarrange(h,f,labels = c("A", "B"), ncol = 1, nrow = 2, common.legend =TRUE, legend="bottom" )

figureoverlap8 <- ggarrange(h.5,f.5, labels = c("C", "D"), ncol = 1, nrow = 2,common.legend = TRUE, legend="bottom" )
#ggexport(figureoverlap1, filename = "flo.8.4/figureoverlap1.pdf")


figureoverlap3 <- ggarrange(PC, OverPCA, labels = c("E", "F"), ncol = 1, nrow = 2,common.legend = T, legend="bottom" )
figureoverlap4 <- ggarrange(PC, OverPCA, labels = c("E", "F"), ncol = 1, nrow = 2,common.legend = F, legend="bottom" )



figureoverlap6 <- ggarrange(Ind.SNP.T,PS, labels = c("K","L"), ncol = 1, nrow = 2,common.legend = F, legend="right")

#ggexport(figure116, filename = "flo.8.4/2272.Single.pdf")

