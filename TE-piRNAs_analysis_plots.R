### this script is to plot the analyses of TE mapping piRNAs

### R libraries
library("ggplot2")
library("reshape2")
library("gridExtra")
pdf.options(useDingbats = FALSE)


### egg piRNAs vs ovary piRNAs from Dmel --- START ---
table <- read.table("Dmel_98TE_mappers.counts", header=F)
table.d <- dcast(table,V1~V3,value.var="V2")

p<-ggplot(table.d)+
geom_point(aes(x=`w1118-ovaries-ox-sRNA-rep2`,
               y=`w1118-egg-ox-sRNA-rep2`),size=2,shape=1,alpha=0.7,colour="gray")+
labs(title="TE piRNAs per 1 Mio genome mappers >22nt, Dmel oxidised libs, ovaries vs eggs",
     x="w1118-ovaries-ox-sRNA-rep2 (ovaries)", y="w1118-egg-ox-sRNA-rep2 (eggs)")+
scale_x_log10(limits=c(50,100000))+
scale_y_log10(limits=c(50,100000))+
coord_fixed()+
theme_bw()

pdf(file="Dmel_ovary-vs-egg_TE_RH_2022-09-07.pdf",width=8,height=8)
p
dev.off()
### egg piRNAs vs ovary piRNAs from Dmel --- END ---


### egg piRNAs vs ovary piRNAs from Dpse --- START ---
table <- read.table("Dpse_109TE_mappers.counts", header=F)
table.d <- dcast(table,V1~V3,value.var="V2")

p<-ggplot(table.d)+
geom_point(aes(x=`Dpse-ovaries-ox-sRNA-rep2`,
               y=`Dpse-egg-ox-sRNA-rep2`),size=2,shape=1,alpha=0.7,colour="gray")+
labs(title="TE piRNAs per 1 Mio genome mappers >22nt, Dpse oxidised libs, ovaries vs eggs",
     x="Dpse-ovaries-ox-sRNA-rep2 (ovaries)", y="Dpse-egg-ox-sRNA-rep2 (eggs)")+
scale_x_log10(limits=c(50,100000))+
scale_y_log10(limits=c(50,100000))+
coord_fixed()+
theme_bw()

pdf(file="Dpse_ovary-vs-egg_TE_RH_2022-05-05.pdf",width=8,height=8)
p
dev.off()
### egg piRNAs vs ovary piRNAs from Dpse --- END ---


### egg piRNAs vs ovary piRNAs from Deug --- START ---
table <- read.table("Deug_94TE_mappers.counts", header=F)
table.d <- dcast(table,V1~V3,value.var="V2")

p<-ggplot(table.d)+
geom_point(aes(x=`Deug-ovaries-ox-sRNA-rep2`,
               y=`Deug-egg-sRNA`),size=2,shape=1,alpha=0.7,colour="gray")+
labs(title="TE piRNAs per 1 Mio genome mappers >22nt, Deug oxidised libs, ovaries vs eggs",
     x="Deug-ovaries-ox-sRNA-rep2 (ovaries)", y="Deug-egg-sRNA (eggs)")+
scale_x_log10(limits=c(30,50000))+
scale_y_log10(limits=c(30,50000))+
coord_fixed()+
theme_bw()

pdf(file="Deug_ovary-vs-egg_TE_RH_2022-09-06.pdf",width=8,height=8)
p
dev.off()
### egg piRNAs vs ovary piRNAs from Deug --- END ---


### TE piRNAs AS S bias --- START ---
table <- read.table("TE_S_AS_bias.stats", header=F)
colnames(table) <- c("TE","AS.ratio","library")

table$library <- factor(table$library,levels = c("w1118-ovaries-ox-sRNA-rep2", "Deug-ovaries-ox-sRNA-rep2", "Dpse-ovaries-ox-sRNA-rep2"))

pdf(file="TE_S_AS_bias_Deug_Dpse_w1118_2022-09-07.pdf", height=6, width=5)
ggplot(table, aes(x=library, y=AS.ratio)) +
  geom_boxplot()+
  geom_jitter(size=2,shape=1, position=position_jitter(0.2)) +
  labs(title="TE AS piRNAs", y="AS/AS+S")
dev.off()

table.Deug <- subset(table, table$library == "Deug-ovaries-ox-sRNA-rep2")
table.w1118 <- subset(table, table$library == "w1118-ovaries-ox-sRNA-rep2")
table.Dpse <- subset(table, table$library == "Dpse-ovaries-ox-sRNA-rep2")

wilcox.test(table.Deug$AS.ratio, table.w1118$AS.ratio, alternative = "two.sided")
W = 7288, p-value < 2.2e-16
wilcox.test(table.Dpse$AS.ratio, table.w1118$AS.ratio, alternative = "two.sided")
W = 6640, p-value = 1.224e-06
### TE piRNAs AS S bias --- END ---



### 5' 5' pingpong linkage --- START ---
table <- read.table("pingpong.stats.selected", header=F)
colnames(table) <- c("TE","pingpong","library")

table$library <- factor(table$library,levels = c("w1118-ovaries-ox-sRNA-rep2", "Deug-ovaries-ox-sRNA-rep2", "Dpse-ovaries-ox-sRNA-rep2"))

pdf(file="pingpong_linkage_Deug_Dpse_w1118_2022-09-07.pdf", height=6, width=5)
ggplot(table, aes(x=library, y=pingpong)) +
  geom_boxplot()+
  ylim(c(0,0.6))+
  geom_jitter(size=2,shape=1, position=position_jitter(0.2)) +
  labs(title="5' 5' pingpong linkage", y="linkage")
dev.off()

table.Deug <- subset(table, table$library == "Deug-ovaries-ox-sRNA-rep2")
table.w1118 <- subset(table, table$library == "w1118-ovaries-ox-sRNA-rep2")
table.Dpse <- subset(table, table$library == "Dpse-ovaries-ox-sRNA-rep2")

wilcox.test(table.Deug$pingpong, table.w1118$pingpong, alternative = "two.sided")
W = 0, p-value = 2.495e-07
wilcox.test(table.Dpse$pingpong, table.w1118$pingpong, alternative = "two.sided")
W = 7172, p-value < 2.2e-16
### 5' 5' pingpong linkage --- END ---



### 3' 5' phasing linkage --- START ---
table <- read.table("phasing.stats.selected", header=F)
colnames(table) <- c("TE","phasing","library")

table$library <- factor(table$library,levels = c("Deug-ovaries-ox-sRNA-rep2", "Deug-Piwi-IP-rep1", "Deug-Aub-IP-rep1"))

pdf(file="phasing_linkage_Deug_IP_2022-09-07.pdf", height=6, width=5)
ggplot(table, aes(x=library, y=phasing)) +
  geom_boxplot()+
  ylim(c(0,0.2))+
  geom_jitter(size=2,shape=1, position=position_jitter(0.2)) +
  labs(title="3' 5' TE AS piRNAs phasing linkage", y="linkage")
dev.off()

table.Deug <- subset(table, table$library == "Deug-ovaries-ox-sRNA-rep2")
table.Deug.Aub <- subset(table, table$library == "Deug-Aub-IP-rep1")
table.Deug.Piwi <- subset(table, table$library == "Deug-Piwi-IP-rep1")

wilcox.test(table.Deug$phasing, table.Deug.Aub$phasing, alternative = "two.sided")
W = 4665, p-value = 8.143e-05
wilcox.test(table.Deug$phasing, table.Deug.Piwi$phasing, alternative = "two.sided")
W = 2326, p-value = 0.0003053
### 3' 5' phasing linkage --- END ---


### quantifying TE piRNAs normalised to miRNA abundance --- START ---
table <- read.table("miRNA_TE-mappers.stats", header=F)
colnames(table) <- c("type","reads","library")
table$type <- factor(table$type, levels = c("others", "TE_S_AS", "TE_S", "TE_AS"))

table$library <- factor(table$library,levels = c("SRR3715418", "D-eugracilis-sRNA", "Dpse-sRNA-rep1"))

pdf(file="TE_S_AS_piRNAs_per-miRNAs_2022-09-07.pdf",width=5,height=6)
ggplot(table,aes(x=library,y=reads,fill=type))+
  labs(title="abundance of TE piRNAs per 1 Mio miRNAs")+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
### quantifying TE piRNAs normalised to miRNA abundance --- END ---


### Aub-IP or Piwi-IP per whole ovary per 0.5kb tiles --- START ---
table <- read.table("deug_0.5kb_tiles.counts.Aub-Piwi-IP.table", header=F)
colnames(table) <- c("tile","IP.ovary.ratio","library")

table$library <- factor(table$library,levels = c("Piwi-IP-per-ovary_total", "Piwi-IP-per-ovary_uni", "Aub-IP-per-ovary_total", "Aub-IP-per-ovary_uni"))

pdf(file="Aub-Piwi-IP_vs_ovary_sRNA-tiles_Deug_2022-09-07.pdf", height=5, width=5)
ggplot(table, aes(x=library, y=IP.ovary.ratio)) +
  geom_boxplot()+
  geom_jitter(size=2,shape=1, position=position_jitter(0.2)) +
  scale_y_continuous(trans='log2', limits=c(0.0625,8))+
  labs(title="Enrichment of somatic piRNA tiles in Aub/Piwi-IP", y="IP/whole ovary")
dev.off()

table.Aub <- subset(table, table$library == "Aub-IP-per-ovary_total")
table.Aub.uni <- subset(table, table$library == "Aub-IP-per-ovary_uni")
table.Piwi <- subset(table, table$library == "Piwi-IP-per-ovary_total")
table.Piwi.uni <- subset(table, table$library == "Piwi-IP-per-ovary_uni")

wilcox.test(table.Aub$IP.ovary.ratio, table.Aub.uni$IP.ovary.ratio, alternative = "two.sided")
W = 36495, p-value = 3.064e-07
wilcox.test(table.Piwi$IP.ovary.ratio, table.Piwi.uni$IP.ovary.ratio, alternative = "two.sided")
W = 1262.5, p-value = 2.196e-11
### Aub-IP or Piwi-IP per whole ovary per 0.5kb tiles --- END ---



### TE mapping reads AS_1st size and nucleotide profile --- START ---
LIB="Deug-sRNA-seq-ox"
v=c(LIB,"_1st-nucleotide_length.TE_AS.table")
vname=paste(v,collapse="")
table_AS_1st=read.table(vname,header=F)

colnames(table_AS_1st)=c("size","nucleotide","reads","library")

mains_AS_1st=c(LIB,", TE_3MM_mappers_AS_1st size profile")
mains_AS_1st=paste(mains_AS_1st,collapse="")
AS_1st <- ggplot(table_AS_1st, aes(x=size,y=reads,fill=nucleotide))+
  labs(title=mains_AS_1st)+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### TE S_1st size and nucleotide profile
v=c(LIB,"_1st-nucleotide_length.TE_S.table")
vname=paste(v,collapse="")
table_S_1st=read.table(vname,header=F)

colnames(table_S_1st)=c("size","nucleotide","reads","library")

mains_S_1st=c(LIB,", TE_3MM_mappers_S_1st size profile")
mains_S_1st=paste(mains_S_1st,collapse="")
S_1st <- ggplot(table_S_1st, aes(x=size,y=reads,fill=nucleotide))+
  labs(title=mains_S_1st)+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

### combine two plots in a single pdf
p=c(LIB,"_TE_3MM_mappers_AS-and-S_1st_size_profile.pdf")
pname=paste(p,collapse="")

pdf(file=pname,width=12,height=6)
grid.arrange(AS_1st,S_1st, ncol=2)
dev.off()
### TE mapping reads AS_1st size and nucleotide profile --- END ---
