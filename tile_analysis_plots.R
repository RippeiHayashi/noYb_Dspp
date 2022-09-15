### this script is to plot the abundance of piRNAs uniquely mapping to the 0.5kb genomic tiles comparing the whole ovary and the embryonic piRNAs

### R libraries
library(ggplot2)
library(gridExtra)
library(reshape2)
pdf.options(useDingbats = FALSE)


### D.melanogaster 0.5kb genomic tiles
table=read.table("dmel_0.5kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0

table_d2 <- subset(table_d, table_d$`w1118-ovaries-ox-sRNA-rep2` > 0.1 & table_d$`w1118-egg-ox-sRNA-rep2` > 0.1)
p<-ggplot(table_d2)+
  geom_point(aes(x=`w1118-ovaries-ox-sRNA-rep2`,
                 y=`w1118-egg-ox-sRNA-rep2`),size=2,shape=1,alpha=0.5,colour="gray")+
  geom_point(data=subset(table_d2, ClusterFlam > 0.5),
             aes(x=`w1118-ovaries-ox-sRNA-rep2`,
                 y=`w1118-egg-ox-sRNA-rep2`),size=2,shape=1,alpha=0.6,colour="red")+
  geom_point(data=subset(table_d2, Cluster80F > 0.5),
             aes(x=`w1118-ovaries-ox-sRNA-rep2`,
                 y=`w1118-egg-ox-sRNA-rep2`),size=2,shape=1,alpha=0.6,colour="blue")+
  labs(title="Dmel piRNAs, ovaries vs eggs per 0.5kb tiles, ClusterFlam(red), Cluster80F(blue)",
       x="w1118-ovaries-ox-sRNA-rep2 / RPKM", y="w1118-egg-ox-sRNA-rep2 / RPKM")+
  scale_x_log10(limits=c(0.1,50000))+
  scale_y_log10(limits=c(0.1,50000))+
  coord_fixed()+
  theme_bw()

pdf(file="Dmel_ovary-vs-egg_0.5kbtiles_2022.pdf",width=8,height=8)
p
dev.off()


### D.pseudoobscura 0.5kb genomic tiles
table=read.table("dpse_0.5kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0

table_d2 <- subset(table_d, table_d$`Dpse-ovaries-ox-sRNA-rep2` > 0.1 & table_d$`Dpse-egg-ox-sRNA-rep2` > 0.1)
p<-ggplot(table_d2)+
  geom_point(aes(x=`Dpse-ovaries-ox-sRNA-rep2`,
                 y=`Dpse-egg-ox-sRNA-rep2`),size=2,shape=1,alpha=0.5,colour="gray")+
  geom_point(data=subset(table_d2, cluster3l_uni > 0.5),
             aes(x=`Dpse-ovaries-ox-sRNA-rep2`,
                 y=`Dpse-egg-ox-sRNA-rep2`),size=2,shape=1,alpha=0.6,colour="red")+
  geom_point(data=subset(table_d2, cluster40l_dual > 0.5),
             aes(x=`Dpse-ovaries-ox-sRNA-rep2`,
                 y=`Dpse-egg-ox-sRNA-rep2`),size=2,shape=1,alpha=0.6,colour="blue")+
  labs(title="Dpse piRNAs, ovaries vs eggs per 0.5kb tiles, cluster3l_uni(red), cluster40l_dual(blue)",
       x="Dpse-ovaries-ox-sRNA-rep2 / RPKM", y="Dpse-egg-ox-sRNA-rep2 / RPKM")+
  scale_x_log10(limits=c(0.1,50000))+
  scale_y_log10(limits=c(0.1,50000))+
  coord_fixed()+
  theme_bw()

pdf(file="Dpse_ovary-vs-egg_0.5kbtiles_2022.pdf",width=8,height=8)
p
dev.off()


### D.eugracilis 0.5kb genomic tiles
table=read.table("deug_0.5kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0

table_d2 <- subset(table_d, table_d$`Deug-ovaries-ox-sRNA-rep2` > 0.1 & table_d$`Deug-egg-sRNA` > 0.1)
p<-ggplot(table_d2)+
  geom_point(aes(x=`Deug-ovaries-ox-sRNA-rep2`,
                 y=`Deug-egg-sRNA`),size=2,shape=1,alpha=0.5,colour="gray")+
  geom_point(data=subset(table_d2, Cluster2449_uni > 0.5),
             aes(x=`Deug-ovaries-ox-sRNA-rep2`,
                 y=`Deug-egg-sRNA`),size=2,shape=1,alpha=0.6,colour="red")+
  geom_point(data=subset(table_d2, Cluster3031_uni > 0.5),
             aes(x=`Deug-ovaries-ox-sRNA-rep2`,
             y=`Deug-egg-sRNA`),size=2,shape=1,alpha=0.6,colour="orange")+
  geom_point(data=subset(table_d2, Cluster3703_uni > 0.5),
             aes(x=`Deug-ovaries-ox-sRNA-rep2`,
             y=`Deug-egg-sRNA`),size=2,shape=1,alpha=0.6,colour="purple")+
  geom_point(data=subset(table_d2, Cluster3378_dual > 0.5),
             aes(x=`Deug-ovaries-ox-sRNA-rep2`,
                 y=`Deug-egg-sRNA`),size=2,shape=1,alpha=0.6,colour="blue")+
  labs(title="Deug piRNAs, ovaries vs eggs per 0.5kb tiles, cluster3l_uni(red), cluster40l_dual(blue)",
       x="Deug-ovaries-ox-sRNA-rep2 / RPKM", y="Deug-egg-sRNA / RPKM")+
  scale_x_log10(limits=c(0.1,50000))+
  scale_y_log10(limits=c(0.1,50000))+
  coord_fixed()+
  theme_bw()

pdf(file="Deug_ovary-vs-egg_0.5kbtiles_2022.pdf",width=8,height=8)
p
dev.off()


### D.bifasciata 0.5kb genomic tiles
table=read.table("dbif_0.5kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0

table_d2 <- subset(table_d, table_d$`Dbif-sRNA-seq-ox-rep2` > 0.1 & table_d$`Dbif-egg-ox-sRNA` > 0.1)
p<-ggplot(table_d2)+
  geom_point(aes(x=`Dbif-sRNA-seq-ox-rep2`,
                 y=`Dbif-egg-ox-sRNA`),size=2,shape=1,alpha=0.5,colour="gray")+
  geom_point(data=subset(table_d2, CM1_137_uni > 0.5),
             aes(x=`Dbif-sRNA-seq-ox-rep2`,
                 y=`Dbif-egg-ox-sRNA`),size=2,shape=1,alpha=0.6,colour="red")+
  geom_point(data=subset(table_d2, CM1_99_uni > 0.5),
             aes(x=`Dbif-sRNA-seq-ox-rep2`,
                 y=`Dbif-egg-ox-sRNA`),size=2,shape=1,alpha=0.6,colour="orange")+
  geom_point(data=subset(table_d2, `CM1_132-135_dual` > 0.5),
             aes(x=`Dbif-sRNA-seq-ox-rep2`,
                 y=`Dbif-egg-ox-sRNA`),size=2,shape=1,alpha=0.6,colour="blue")+
  labs(title="Dbif piRNAs, ovaries vs eggs per 0.5kb tiles, CM1_137_uni(red), CM1_99_uni(orange), CM1_132-135_dual(blue)",
       x="Dbif-sRNA-seq-ox-rep2 / RPKM", y="Dbif-egg-ox-sRNA / RPKM")+
  scale_x_log10(limits=c(0.1,50000))+
  scale_y_log10(limits=c(0.1,50000))+
  coord_fixed()+
  theme_bw()

pdf(file="Dbif_ovary-vs-egg_0.5kbtiles_2022.pdf",width=8,height=8)
p
dev.off()


### D.mojavensis 0.5kb genomic tiles
table=read.table("dmoj_0.5kb_tiles.counts", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0

table_d2 <- subset(table_d, table_d$`Dmoj-ovaries-ox-sRNA` > 0.1 & table_d$`Dmoj-egg-ox-sRNA` > 0.1)
p<-ggplot(table_d2)+
  geom_point(aes(x=`Dmoj-ovaries-ox-sRNA`,
                 y=`Dmoj-egg-ox-sRNA`),size=2,shape=1,alpha=0.5,colour="gray")+
  geom_point(data=subset(table_d2, uni_667 > 0.5),
             aes(x=`Dmoj-ovaries-ox-sRNA`,
                 y=`Dmoj-egg-ox-sRNA`),size=2,shape=1,alpha=0.6,colour="red")+
  labs(title="Dbif piRNAs, ovaries vs eggs per 0.5kb tiles, uni_667(red)",
       x="Dmoj-ovaries-ox-sRNA / RPKM", y="Dmoj-egg-ox-sRNA / RPKM")+
  scale_x_log10(limits=c(0.1,50000))+
  scale_y_log10(limits=c(0.1,50000))+
  coord_fixed()+
  theme_bw()

pdf(file="Dmoj_ovary-vs-egg_0.5kbtiles_2022.pdf",width=8,height=8)
p
dev.off()
