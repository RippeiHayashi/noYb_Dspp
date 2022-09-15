### this script is to measure the in trans linkage values between the whole ovary small RNAs and Piwi and Aubergine bound small RNAs from D. eugracilis ovaries
library("reshape2")

analysis="<directory where the output files are kept>"
CLASS="<all or unique>"

table <- read.table("${analysis}/Deug/Deug_${CLASS}_all_m9_dump.tab", header=F)
table_d=dcast(table,V1~V3+V4,value.var="V2")
table_d[is.na(table_d)] <- 0

PRODUCTSUM_Aub <- data.frame()
PRODUCTSUM_Aub[1,1] <- sum(table_d$'g1g9_Deug-ovaries-ox-sRNA-rep2'*table_d$'g2g10_revComp_Deug-Aub-IP-rep1')
PRODUCTSUM_Aub[2,1] <- sum(table_d$'g2g10_revComp_Deug-Aub-IP-rep1'*table_d$'last9_Deug-ovaries-ox-sRNA-rep2')
PRODUCTSUM_Aub[1,2] <- "Deug_ovary-Aub-IP"
PRODUCTSUM_Aub[2,2] <- "Deug_ovary-Aub-IP"
rownames(PRODUCTSUM_Aub) <- c("g1g9_g2g10_revComp","g2g10_revComp_last9")
write.table(PRODUCTSUM_Aub, file="${analysis}/Deug/Deug_${CLASS}_ovary-Aub-IP_linkage.txt", quote=F, col.names=F, row.names = T)

PRODUCTSUM_Piwi <- data.frame()
PRODUCTSUM_Piwi[1,1] <- sum(table_d$'g1g9_Deug-ovaries-ox-sRNA-rep2'*table_d$'g2g10_revComp_Deug-Piwi-IP-rep1')
PRODUCTSUM_Piwi[2,1] <- sum(table_d$'g2g10_revComp_Deug-Piwi-IP-rep1'*table_d$'last9_Deug-ovaries-ox-sRNA-rep2')
PRODUCTSUM_Piwi[1,2] <- "Deug_ovary-Piwi-IP"
PRODUCTSUM_Piwi[2,2] <- "Deug_ovary-Piwi-IP"
rownames(PRODUCTSUM_Piwi) <- c("g1g9_g2g10_revComp","g2g10_revComp_last9")
write.table(PRODUCTSUM_Piwi, file="${analysis}/Deug/Deug_${CLASS}_ovary-Piwi-IP_linkage.txt", quote=F, col.names=F, row.names = T)
