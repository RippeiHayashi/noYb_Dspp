### this script is to measure the in trans linkage values
### it is designed to be used in a bash script with the command 'Rscript'

#get in variables from bash
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

library(reshape2)

v=c(DIRECTORY,"/",LIB,"/jellyfish/",LIB,"_",CLASS,"_all_m9_dump.tab")
vname=paste(v,collapse="")
table=read.table(vname, header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0

PRODUCTSUM <- data.frame()
PRODUCTSUM[1,1] <- sum(table_d$'g1g9'*table_d$'last9')
PRODUCTSUM[2,1] <- sum(table_d$'g1g9'*table_d$'g2g10_revComp')
PRODUCTSUM[3,1] <- sum(table_d$'g2g10_revComp'*table_d$'last9')

rownames(PRODUCTSUM) <- c("g1g9_last9","g1g9_g2g10_revComp","g2g10_revComp_last9")

t=c(DIRECTORY,"/",LIB,"/jellyfish/",LIB,"_",CLASS,"_all_m9.linkage.txt")
tname=paste(t,collapse="")
write.table(PRODUCTSUM, file=tname, quote=F, row.names=T, col.names=F)
