### this script is to measure the phasing linkage values
### it is designed to be used in a bash script with the command 'Rscript'

### get in variables from bash
args  =  commandArgs(TRUE);
argmat  =  sapply(strsplit(args, "="), identity)

for (i in seq.int(length=ncol(argmat))) {
  assign(argmat[1, i], argmat[2, i])
}

### R libraries
library(reshape2)
pdf.options(useDingbats = FALSE)

### import the count table
v=c(DIRECTORY,"/",LIB,"_",TE,"_3MM_mappers.counts")
vname=paste(v,collapse="")
table=read.table(vname,header=F)
table_d <- dcast(table,V1~V3, value.var="V2")

### set up the linkage position
LENGTH=length(table_d[,1])
START=-10
END=9
POSITION=-1
table_a=table_d$minus_3end
table_b=table_d$minus_5end

### make variables numeric
length=as.numeric(LENGTH)
start=as.numeric(START)
end=as.numeric(END)
link=as.numeric(POSITION)
start.position_B=1-start
end.position_B=length-end
link.position=link-start+1
span=end-start+1

if (link > 0){
  start.position_R=1
  end.position_R=length-link
} else {
  start.position_R=1-link
  end.position_R=length
}

### define a data.frame and a vector
offset <- data.frame()
sum_offset = NULL

### core of the linkage calculation
for (i in 1:span) {
  for (j in start.position_B:end.position_B) {
    offset[j,i] = table_a[j]*table_b[j+start+i-1]/sum(table_b[(j+start):(j+end)])
  }
}

offset[is.na(offset)] <- 0

for (i in 1:span) {
  sum_offset[i] = sum(offset[,i])/sum(offset)
}

mainname="pingpong, minus_3end vs minus_5end"
p=c(DIRECTORY,"/",LIB,"_",TE,"_3MM_mappers.pingpong.minus_3end.minus_5end.pdf")
pname=paste(p,collapse="")

pdf(file=pname,width=7,height=6)
plot(sum_offset,
     main=mainname,
     ylim=c(0,0.2))
dev.off()

#### calculate the Z-score
sum_offset2=sum_offset[c(1:(link.position-1),(link.position+1):span)]
m=mean(sum_offset2)
s=sd(sum_offset2)
z=(sum_offset[link.position]-m)/s
sum_offset3=c(z,sum_offset)

t=c(DIRECTORY,"/",LIB,"_",TE,"_3MM_mappers.pingpong.minus_3end.minus_5end.table")
tname=paste(t,collapse="")
write.table(sum_offset3, file=tname, quote=F, row.names=F)
