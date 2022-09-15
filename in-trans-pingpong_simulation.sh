### this script is to simulate the in trans pingpong linkage values using artificially created libraries with varying number of pingpong pairs.
analysis="<directory where the output files are kept>"

mkdir -p ${analysis}/simulation

### step 1: generate libraries of 1 million numbers as a substitute of piRNA sequence, and incorporate the same numbers into different libraries to regard them as pingpong pairs --- START ---

### for x % of pingpong pairs, repeat this process for x times 100 times
### this is for 1% --- START ---
for j in {1..100}; do
read=`(shuf -i 1-262144 -n1)`
### for 10 groups of numbers
for k in {1..10}; do
### print the same number for 100 times, which means that each pingpong pair has 100 counts per 1 million reads
for i in {1..100}; do
echo $read >> ${analysis}/simulation/simulated.numbers.1Mio.1percent.${k}
done
done
done

### for the 10 groups generated above
for k in {1..10}; do
### add the rest of numbers that are randomly produced
for i in {1..990}; do
shuf -i 1-262144 -n1000 >> ${analysis}/simulation/simulated.numbers.1Mio.1percent.${k}
done
done
### this is for 1% --- END ---

### this is for 3% --- START ---
for j in {1..300}; do
read=`(shuf -i 1-262144 -n1)`
### for 10 groups of numbers
for k in {1..10}; do
### print the same number for 100 times, which means that each pingpong pair has 100 counts per 1 million reads
for i in {1..100}; do
echo $read >> ${analysis}/simulation/simulated.numbers.1Mio.3percent.${k}
done
done
done

### for the 10 groups generated above
for k in {1..10}; do
### add the rest of numbers that are randomly produced
for i in {1..970}; do
shuf -i 1-262144 -n1000 >> ${analysis}/simulation/simulated.numbers.1Mio.3percent.${k}
done
done
### this is for 3% --- END ---

### this is for 5% --- START ---
for j in {1..500}; do
read=`(shuf -i 1-262144 -n1)`
### for 10 groups of numbers
for k in {1..10}; do
### print the same number for 100 times, which means that each pingpong pair has 100 counts per 1 million reads
for i in {1..100}; do
echo $read >> ${analysis}/simulation/simulated.numbers.1Mio.5percent.${k}
done
done
done

### for the 10 groups generated above
for k in {1..10}; do
### add the rest of numbers that are randomly produced
for i in {1..950}; do
shuf -i 1-262144 -n1000 >> ${analysis}/simulation/simulated.numbers.1Mio.5percent.${k}
done
done
### this is for 5% --- END ---

### this is for 10% --- START ---
for j in {1..1000}; do
read=`(shuf -i 1-262144 -n1)`
### for 10 groups of numbers
for k in {1..10}; do
### print the same number for 100 times, which means that each pingpong pair has 100 counts per 1 million reads
for i in {1..100}; do
echo $read >> ${analysis}/simulation/simulated.numbers.1Mio.10percent.${k}
done
done
done

### for the 10 groups generated above
for k in {1..10}; do
### add the rest of numbers that are randomly produced
for i in {1..900}; do
shuf -i 1-262144 -n1000 >> ${analysis}/simulation/simulated.numbers.1Mio.10percent.${k}
done
done
### this is for 10% --- END ---

### this is for 0% --- START ---
### for the 10 groups generated above
for k in {1..10}; do
### add the rest of numbers that are randomly produced
for i in {1..1000}; do
shuf -i 1-262144 -n1000 >> ${analysis}/simulation/simulated.numbers.1Mio.0percent.${k}
done
done
### this is for 0% --- END ---
### step 1: generate libraries of 1 million numbers as a substitute of piRNA sequence, and incorporate the same numbers into different libraries to regard them as pingpong pairs --- END ---


### step 2: summarising the tables --- START ---
for TYPE in 0percent 1percent 3percent 5percent 10percent; do
for k in {1..10}; do
awk -v GROUP=${k} '{READS[$1]++} END {for(var in READS) print var,READS[var]/1000,"group_"GROUP
}' ${analysis}/simulation/simulated.numbers.1Mio.${TYPE}.${k} >> ${analysis}/simulation/simulated.numbers.1Mio.${TYPE}.10groups
done
done

### step 2-2: dcast the table
analysis="/g/data/lf10/rh1772/noYb_paper/analysis"
#R/4.0.0
library(reshape2)
table <- read.table("${analysis}/simulation/simulated.numbers.1Mio.${TYPE}.10groups", header=F)
table_d=dcast(table,V1~V3,value.var="V2")
table_d[is.na(table_d)] <- 0
write.table(table_d, file="${analysis}/simulation/simulated.numbers.1Mio.${TYPE}.10groups.table", quote=F, col.names=T, row.names = F)
### step 2: summarising the tables --- END ---


### step 3: measuring the linkage --- START ---
for TYPE in 0percent 1percent 3percent 5percent 10percent; do
awk -v TYPE=${TYPE} '{if(NR>1) {count1+=$2*$3; count2+=$4*$5; count3+=$6*$7; count4+=$8*$9; count5+=$10*$11}
} END {print count1,TYPE"\n"count2,TYPE"\n"count3,TYPE"\n"count4,TYPE"\n"count5,TYPE
}' ${analysis}/simulation/simulated.numbers.1Mio.${TYPE}.10groups.table >> ${analysis}/simulation/simulated.numbers.1Mio.stats
done

### in-trans-pingpong_plots.R was used to make a barchart in Extended Data Figure 9.
