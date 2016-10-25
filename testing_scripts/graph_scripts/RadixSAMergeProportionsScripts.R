setwd("~/Desktop")
data = read.csv("radix_merge_proportions.csv")
data
bsadata =factor(data$func)
plot(data$data, data$prop, type="line")

f= split(data, data$func)

plot(f$merge$data, f$constructTotalRadixSA$prop, type="l", ylim=c(0,1.4),
xlab="Input data (thousands)", ylab="Proportion of total execution",
main="Proportion of time spent executing sub-processes of RadixSAMergeConstruction")
lines(f$merge$data, f$bsatimer$prop, type="l", col="green")
lines(f$merge$data, f$radixtimer$prop, type="l", col="purple")
lines(f$merge$data, f$maptimer$prop, type="l", col="blue")
lines(f$merge$data, f$merge$prop, type="l", col="red")
legend(locator(1),c("BSA","RadixSA", "Map", "Merge"),pch=c(1,16), col=c("green","purple", "blue", "red"))

