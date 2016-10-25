setwd("~/Documents/Computing/MscComputingsci/thesis/")
data = read.csv("memusage.csv")
data
eq25 = function(x){x*25}
eq26 = function(x){x*26}
per_dat = split(data, data$datsize)
per_class = split(data, data$class)
per_dat$`32`$class
per_class
gbin = per_class$` main`$kb/2^20
gbout = per_class$` main`$input_kb/2^20

per_dat$`70554`
par(mfrow=c(3,3))
barplot(per_dat$`128`$kb/2^20, names.arg = per_dat$`32`$class, ylab="RSS (GB) after class constructor",
xlab="Order of class instantiation", main="128k reads (0.009 GB)")
barplot(per_dat$`256`$kb/2^20, names.arg = per_dat$`32`$class, ylab="RSS (GB) after class constructor",
xlab="Order of class instantiation", main="256k reads (0.018 GB)")
barplot(per_dat$`512`$kb/2^20, names.arg = per_dat$`32`$class, ylab="RSS (GB) after class constructor",
xlab="Order of class instantiation", main="512k reads (0.037 GB)")
barplot(per_dat$`1024`$kb/2^20, names.arg = per_dat$`32`$class, ylab="RSS (GB) after class constructor",
xlab="Order of class instantiation", main="1024k reads (0.074 GB)")
barplot(per_dat$`2048`$kb/2^20, names.arg = per_dat$`32`$class, ylab="RSS (GB) after class constructor",
xlab="Order of class instantiation", main = "2048k reads (0.148 GB)")
barplot(per_dat$`4096`$kb/2^20, names.arg = per_dat$`32`$class, ylab="RSS (GB) after class constructor",
xlab="Order of class instantiation", main="4096k reads (0.296 GB)")
barplot(per_dat$`70554`$kb/2^20, names.arg = per_dat$`70554`$class, ylab="RSS (GB) after class constructor",
xlab="Order of class instantiation", main="70554k reads (1.283 GB)")
barplot(per_dat$`104522`$kb/2^20, names.arg = per_dat$`32`$class, ylab="RSS (GB) after class constructor",
xlab="Order of class instantiation", main="In silico data set (1.890 GB)")
plot(per_class$` main`$input_kb/2^20, per_class$` main`$kb/2^20, col="red", ylab="Peak memory usage (GB)",
xlab="Input data (GB)", main="Peak memory per input size")
abline(lm(gbin ~ gbout), col="red")
lines(per_class$` main`$input_kb/2^20,eq25(per_class$` main`$input_kb/2^20), col="green")
lines(per_class$` main`$input_kb/2^20,eq26(per_class$` main`$input_kb/2^20), col="blue")
 legend(locator(1),c("Regression","f(x) = x*25", "f(x) = x*26"),pch=c(1,16), col=c("red","green", "blue"))





