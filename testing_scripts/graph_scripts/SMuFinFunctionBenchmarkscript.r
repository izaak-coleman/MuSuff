data <- read.csv("~/Documents/Computing/MscComputingsci/thesis/src/SMuFin/timestamps/std_clock/suntimes.csv")
data

plot(data$function., data$time.s., las=2, cex.axis=0.7, main="Total execution time of function calls in SMuFin", xlab = "Function names", ylab="Time (s)")
par(xaxt="n")  % remove x axis plot

