
data  = read.csv("ICSMuFinRuntime.csv")
data
data$ms/data$ms$main

class_split = split(data, data$class)
class_split
par(mfrow=c(1,3))
plot(class_split$'main()'$dats, class_split$'main()'$ms/1000/60,  main="IC-SMuFin runtime",
ylab="Time (min)", xlab = "Input data (thousands)")
plot(class_split$'bpg()'$dats, class_split$'bpg()'$ms/1000/60, col="blue",
main="Runtime of IC-SMuFin classes", xlab="Input data (thousands)", ylab="Time (min)")
points(class_split$'suf_array()'$dats, class_split$'suf_array()'$ms/1000/60,  col="red")
points(class_split$'reads_manip()'$dats, class_split$'reads_manip()'$ms/1000/60,  col="green")
points(class_split$'gm()'$dats, class_split$'gm()'$ms/1000/60, col="purple")
legend(locator(1),c("bpg","suf", "reads", "gm"),pch=c(1,16), col=c("blue","red", "green", "purple"))

plot(class_split$'main()'$dats, class_split$'main()'$prop, type="l", col="black", ylim=c(0,1.4),
main="Class runtime : total runtime ratio", xlab="Input data (thousands)", ylab="Proportion")
lines(class_split$'bpg()'$dats, class_split$'bpg()'$prop, type="l", col="blue")
lines(class_split$'suf_array()'$dats, class_split$'suf_array()'$prop, type="l", col="red")
lines(class_split$'reads_manip()'$dats, class_split$'reads_manip()'$prop, type="l", col="green")
lines(class_split$'gm()'$dats, class_split$'gm()'$prop, type="l", col="purple")
 legend(locator(1),c("bpg","suf", "reads", "gm"),pch=c(1,16), col=c("blue","red", "green", "purple"))

length(class_split$'bpg()'$dats)
class_split$'bpg()'[3]
as(class_split$'gm()'$dats)
length(as.vector(x))