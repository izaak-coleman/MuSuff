# raw data for reference

#cns_position, occurences
#0, 43533
#1, 14776
#2, 13177
#3, 13103
#4, 12714
#5, 12524
#6, 12213
#7, 12167
#8, 11741
#9, 11582
#10, 11192
#11, 11009
#12, 10721
#13, 10538
#14, 10325
#15, 9863
#16, 9745
#17, 9492
#18, 9118
#19, 9091
#20, 8799
#21, 8540
#22, 8597
#23, 8504
#24, 8189
#25, 8122
#26, 7810
#27, 7693
#28, 7333
#29, 7439
#30, 7177
#31, 7113
#32, 6835
#33, 6640
#34, 6435
#35, 6291
#36, 6097
#37, 5957
#38, 5753
#39, 5449



 setwd("~/Documents/Computing/MscComputingsci/thesis/writeups/msc-template/figures/")
 data = read.csv("cns_bias.csv")
 data
 
 
 
 plot(data$cns_position, data$occurences, type = "h", col="red", xlab="Position from either end of consensus sequence", ylab="Frequency of identified SNVs", main="Bias relationship between consensus sequence position and SNV identification")
 
 barplot(data$occurences, names.arg = data$cns_position)
 
 
 
 
 