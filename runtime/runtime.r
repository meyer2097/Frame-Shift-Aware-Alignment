library(reshape2)

d = read.csv("median.csv", header = TRUE)
meltedData = melt(d, id.vars = "X")
names(meltedData) = c("DNA", "AA", "Time")

meltedData["AA"] = lapply(meltedData["AA"], function(x) substring(x, 2))
meltedData$Time = meltedData$Time*1000
meltedData$AA = as.numeric(meltedData$AA)
meltedData$DNA = as.numeric(meltedData$DNA)

l = lm(Time ~ AA:DNA, data = meltedData)
l
