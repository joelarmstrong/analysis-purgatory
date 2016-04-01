# Plot the size of each ancestor relative to the average size of its
# children, relative to the divergence of its children.
library(ggplot2)

laurasiatheria = read.table("laurasiatheriaSeqLossRatio", sep="\t", header=T)
primates = read.table("primateSeqLossRatio", sep="\t", header=T)
atlantogenata = read.table("atlantogenataSeqLossRatio", sep="\t", header=T)
crocobirds = read.table("crocobirdSeqLossRatio", sep="\t", header=T)
glires = read.table("gliresSeqLossRatio", sep="\t", header=T)
fish = read.table("fishSeqLossRatio", sep="\t", header=T)

laurasiatheria$Clade = rep("laurasiatheria", length(laurasiatheria$name))
primates$Clade = rep("primates", length(primates$name))
atlantogenata$Clade = rep("atlantogenata", length(atlantogenata$name))
crocobirds$Clade = rep("crocobirds", length(crocobirds$name))
glires$Clade = rep("glires", length(glires$name))
fish$Clade = rep("fish", length(fish$name))

df = rbind(laurasiatheria, primates, atlantogenata, crocobirds, glires, fish)

pdf("seqLossPlot.pdf")
print(ggplot(df, aes(x=distBtwnChildren, y=sequenceLossRatio, color=Clade)) + geom_point() + theme_bw())
dev.off()

# hundredWay vs hundredWayv2
hundredWayOld = read.table("100way-oldSeqLossRatio", sep="\t", header=T)
hundredWayForBakeOff = read.table("100way-forBakeOffSeqLossRatio", sep="\t", header=T)
hundredWayOld$Run = rep("Old", length(hundredWayOld$name))
hundredWayForBakeOff$Run = rep("multipleOutgroups & realisticRef", length(hundredWayForBakeOff$name))
hundredWays = rbind(hundredWayOld, hundredWayForBakeOff)
pdf("seqLossPlot-2-hundredWays.pdf", width=7, height=6)
print(ggplot(hundredWays, aes(x=distBtwnChildren, y=sequenceLossRatio, color=Run)) + geom_point() + theme_classic() + geom_smooth(method="lm"))
dev.off()
