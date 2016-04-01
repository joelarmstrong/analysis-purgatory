# Plot coverage vs. distance to the nearest outgroup.
library(ggplot2)

laurasiatheria = read.table("laurasiatheriaCoverage", sep="\t", header=T)
primates = read.table("primateCoverage", sep="\t", header=T)
atlantogenata = read.table("atlantogenataCoverage", sep="\t", header=T)
crocobirds = read.table("crocobirdsCoverage", sep="\t", header=T)
glires = read.table("gliresCoverage", sep="\t", header=T)
fish = read.table("fishCoverage", sep="\t", header=T)

laurasiatheria$Clade = rep("laurasiatheria", length(laurasiatheria$ingroup))
primates$Clade = rep("primates", length(primates$ingroup))
atlantogenata$Clade = rep("atlantogenata", length(atlantogenata$ingroup))
crocobirds$Clade = rep("crocobirds", length(crocobirds$ingroup))
glires$Clade = rep("glires", length(glires$ingroup))
fish$Clade = rep("fish", length(fish$ingroup))


df = rbind(laurasiatheria, primates, atlantogenata, crocobirds, fish, glires)

pdf("coveragePlot.pdf")
print(ggplot(df, aes(x=firstOutgroupDist, y=coverage, shape=Clade, label=ingroup, color=numOutgroupsTotal)) + geom_point() + theme_bw() + ylab("coverage (cumulative)") + scale_color_gradient(low="red", high="black")) # + geom_text(size=rel(0.95), hjust=-0.15, vjust=0)
dev.off()
