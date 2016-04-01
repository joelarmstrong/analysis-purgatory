# Analyse mafCoverage files produced from the different 100ways and
# produce comparison plots. The files should have been created by the
# mafCoverage tool from mafTools.
require(ggplot2)
require(reshape)

cladeMap <- read.table("cladeMap", sep="\t", header=T)
genomeToCommonName <- read.table("genomeToCommonName", sep="\t", header=T)

bakeOff <- read.table("mafCoverage-bakeOff", sep="\t", header=T)
multiz <- read.table("mafCoverage-multiz", sep="\t", header=T)
old100way <- read.table("mafCoverage-old100way", sep="\t", header=T)
tenWay <- read.table("mafCoverage-tenWay", sep="\t", header=T)

old100way$Alignment <- rep("Cactus 100-way (Old)",length(old100way$"querySpecies.Chr"))
bakeOff$Alignment <- rep("Cactus 100-way (New)",length(bakeOff$"querySpecies.Chr"))
multiz$Alignment <- rep("MultiZ",length(multiz$"querySpecies.Chr"))
tenWay$Alignment <- rep("Cactus 10-way",length(tenWay$"querySpecies.Chr"))

df <- rbind.fill(old100way, bakeOff, multiz)
# Probably hacky
df$Clade = factor(unlist(sapply(df$"querySpecies.Chr", function(x) { 


         name <- cladeMap[cladeMap$Genome == as.character(x),]$Clade
         if(length(name) == 0) {
                         name = "NA"
         }
         return(as.character(name))


 })))

# Get rid of stuff that isn't shared between all 100-ways
df = subset(df, Clade != "NA")

pdf("mafCoverage-100ways.pdf", width=7.5, height=5)
print(ggplot(df, aes(x=Alignment, y=coverage, group=querySpecies.Chr, color=Clade)) + geom_line() + theme_classic() + scale_x_discrete(limits=c("Cactus 100-way (Old)", "Cactus 100-way (New)", "MultiZ")) + scale_color_brewer(type="qual", palette=6) + ylab("Coverage on hg19"))
dev.off()

### Now the plot comparing the 10-way to the 100-way and multiZ,
### showcasing that the bad assemblies are dragging us down.

df <- rbind.fill(bakeOff, multiz, tenWay)
# Probably hacky
df$Clade = factor(unlist(sapply(df$"querySpecies.Chr", function(x) { 


         name <- cladeMap[cladeMap$Genome == as.character(x),]$Clade
         if(length(name) == 0) {
                         name = "NA"
         }
         return(as.character(name))


 })))

# Get rid of stuff that isn't in the 10-way
df = subset(df, Clade != "NA")
df = subset(df, querySpecies.Chr %in% tenWay$querySpecies.Chr)

# Fill in the common names.
df$CommonName = factor(unlist(sapply(df$querySpecies.Chr, function(x) {return(genomeToCommonName[genomeToCommonName$Genome == as.character(x),]$CommonName)})), levels=c("Rhesus", "Horse", "Cat", "Dog", "Cow", "Pig", "Mouse", "Rat", "Human"))

pdf("mafCoverage-10way-vs-100ways.pdf", width=7.5, height=5)
print(ggplot(df, aes(x=Alignment, y=coverage, group=querySpecies.Chr, color=CommonName)) + geom_line() + theme_classic() + scale_color_brewer(type="qual", palette = 3) + ylab("Coverage on hg19"))
dev.off()

pdf("mafCoverage-10way-vs-100ways-bar.pdf", width=7.5, height=5)
print(ggplot(df, aes(x=CommonName, y=coverage, fill=Alignment)) + geom_bar(stat="identity", position="dodge") + theme_classic() + ylab("Coverage on hg19"))
dev.off()
