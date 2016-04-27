#!/usr/bin/env Rscript
if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org/")
pacman::p_load(argparse, ggplot2)

parser <- ArgumentParser(description = 'Plot a look at where coverage has been lost.')
parser$add_argument('coverageLostAtEachStep')
parser$add_argument('--output', default='output.pdf')
args <- parser$parse_args()

df <- read.table(args$coverageLostAtEachStep, header=T)
p <- ggplot(df, aes(x=ParentGenome, y=CoverageLost))
p <- p + scale_x_discrete(limits=df$ParentGenome)
p <- p + geom_bar(stat="identity")
p <- p + theme_classic()
p <- p + theme(axis.text.x=element_text(angle=45, hjust=1))
ggsave(args$output, p)
