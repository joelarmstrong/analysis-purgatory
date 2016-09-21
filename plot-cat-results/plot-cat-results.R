#!/usr/bin/env Rscript
# Plot Comparative-Annotation-Toolkit results
if (!require("pacman")) install.packages("pacman", repos = "https://cloud.r-project.org/")
pacman::p_load(argparse, ggplot2, tidyr, RSQLite, DBI, dplyr)

parser <- ArgumentParser(description = 'Plot statistics from a ComparativeAnnotationToolkit db')
parser$add_argument('database', nargs = '+', help = 'path to database')
parser$add_argument('--output', help = 'Output plot file prefix')
parser$add_argument('--labels', nargs = '+', help = 'Labels for each input file')
args <- parser$parse_args()
dfs <- sapply(args$database, function (x) {
    con = dbConnect(RSQLite::SQLite(), dbname=x)
    dbGetQuery(con, 'select * from mRNA_transMap_Metrics')
}, simplify=F)

df <- bind_rows(dfs, .id = 'Run')

# Apply the labels, if any
if (!is.null(args$labels)) {
    file_to_label <- data.frame(filename=args$database, label=args$labels)
    matching_rows <- match(df$Run, file_to_label$filename)
    df$Run <- file_to_label[matching_rows, ]$label
}

df <- spread(df, classifier, value)

df$AlnCoverage <- as.numeric(df$AlnCoverage)
df$AlnIdentity <- as.numeric(df$AlnIdentity)
df$Badness <- as.numeric(df$Badness)
df$NumMissingIntrons <- as.numeric(df$NumMissingIntrons)
df$NumMissingExons <- as.numeric(df$NumMissingExons)
df$PercentUnknownBases <- as.numeric(df$PercentUnknownBases)
df$RefTranscriptId <- sapply(df$TranscriptId, function(x) { split = strsplit(x, "-")[[1]]; return(paste(split[1:(length(split)-1)], sep="-"))})

ref_df <- df %>% group_by(Run, RefTranscriptId) %>% summarize(maxCov=max(AlnCoverage), paralogy=n())

glimpse(df)

ref_df$maxCovBin <- cut(ref_df$maxCov, c(0,.75,.90,.95,.98,.99999999999999,1.00), c("0-75", "75-90", "90-95", "95-98", "95-100", "100"))

## p <- ggplot(df, aes(x=AlnCoverage)) + theme_classic()
## p <- p + geom_histogram(breaks=c(0,50,70,90,95,99,100), position='identity')
p <- ggplot(ref_df, aes(x=paralogy, fill=Run)) + theme_classic()
p <- p + geom_bar(position="dodge") + xlab("Number of copies")

ggsave(paste(args$output, 'paralogy', 'pdf', sep='.'), p)

p <- ggplot(ref_df, aes(x=Run, fill=maxCovBin)) + theme_classic() + theme(axis.text.x=element_text(angle=45, hjust=1))
p <- p + geom_bar(position="stack")

ggsave(paste(args$output, 'coverage.stacked', 'pdf', sep='.'), p)


p <- ggplot(subset(ref_df, maxCov >= .9), aes(x=maxCov, fill=Run)) + theme_classic()
p <- p + geom_histogram(position="dodge")

ggsave(paste(args$output, 'coverage', 'pdf', sep='.'), p, width=8*16/9, height=8)
