library(cluster)
library(Biobase)
library(qvalue)
library(fastcluster)
library(factoextra)
library(ggplot2)
library(reshape2)
options(stringsAsFactors = FALSE)
NO_REUSE = F

# This script is based on the scripts generated use https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Differential-Expression
expression.tab <- "/PATH/TO/tab-delimited_expression.matrix"
samples.txt <- "/PATH/TO/tab-delimited_samples.txt"
plot.out <- "/PATH/TO/output-folder"

data <- read.delim(expression.tab, row.names=1)

myheatcol = colorpanel(75, 'blue','white','red')
samples_data = read.table(samples.txt, header=F, check.names=F, fill=T)
samples_data = samples_data[samples_data[,2] != '',]
colnames(samples_data) = c('sample_name', 'replicate_name')
sample_types = as.character(unique(samples_data[,1]))
rep_names = as.character(samples_data[,2])
data = data[, colnames(data) %in% rep_names, drop=F ]
nsamples = length(sample_types)
sample_colors = rainbow(nsamples)
names(sample_colors) = sample_types
sample_type_list = list()
for (i in 1:nsamples) {
    samples_want = samples_data[samples_data[,1]==sample_types[i], 2]
    sample_type_list[[sample_types[i]]] = as.vector(samples_want)
}
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
initial_matrix = data # store before doing various data transformations
data = log2(data+1)
sample_factoring = colnames(data)
for (i in 1:nsamples) {
    sample_type = sample_types[i]
    replicates_want = sample_type_list[[sample_type]]
    sample_factoring[ colnames(data) %in% replicates_want ] = sample_type
}
sampleAnnotations = matrix(ncol=ncol(data),nrow=nsamples)
for (i in 1:nsamples) {
  sampleAnnotations[i,] = colnames(data) %in% sample_type_list[[sample_types[i]]]
}
sampleAnnotations = apply(sampleAnnotations, 1:2, function(x) as.logical(x))
sampleAnnotations = sample_matrix_to_color_assignments(sampleAnnotations, col=sample_colors)
rownames(sampleAnnotations) = as.vector(sample_types)
colnames(sampleAnnotations) = colnames(data)
data = as.matrix(data) # convert to matrix

# Centering rows
data = t(scale(t(data), scale=F))

if (nrow(data) < 2) { stop(

**** Sorry, at least two rows are required for this matrix.

);}
if (ncol(data) < 2) { stop(

**** Sorry, at least two columns are required for this matrix.

);}

if (nrow(data) <= 1) { message('Too few genes to generate heatmap'); quit(status=0); }

# Setting up the random seed.
set.seed(42)
# Define Kmeans algorithm.
MyKmeansFUN <- function(x, k) kmeans(x, k, algorithm=MacQueen, iter.max=200)

# Plots to determine K.
png(file.path(plot.out, "elbow.png"))
fviz_nbclust(data, FUNcluster=MyKmeansFUN, method = wss, k.max=20)+
  labs(subtitle = Elbow method)
dev.off()
png(file.path(plot.out, "silhouette.png"))
fviz_nbclust(data, FUNcluster=MyKmeansFUN, method = silhouette, k.max=20)+
  labs(subtitle = Silhouette method)
dev.off()
png(file.path(plot.out, "gap.png"))
fviz_nbclust(data, FUNcluster=MyKmeansFUN, method = gap_stat, k.max=20, nboot=500)+
  labs(subtitle = Gap statistics method)
dev.off()