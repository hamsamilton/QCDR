#!/usr/bin/env Rscript

# install / load packages 
packages <- c("optparse", "tidyverse","rio")

lapply(packages,function(package){
    if(!require(package,character.only = T)){
        install.packages(package)}
    require(package, character.only = T)})

option_list = list(
    make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name, a counts table"),
    make_option(c("-o", "--out"), type="character", default="out.csv", 
              help="output file name"),
    make_option(c("-w", "--width"), type="numeric", default=.25, 
              help="the width of each bin (log2(CPM))"),
    make_option(c("-m", "--max"), type="numeric", default=18.5, 
              help="The maximum CPM"),
    make_option(c("-r", "--rownames"),type ="numeric",default = 2,
              help="are there preceding rows specifying gene names? If so, how many?[default= %default]",metavar="numeric")
)
 
opt_parser = OptionParser(option_list=option_list)
opt        = parse_args(opt_parser)

if (is.null(opt$file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# input count df 
df = rio::import(opt$file)

#remove rownames
df = df[,-c(1:opt$rownames)]

create_bin_hist <- function(df, binsize = .25,maxdepth = 18.5){

# transform to matrix
df = as.matrix(df)
class(df) = "numeric"
# CPM transform matrix

readDepth = colSums(df)
readDepth = readDepth / 1e6
df = sweep(df, 2, readDepth, FUN = '/')

# log2 transform data
df <- log2(df + 1)

# make bins
bins = seq(0,maxdepth - binsize,binsize)

bin_counts = sapply(bins,function(bin){

    # get the top of the bin 
    topbin = bin + binsize

    # count how many genes are within the top and bottom of the bin
    countdf = (df > bin & df <= topbin)
    countdf2 = colSums(countdf)

    return(countdf2)

})

# flip matrix
bin_counts = t(bin_counts) %>% as.data.frame()

# add column for xaxis
topbins = bins + binsize
Xaxis_names = str_c("(",bins,",",topbins,"]")
bin_counts = bin_counts %>% mutate(Xaxis = Xaxis_names, .before = 1)

return(bin_counts)}

hist_out = create_bin_hist(df = df,binsize = opt$width,maxdepth = opt$max)


write.csv(hist_out,opt$out,row.names = F)
