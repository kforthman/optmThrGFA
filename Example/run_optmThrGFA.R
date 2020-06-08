# Before running this script...
# (1) Create a new folder for your project. Save this script and the .Rmd file to that folder.
# (2) Format your raw data.
#     Your raw data is expected as an R data.frame object named "MYdf". You need to save this object as a .RData file. If your data is in the form of a .csv, you will need to load it into R and format it this way yourself. Remove any columns that are not predictors (i.e. ID or session). Commented below is an example of how you might do that. You only need to do this once, so do not put the following code into your main script.
# ```
# data <- read.csv("mydata.csv")
# data.df <- as.data.frame(data)
# MYdf <- data.df[!(names(data.df) %in% c("ID", "session"))]
# save(MYdf, file = "mydata.RData")
# ```
# (3) Create a folder in your project directory named 'Data'. Move the .RData file containing your formatted raw data into this folder.
# (4) Create a csv file containing block assignments for each variable. Name the first column 'block' and the second column 'variable'. Name the file "var_blocks.csv". An example is in the repository in the Example directory. Place this file in you project directory.
# (5) Fill in the blanks in this script marked by '__________'
# (6) Save the example bash scripts named 'run_R.ssub', 'run_run_R.ssub', 'xvfb-run-safe' to your project folder. Fill in the blanks in the bash script 'run_run_R.ssub'. 
#
# Running this script...
# It is highly recommended that you run this script on submit0. It will run much faster and not take up resources on your local computer. The files 'run_R.ssub', 'run_run_R.ssub', and 'xvfb-run-safe' are provided in this example for running this script on the server. To run the job,
# (1) Open the terminal.
# (2) Start submit0 by executing the command
#     `ssh submit0`
#     Type in your password.
# (3) cd into the directory where this script is located.
#     `cd /your/project/directory`
# (4) Run the job.
#     `sbatch < run_run_R.ssub`
# Pro tip: you can monitor progress in the terminal using the command
#     `watch tail out.log -n 50`
#
# This example script was created so that you can get the package optmThrGFA up and running quickly and painlessly! Please keep in mind that you can make edits to the code in this script for your own project as you see fit. Please read through the comments for various tips and suggestions.

bgn_time <- Sys.time()
message("Begin script: run_optmThrGFA.R")
message(bgn_time)

# ---- Set working directory ----
# This sets the working directory to your local project path
# ONLY if it is being run on your computer. That way, you do 
# not have to make any changes when you are running this code
# on your computer vs on the server.
local.nodename <- "__________" # to figure out your computers nodename,
# type `Sys.info()["nodename"]` into the console. The returned string
# is the nodename of your computer.
local.ed <- "__________" # The path to the working directory on your computer.

# If you run the script locally, it will set the working directory.
if(Sys.info()["nodename"] == local.nodename){setwd(local.ed)}


## Load optmThrGFA
# The packages 'devtools' and 'rmarkdown' are available on the submit0 shared R library.
# There is no need to install these packages in your personal library on submit0.
library(devtools)
library(rmarkdown)
devtools::install_github("kforthman/optmThrGFA") # Do not comment out this line. This code updates frequently, you should check for updates with every run. If there are no updates, the package will NOT re-install, so this will not slow you down.
library(optmThrGFA)

# ---- Set up parallel computing ----
# Set how many cores you would like to use.
# If you would like to know how many cores are available,
# use the function `detectCores()`. I would recommend using
# 5-10 cores, if available.
cl = __________; registerDoParallel(cl)
getDoParWorkers()


# ---- Set GFA Options ----
opts <- getDefaultOpts()
opts$convergenceCheck <- T
# opts$verbose <- 0 # Set to 0 if you do not 
# want to print results of each gfa replicate


# ---- Set file paths ----
out_dir <- "Output"
if(!file.exists(out_dir)){dir.create(out_dir)} # Will create the output folder if it does not exist
# If you would like to name the results and replicates folders differently,
# go ahead! You can also run another analysis without over
# -writing previous results by changing the folder names here. Note that 
# you do not have to create this folder before running
# the script. R will look to see if the folder exists, and if it doesn't
# exist the folder will be created.
folder_res_name <- "res" # The folder where you want to put the results
folder_rep_name <- "rep" # The folder where you want to put the GFA replicates
folder_res <- paste0(out_dir,"/",folder_res_name)
folder_rep <- paste0(out_dir,"/",folder_rep_name)
if(!file.exists(folder_res)){dir.create(folder_res)} # Will create the results folder if it does not exist
if(!file.exists(folder_rep)){dir.create(folder_rep)} # Will create the results folder if it does not exist

data_dir <- "Data" # The folder that contains your raw data

# ---- Load raw data ----
# Raw data should be a matrix or data frame with a row 
# for each observation/participant,
# and a column for each variable
message("\nLoading original dataset.")
raw_data_filename <- "__________.RData"
load(paste0(data_dir, "/", raw_data_filename))
message("Original dataset loaded.\n")

# Impute missing data.
Y <- as.matrix(MYdf)
if(sum(is.na(Y)) > 0){
  Y <- knnImputation(Y)
}

# ---- Setting up the blocks ----  

Y_grouped <- list()

var_blocks <- read.csv(paste0(data_dir, "/", "var_blocks.csv"), header = T)
block.names <- unique(var_blocks$block)

for (k in block.names){
  this_block_vars <- var_blocks[var_blocks$block == k, "variable"]
  Y_grouped[[k]] <- Y[,this_block_vars]
}


# Normalize the data
mynorm <- normalizeData(Y_grouped, type="scaleFeatures")
Y_norm <- mynorm$train # pull out the normalized data
Y_norm_bound <- do.call(cbind, Y_norm) # bind the groups into matrix


varIdx.by.block <- list()
for (k in 1:length(block.names)){
  this_block_vars <- var_blocks[var_blocks$block == block.names[k], "variable"]
  this_block_idx <- which(is.element(colnames(Y_norm_bound), this_block_vars))
  varIdx.by.block[[k]] <- this_block_idx
}

block.length <- unlist(lapply(varIdx.by.block, length))
n.blocks <- length(block.length)

block4vars <-unlist(
  lapply(1:n.blocks, function(x){
    rep(block.names[x], block.length[x])
  }))

allvarlabs <- var_blocks$variable

save(block.names, file = paste0(folder_res, "/block.names.rda"))
save(varIdx.by.block, file = paste0(folder_res, "/varIdx.by.block.rda"))
save(block4vars, file = paste0(folder_res, "/block4vars.rda"))

# Set the default K (the initial guess to the number of factors) 
# to the number of variables in the dataset.
startK <- dim(Y)[2]

# ---- Overwrite  settings ---
# Indicate if you would like to overwrite files if they already exist.
# This is useful for testing, if you made a mistake and need to
# re-run a step in the process.
overwrite_rep <- T
overwrite_gfaList <- T
overwrite_xw <- T
overwrite_match <- T

# ---- Create and load replicates ---- 

# Number of GFA replicates
# If you don't know, just set this to 10.
R <- __________

# Writes output from gfa function to a .rda file.
foreach (r = 1:R) %dopar% {
  this.filename <- paste0(folder_rep, "/GFA_rep_", r, ".rda")
  if(!file.exists(this.filename) | overwrite_rep){
    message(paste0("Creating replicate ", r, " of ", R))
    set.seed(r)
    res <- gfa(Y_norm, opts=opts, K=startK)
    save(res, file = this.filename)
    remove(res)
  }
}

# Load in the generated GFA results. Recorded in the data.frame `rep.summ` are the values
# - conv: An estimate of the convergence of the model's reconstruction based on Geweke diagnostic.
# Values significantly above 0.05 imply a non-converged model,
# and hence the need for a longer sampling chain.
# - K: The number of components inferred.
rep.summ <- data.frame(Replicate=1:R, conv=rep(NA, R), K=rep(NA, R))
gfaList_full <- list()
for(r in 1:R){
  message(paste0("Loading replicate ", r, " of ", R))
  load(paste0(folder_rep, "/GFA_rep_", r, ".rda"))
  gfaList_full[[r]] <- res
  rep.summ$conv[r] <- res$conv
  rep.summ$K[r] <- res$K
  message(paste0("Replicate ", r, " has ", rep.summ$K[r], " factors."))
}
remove(res)

message("\nAll replicates loaded.\n")


# ---- Analysis ---- 

# gfaList
gfaList.filename <- paste0(folder_res, "/gfaList_p50.rda")
if(!file.exists(gfaList.filename) | overwrite_gfaList){
  message("\nCreating GFA list.")
  gfaList_p50 <- list()
  for (r in 1:R){ 
    message(paste0("running replicate ", r, " of ", R))
    gfaList_p50[[r]] <- psSummary(gfa.res=gfaList_full[[r]], credible.lv=0.95)
  }
  save(gfaList_p50, file = gfaList.filename)
  message("GFA list created.")
}else{
  message("\nLoading GFA list.")
  load(gfaList.filename)
  message("GFA list loaded.")
}

# xw
xw_by_rep_comp.filename <- paste0(folder_res, "/xw_by_rep_comp.rda")
if(!file.exists(xw_by_rep_comp.filename) | overwrite_xw){
  message("\nCreating xw_by_rep_comp.")
  xw <- list()
  for(r in 1:R){
    message(paste(c("Loading|", rep("-", r), rep(" ", R-r), "|"), collapse = ""))
    load(paste0(folder_rep, "/GFA_rep_", r, ".rda"))
    xw[[r]] <- pmXW_by_factor(res)
  }
  
  save(xw, file = xw_by_rep_comp.filename)
  message("xw_by_rep_comp created.")
}else{
  message("\nLoading xw_by_rep_comp.")
  load(xw_by_rep_comp.filename)
  message("xw_by_rep_comp loaded.")
}

# corGrids and matchGrids are the threshold parameters.
# corThr: How close two components are required to be, 
# in terms of correlation, in order to match them.
corGrids   <- c(seq(0.1, 0.5, by=0.2), seq(0.6, 0.9, 0.1))
# matchThr: The proportion of sampling chains that need to contain a component 
# in order to include it in the robust components.
matchGrids <- c(seq(0.1, 0.5, by=0.2), seq(0.6, 0.9, 0.1))

match.xw.filename <- paste0(folder_res, "/match.xw_pm.rda")
if(!file.exists(match.xw.filename) | overwrite_xw){
  message("\nCreating match.xw")
  
  match.xw <- match_dopar(comps=xw, corGrids, matchGrids)
  
  save(match.xw, file = match.xw.filename)
  message("match.xw created.")
}else{
  message("\nLoading match.xw.")
  load(match.xw.filename)
  message("match.xw loaded.")
}

match.filename <- paste0(folder_res, "/match_orig.rda")
if(!file.exists(match.filename) || overwrite_match){
  # 3. Run MSE.Grids
  match.mse <- MSE.Grids(
    Ymtx=Y_norm_bound,            ## the observed (normalized) data matrix (N x D)
    comps=xw,                     ## a list of GFA replicates with posterior medians
    match.res = match.xw,         ## use the matching results from match_dopar()
    corGrids=corGrids,            ## the grids of corThr values to be assessed
    matchGrids=matchGrids)        ## the grids of matchThr values to be assessed
  save(match.mse, file = match.filename)
  
}else{
  load(match.filename)
}

# ---- Summary data ---- 
# These data frames will be used in the Rmarkdown script to display
# the results.
tmp.filename <- paste0(folder_res, "/opt.par.rda")
message(tmp.filename)
opt.par <- optimizeK(K.grids=match.mse$K.grid, mse.array=match.mse$mse$all)
save(opt.par, file = tmp.filename)

last.idx <- dim(opt.par$par.1se)[1]
opt.cor <- opt.par$par.1se[last.idx, "opt.corThr"]
opt.match <- opt.par$par.1se[last.idx, "opt.matchThr"]
opt.K <- opt.par$par.1se[last.idx, "optK"]
opt.cor.idx <- which(corGrids == opt.cor)
opt.match.idx <- which(matchGrids == opt.match)

Kgrids.filename <- paste0(folder_res, "/Kgrids.rda")
message(Kgrids.filename)
Kgrids <- match.mse$K.grids
save(Kgrids, file = Kgrids.filename)

optParams.filename <- paste0(folder_res, "/optParams.rda")
message(optParams.filename)
optParams <- Kgrids
optParams[opt.par$mse.m > opt.par$mseThr | Kgrids != opt.par$Krobust.1se] <- NA
save(optParams, file = optParams.filename)

varexp.filename <- paste0(folder_res, "/varexp.rda")
message(varexp.filename)

rob.ind <- list(indices=match.xw[[opt.cor.idx]][[opt.match.idx]]$indices)
varexp <- rob.var.exp(models=gfaList_p50,
                      indices=rob.ind,
                      block.names=block.names,
                      varIdx.by.block=varIdx.by.block,
                      use.unmatched=T,
                      by.block=T)
save(varexp, file = varexp.filename)

robust.xw.filename <- paste0(folder_res, "/robust.xw.rda")
message(robust.xw.filename)
robust.xw <- rob_wx(models=gfaList_p50, indices=varexp$indices, block.labs=block4vars, var.labs=allvarlabs)
save(robust.xw, file=robust.xw.filename)

# ---- Robust Components ----
rcomp.filename <- paste0(folder_res, "/rcomp.rda")
message(rcomp.filename)
if(!file.exists(rcomp.filename)){
  rcomp <- robustComponents(gfaList_full, 
                            corThr=opt.par$par.1se[1, "opt.corThr"], 
                            matchThr=opt.par$par.1se[1, "opt.matchThr"]
  )
  save(rcomp, file = rcomp.filename)
}else{
  load(rcomp.filename)
}

# ---- Create the summary report ----
# Global environment is passed to the Rmarkdown script.
message("Creating report.")
rmarkdown::render("createReport.Rmd", output_file = "Report.pdf", clean = T)

message("Done.")
end_time <- Sys.time()
message(end_time)
message(paste0("Script finished in ", round((end_time - bgn_time)/60, 3), " minutes."))