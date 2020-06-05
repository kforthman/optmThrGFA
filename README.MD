Optimization of threshold parameters for GFA replicates

# Example Script

The example script "run_optmThrGFA.R" was created so that you can get the package optmThrGFA up and running quickly and painlessly! Please keep in mind that you can make edits to the code in this script for your own project as you see fit. Please read through the comments for various tips and suggestions.

## To download the Example files...
Open the terminal and cd to the directory where you would like to save the example files. Then execute the command

`svn export https://github.com/kforthman/optmThrGFA/trunk/Example`

## Before running this script...

1. Create a new folder for your project. Save this script and the .Rmd file to that folder.
2. Format your raw data.
    Your raw data is expected as an R data.frame object named "MYdf". You need to save this object as a .RData file. If your data is in the form of a .csv, you will need to load it into R and format it this way yourself. Remove any columns that are not predictors (i.e. ID or session). Commented below is an example of how you might do that. You only need to do this once, so do not put the following code into your main script.
    
```
data <- read.csv("mydata.csv")
data.df <- as.data.frame(data)
MYdf <- data.df[!(names(data.df) %in% c("ID", "session"))]
save(MYdf, file = "mydata.RData")
```

3. Create a folder in your project directory named 'Data'. Move the .RData file containing your formatted raw data into this folder.
4. Create a csv file containing block assignments for each variable. Name the first column 'block' and the second column 'variable'. Name the file "var_blocks.csv". You can name the blocks however you like, just don't use any special characters or spaces in the block names. In the variable column, make sure the variable names match the column names in your raw dataframe. An example is in the repository in the Example directory. Place this file in you project directory.
5. Fill in the blanks in the example script marked by '__________'
6. Save the example bash scripts named 'run_R.Rscript', 'run_run_R.ssub', 'xvfb-run-safe' to your project folder. Fill in the blanks in the bash script 'run_run_R.ssub'.

## Running this script...

It is highly recommended that you run this script on submit0. It will run much faster and not take up resources on your local computer. The files 'run_R.Rscript', 'run_run_R.ssub', and 'xvfb-run-safe' are provided in this example for running this script on the server. To run the job,

1. Open the terminal.
2. Start submit0 by executing the command

    `ssh submit0`
    
    Type in your password.
3. cd into the directory where this script is located.

    `cd /your/project/directory`
    
4. Run the job.

    `sbatch < run_run_R.ssub`
    
Pro tip: you can monitor progress in the terminal using the command

`watch tail out.log -n 50`
