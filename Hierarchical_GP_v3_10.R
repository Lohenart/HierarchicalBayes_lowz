
# Hierarchical Bayesian model
# Author: Arturo Avelino

# Inspired from Gelman, Andrew (2013-11-01). Bayesian Data Analysis, 
# Third Edition (Chapman & Hall/CRC Texts in Statistical Science) (Page 595). 
# Chapman and Hall/CRC. Kindle Edition. 

# DESCRIPTION

# The code will output 4 files:

# Template_phase_mu_stdError_FromR.dat:
# Template_phase_mu_stdError_FromR_Norma.dat:
# Template_phase_mu_tau_FromR.dat:
# Template_phase_mu_tau_FromR_Norma.dat:

# Note that this code doesn't say or know anything about the peculiar velocity uncertainty.

################################################################

#     USER INPUTS

#--  Choosing the directory

# Band <- 'J'
# Band <- 'Y'
# Band <- 'H'
Band <- 'K'

# KindOfData <- 'CfA'
# KindOfData <- 'CSP' 
# KindOfData <- 'Others'
KindOfData <- 'AllSamples'

# Compute a normalized template?:
NormalizedTemp <- TRUE # Options: (TRUE, FALSE)

#-- Redshift cutoff. I've set 3 options (z=0, 0.01, anything else).
z_lowerLimit = 0.0  
# z_lowerLimit = 0.01
# z_lowerLimit =  # <-- For any other value of redshift cutoff


# Indicate the technique ALREADY USED to compute the GP fitting in the '1_AllData_InitialFit' step.
# This information is just to label the '3_Template' folder where I will save the output.
# FlatPrior= Assuming a flat prior at ~ -17 Abs mag, then computing the hyperparameters using all the LCs simultaneously.
# TempPrior_Hyper = Using a template prior (computed from the Moving Windows Average template) for the Gaussian Process fitting, then determining the hyperparameters using all the LCs simultaneously.
Technique = 'FlatPrior' # Options: ('FlatPrior', 'TempPrior_Hyper'). 'FlatPrior' is the option I will be using for the paper.

#-----(Fixed values)----------------------------

z_upperLimit = 0.09 # 0.06 # upper limit cutoff in redshift

dm15LowerLim = 0.78
dm15UpperLim = 1.6

EBVhost_Min <- -0.4
EBVhost_Max <- 0.4
EBV_mw_Max <- 1

# Considering that the data counts every 1/2 phase and that all of them stars at phase =-20, then to cover the range phase = [-20, 50] I have to account for 140 steps. I count a little less (i.e., 138) to allow for some light-curve data without data very close to epoch=50.
# With steps of phase of 0.2, then I consider the first 300 steps
# With steps of phase of 0.2, and interval [-20, 60]  then I consider the first 400 steps
# With steps of phase of 0.5, and interval [-35, 60]  then I consider the first 190 steps
# With steps of phase of 1, and interval [-35, 60]  then I consider the first 96 steps, to include the zero.
numPhases = 190
numPhases
PhaseIndexZero <- 71 # row when phase = 0

FilterSyst = 'Std_filters'
# FilterSyst = 'CSP_filters'

#   <---- END OF USER INPUTS.
################################################################

#             AUTOMATIC
# The following part of the code does not need user manipulations.

#---------------------------------------------------------------
#- Create the path of the main directory for a given band

MainDir <- file.path('/Users/arturo/Dropbox/Research/Articulos/10_AndyKaisey/10Compute/TheTemplates',paste(Band,'_band', sep=''))
MainDir

#---------------------------------------------------------------

#- Create the 3_template folder and subfolders

# Create the path of the '3_Template' directory:
TemplatePath <- paste('3_Template_', Technique, sep= '')
TemplatePath
dir.create(file.path(MainDir, FilterSyst, TemplatePath))
dir.create(file.path(MainDir, FilterSyst, TemplatePath, KindOfData))

# Defining the subfolder in '3_Template' depending on z_lowerLimit
if(z_lowerLimit == 0.0){GoodPath2 = 'z_gr_0'
} else if (z_lowerLimit == 0.01){GoodPath2 = 'z_gr_001'
} else {GoodPath2 = 'z_gr_'}

DirTemplate1 <- file.path(MainDir, FilterSyst, TemplatePath, KindOfData, GoodPath2)
DirTemplate1
dir.create(DirTemplate1)

# Dir to save the plots
if (NormalizedTemp == TRUE){
  dir.create(file.path(DirTemplate1, 'Plots_histo_Norma'))
} else {
  dir.create(file.path(DirTemplate1, 'Plots_histo'))
}

#---------------------------------------------------------------

#   CREATE THE LIST OF SNe TO BE USED TO CONSTRUCT THE TEMPLATE ('SN_list_template_Norma.txt')

#- Create the path to the GP data:
SelectionPath = paste('2_Selection_', Technique, sep= '') 
GoodPath = 'Goods/'
DirLCData <- file.path(MainDir, FilterSyst, SelectionPath, KindOfData, GoodPath) # OK
DirLCData
# Set as working directory the folder "NIR_band/2_Selection/[sample]/Goods/"
setwd(DirLCData)
getwd() # Show the current directory

#--  Read the data files located in the working directory
# ExtensionName1 <- '*.txt' # old
ExtensionName1 <- paste('*',Band,'.txt', sep='')

list_SNe_GP = list.files(pattern=ExtensionName1, all.files=TRUE) 
# list_SNe_GP
list_SNe_GP[1] # Check that 'list_SNe_GP' has data by showing the first entry.

numSNe = length(list_SNe_GP)
print('Number of SNe in the 2_Selection/[ ]/Goods/ folder:')
numSNe

# Text to be put in the header of 'SN_list_template_*.txt':
text1 = '# List of supernovae used to construct this template.'
text2 = '# Created from data in directory:'
text3 = file.path(FilterSyst, SelectionPath, KindOfData, GoodPath)
text4 = sprintf('# Cutoff:  %.2f < z_cmb < %.2f', z_lowerLimit, z_upperLimit)
text5 = sprintf('#  %.2f < dm15 < %.2f | %.2f < EBVhost < %.2f | EBV_mw < %.2f',dm15LowerLim, dm15UpperLim,EBVhost_Min,EBVhost_Max,EBV_mw_Max)
text6 = sprintf('# Gaussian-Processes LC functions with data at phase=0 only, for the normalized template.')
text7 = '#	SN								COMMENTS'
  
# Go back to the '3_template' folder to write the 'SN_list_template_Norma.txt' text file
setwd(DirTemplate1)
getwd()

if (NormalizedTemp == TRUE) {
  # Write to the text file the headers created above:
  write.table(text1, file='SN_list_template_Norma_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(text2, file='SN_list_template_Norma_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text3, file='SN_list_template_Norma_.txt', quote=FALSE, row.names='#',   col.names=FALSE, append=TRUE)
  write.table(text4, file='SN_list_template_Norma_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text5, file='SN_list_template_Norma_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text6, file='SN_list_template_Norma_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text7, file='SN_list_template_Norma_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
} else {
  # Write to the text file the headers created above:
  write.table(text1, file='SN_list_template_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE)
  write.table(text2, file='SN_list_template_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text3, file='SN_list_template_.txt', quote=FALSE, row.names='#',   col.names=FALSE, append=TRUE)
  write.table(text4, file='SN_list_template_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text5, file='SN_list_template_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  # write.table(text6, file='SN_list_template_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
  write.table(text7, file='SN_list_template_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
}

#--------------

#     MAKING CUTOFFS: z_cmb, dm15, EBV_host, EBV_mw, T_Bmax 

NumSNeTemplate <-0

for(j in 1:numSNe){ 
  # Set the condition on the redshift: consider SNe with a z_CMB > cutoff only.
  z_cmb <- read.table(paste(DirLCData,list_SNe_GP[j], sep=''))$V1[1] # 
  dm15par <-  read.table(paste(DirLCData,list_SNe_GP[j], sep=''))$V1[2] #
  EBV_mw <-  read.table(paste(DirLCData,list_SNe_GP[j], sep=''))$V1[4] #
  EBVhost <-  read.table(paste(DirLCData,list_SNe_GP[j], sep=''))$V1[5] #
  T_Bmax <- read.table(paste(DirLCData,gsub('.txt', '', list_SNe_GP[j]), '_GP_mean_sigma_Filled.dat', sep=''))$V2[PhaseIndexZero] #
  # T_Bmax 
  # list_SNe_GP[j]
  
  if (z_cmb > z_lowerLimit & z_cmb < z_upperLimit & dm15par>dm15LowerLim & dm15par<dm15UpperLim & EBVhost>EBVhost_Min & EBVhost < EBVhost_Max  & EBV_mw < EBV_mw_Max) { 
    if (NormalizedTemp == TRUE) { 
      if (T_Bmax != 20) { # For Normalized template.
        write.table(list_SNe_GP[j], file='SN_list_template_Norma_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
        NumSNeTemplate <- NumSNeTemplate + 1 # Simple counter
      }
    } else { # For unnormalized template.
        write.table(list_SNe_GP[j], file='SN_list_template_.txt', quote=FALSE, row.names=FALSE, col.names=FALSE, append=TRUE)
        NumSNeTemplate <- NumSNeTemplate + 1 # Simple counter
      }
  } 
} 

cat('Number of SNe for this Normalized template:', NumSNeTemplate)

################################################################


# Set as working directory the folder "NIR_band/2_Selection/[sample]/Goods/"
setwd(DirLCData)
getwd() # Show the current directory

# Set a seed for repeatable plots
set.seed(12345)

# Some useful segments of paths to files
FilledPath = '_GP_mean_sigma_Filled.dat'
DiagonalSymb = ''
# txt_Extension = '.txt' #

# Defining the functions needed to compute mu.

# Eq. (5.20a) in Gelman's book
mu_hat <- function (tau, y, sigma){sum(y/(sigma^2 + tau^2))/sum(1/(sigma^2 + tau^2)) }
# Eq. (5.20b) in Gelman's book. This is the standard error^2 of the mean hyperparameter mu
V_mu <- function (tau, y, sigma){1/sum(1/(tau^2 + sigma^2)) } 

# Creation of a grid of values of tau (the standard deviation of the hyperpopulation):
n_grid <- 1500
# tau_grid <- seq (0.0005, 4, length=n_grid) # OLD 1
tau_grid <- seq (0.00001, 3, length=n_grid)
# tau_grid
# class(tau_grid)

#===============================================================



#         MAIN LOOP

if (NormalizedTemp == TRUE){
  txtFiles <- 'SN_list_template_Norma'
} else { txtFiles <- 'SN_list_template' }

# txtFiles
# txtFiles[2]

# OLD. Loop over the normalized and unnormalized lists
# OLD. This loop finishes at the end of this notebook at the bottom of the page.
# OLD. for (ii in 1:2){
  # ii<-1
  # cat(txtFiles[ii])

# Initializing
numSNeFinal <- 0

#   READ THE 'SN_list_template_Norma.txt' LIST.
# First check if 'SN_list_template_Norma_Notes_.txt' file exist, if so read that file otherwise read the one just created above.
try1 <- try(list_SN_names <- read.table(paste(DirTemplate1,'/', txtFiles,'_Notes_.txt', sep = '')))
if(inherits(try1, "try-error")) {
  list_SN_names <- read.table(paste(DirTemplate1,'/', txtFiles,'_.txt', sep = ''))
  # print('Reading data from SN_list_template_Norma_.txt instead')
  list_SN_names
}
numSNeFinal <- dim(list_SN_names)[1]
cat('Number of SNe:',numSNeFinal, ':')

# It runs over all the different epochs.
#-- IMPORTANT: EVERY TIME THAT I RUN THIS LOOP I HAVE TO RUN FIRST THE FOLLOWING LINES TO INITIALIZE THE VALUES OF "mu_tau_all_df" BECAUSE THE LOOP IS GOING TO APPEND DATA TO IT.
# IF I DO NOT INITIALIZE "mu_tau_all_df" THEN IT MAY ALREADY CONTAIN SOME DATA AND THEN THE LOOP IS GOING TO SIMPLE ADD MORE DATA AT THE END OF THE EXISTING ARRAY!

#--------

# FOR THE STANDARD ERROR OF MU
# Creation of a (0, 0) data frame to initialize the append in the loop below
mu_stdError_all_df <- data.frame(0, 0)
# mu_stdError_all_df
# Changing the headers of the data frame to match the headers of "tau_mu_all_df" below.
names(mu_stdError_all_df) <- c("median.mu.", "sd.mu.")
# mu_stdError_all_df

# FOR THE HYPERPARAMETER TAU
# Creation of a (0, 0) data frame to initialize the append in the loop below
mu_tau_all_df <- data.frame(0, 0)
# mu_tau_all_df
# Changing the headers of the data frame to match the headers of "tau_mu_all_df" below.
names(mu_tau_all_df) <- c("median.mu.", "tau_grid")
# mu_tau_all_df

#--------

for(k in 1:numPhases) # Loop over -phases-
# for(k in 1:65) # Loop over -phases-
{ 
  # Initializing some values for sanity:
  mu <- 0
  V <- 0
  log_p_tau <- rep (NA, n_grid)
  p_tau <- 0
  
  # Creating the arrays for y and sigma for a single epoch k:
  y_NA <- numeric(numSNeFinal)
  sigma_NA <- numeric(numSNeFinal)
  
  #---------------------------
  # LOOP OVER SNe.
  # Loop over the abs. magnitude for different SN at fixed epoch
  for(j in 1:numSNeFinal){ 
    
    # substr(list_SN_names$V1[10],1,50)
    CharacterSize <- nchar(substr(list_SN_names$V1[j],1,50))
    # CharacterSize
    trim1 <- CharacterSize-4
    SN_Name <- substring(list_SN_names$V1[j],1,trim1)
    SN_Name
    
    # Read the value of dm15 for a given SNe:
    # dm15Int <- read.table(paste(SN_Name, '.txt', sep=""))$V1[2]
    
    # Absolute magnitude at phase k:
    ytemp <- read.table(paste(SN_Name, FilledPath, sep=""))$V2[k] 
    
    # Absolute magnitude at phase = 0 (i.e., at T_Bmax):
    if (NormalizedTemp == TRUE) { # For normalized template
      absmag_TBmax <- read.table(paste(SN_Name, FilledPath, sep=""))$V2[PhaseIndexZero]
    } else {  # For unnormalized template
      absmag_TBmax <- 0
    }
    
    #- "ytemp<10" = if normalized magnitude is smaller than 10
    # I gave "20" to ranges without data in the light curve.
    if (ytemp<19) { # if "True" then the datum at this phase for the given SNe is GOOD.
      y_NA[j] <- ytemp - absmag_TBmax
      sigma_NA[j] <- read.table(paste(SN_Name, FilledPath, sep=""))$V3[k]
      # sigma_NA[j] <- read.table(paste(list_SN_names[j,1], FilledPath, sep=""))$V3[k] # OLD
    } else { # "else" = then the datum at this phase for the given SNe is BAD.
      y_NA[j] <- NA
      sigma_NA[j] <- NA
    }
  } # <--- END LOOP OVER SNe FOR A GIVEN PHASE
  
  # create new datasets without the missing data NA:
  y <-na.omit(y_NA)
  # cat(y)
  # cat("  ")
  sigma <- na.omit(sigma_NA)
  
  #-----
  #- "length(y) > 2" = "if there are data in y", i.e., "if for a given phase k there are magnitude data for at least 3 SNe, then":
  if (length(y) > 2) {
  
  # Loop to compute tau (the standard deviation hyperparameter) from the Eq. (5.21). For that, it is needed to  compute values of "mu_hat" and "V_mu"
  for (i in 1:n_grid){
    mu <- mu_hat (tau_grid[i], y, sigma)
    V <- V_mu (tau_grid[i], y, sigma)
    # Equation (5.21):
    log_p_tau[i] <- .5*log(V) - .5*sum(log(sigma^2 + tau_grid[i]^2)) -
      .5*sum((y-mu)^2/(sigma^2 + tau_grid[i]^2)) # Eq. (5.21)
  }
  
  #-----
  # Now, we compute the posterior density for tau on the log scale and rescale 
  # it to eliminate the possibility of computational overflow or 
  # underflow that can occur when multiplying many factors:
  log_p_tau <- log_p_tau - max(log_p_tau) #
  p_tau <- exp(log_p_tau)
  
  # The final PDF for tau:
  p_tau <- p_tau/sum(p_tau) # Normalizing the posterior PDF.
  # plot(tau_grid, p_tau)
  # dev.copy(png, 'p_tau_PDF.png') 
  # dev.off()
  
  # OLD  --->>
  # Creation of a data frame combining (tau_grid, p_tau), to look for the maximum and its corresponding value of tau.
  # p_tau_df <- data.frame(tau_grid, p_tau)
  # Plotting the PDF of tau at a given phase
  # uptodata <- 12
  # plot(tau_grid[1:uptodata], p_tau[1:uptodata])
  # plot(p_tau_df[1:uptodata]) 
  # dev.copy(png, 'p_tau_PDF.png') 
  # dev.off()
  # class(p_tau_df)
  # p_tau_df 
  # Finding the maximum likelihood estimate (MLE) (the mode of the PDF of tau) and its corresponding value of tau:
  # tau_MLE <- p_tau_df[which(p_tau_df$p_tau == max(p_tau_df$p_tau)), ]
  # tau_MLE
  # tau_MLE[1]
  # dim(tau_MLE[1])
  # class(tau_MLE[1])
  # <<--- End OLD
  
  #------- So far, all the calculations have been deterministics -------
  
  n_sims <- 3000 # Number of simulations per MJD day
  # Sampling from tau's PDF.
  tau <- sample (tau_grid, n_sims, replace=TRUE, prob=p_tau)  
  
  #-- Plotting the histogram of tau PDF --
  phase_int <- ((k-PhaseIndexZero)/2)
  # phase_int
  hist(tau, 25, xlab='Sample std dev (mag)') # Histogram of tau
  
  if (NormalizedTemp == TRUE){
    dev.copy(png, paste(DirTemplate1,'/','Plots_histo_Norma/','Plot_hist_phase_', phase_int, '.png', sep=''))
  } else {
    dev.copy(png, paste(DirTemplate1,'/','Plots_histo/','Plot_hist_phase_', phase_int, '.png', sep=''))
  }
  
  dev.off()
  
  # The median of the PDF of tau: this is the value that I write down in the text file and that I will use as the standard deviation for my NIR templates.
  tau_Median <- data.frame(median(tau))
  # Rename the column of data for consistency with 'mu_tau_all_df' (see the last definitions before the main loop). Otherwise I obtain an error and the loops breaks
  names(tau_Median) <- c('tau_grid') 
  tau_Median
  # dim(tau_Median)
  # class(tau_Median)

  J <-numSNeFinal
  
  # Computing the mean hyperparameter "mu":
  mu <- rep (NA, n_sims)
  # theta <- array (NA, c(n_sims,J))
  for (i in 1:n_sims){
    # PDF of mu. Equation above of Eq. (5.20): Sampling from a normal distribution for mu.
    mu[i] <- rnorm (1, mu_hat(tau[i],y,sigma), sqrt(V_mu(tau[i],y,sigma)))
    
    # PDF of theta's. Eqs. (5.17)
    # theta_mean <- (mu[i]/tau[i]^2 + y/sigma^2)/(1/tau[i]^2 + 1/sigma^2)
    # theta_sd <- sqrt(1/(1/tau[i]^2 + 1/sigma^2))
    # theta[i,] <- rnorm (J, theta_mean, theta_sd)
  }
  #-----
  
  # Creation of some data frames to append the data (median of "mu", standard error of "mu")
  mu_stdError_df <- data.frame(median(mu), sd(mu))
  # Creation of some data frames to append the data (mean "mu", tau)
  # mu_tau_df <- data.frame(median(mu), tau_MLE[1]) # OLD
  mu_tau_df <- data.frame(median(mu), tau_Median)
  
  median(mu)
  sd(mu)
  tau_Median
  
  # Relabeling the index of the data frame so that the first entry is equal to 1.
  # rownames(mu_stdError_df) <- 1:nrow(mu_stdError_df)
  rownames(mu_tau_df) <- 1:nrow(mu_tau_df)
  
  # Final data frame where the (mu, tau) is appending
  mu_stdError_all_df <- rbind(mu_stdError_all_df, mu_stdError_df)
  mu_tau_all_df <- rbind(mu_tau_all_df, mu_tau_df)
  
  # if (k<2){cat("Print loop with data")}
  cat(k) # Print on the screen the current loop
  cat(" ")
  # End of the "if" that is true when length(y) > 0.
  } else { # If there are phases with no data (the edges basically),  then create a line in the main file with Absolute magnitude = 0 and standar error/deviation = 0.1
    cat('-') # Print on the screen a dash for phases with no data.
    cat('   ')
     # Create a data frames with mean = 0 and variance = 0.1 .
     mu_stdError_df <- data.frame(0, 0.1)
     mu_tau_df <- data.frame(0, 0.1)
     # Changing the headers of the data frame to match the headers of "tau_mu_all_df" below.
     names(mu_stdError_df) <- c("median.mu.", "sd.mu.")
     names(mu_tau_df) <- c("median.mu.", "tau_grid")
     # Final data frame where the (mu, tau) is appending
     mu_stdError_all_df <- rbind(mu_stdError_all_df, mu_stdError_df)
     mu_tau_all_df <- rbind(mu_tau_all_df, mu_tau_df)
     # cat(k) # Print the current loop
   }
}
# end loop over phases
#===============================================================

# Adding the column of -phase- to the data frame of (mu, tau) so that I will 
# have an array with columns (phase, mu, tau).
# Create array of phases using the first SN data:
phases <- read.table(paste(SN_Name,FilledPath,sep=""))$V1[1:numPhases]
# phases <- read.table(paste(list_SN_names[1,1],FilledPath,sep=""))$V1[1:numPhases] # OLD
length(phases)

# Creating the combined data frame (phase, mu, stdError) and (phase, mu, tau).
# The first row in mu_tau_all_df is zero by construction (when I defined "mu_tau_all_df" before the loop), so I have to skip this row.
phase_mu_stdError <- cbind(phases, mu_stdError_all_df[c(1:numPhases+1),])
phase_mu_tau <- cbind(phases, mu_tau_all_df[c(1:numPhases+1),])
dim(phase_mu_stdError)
# phase_mu_tau

# The row 72 in 'mu_stdError_all_df' and 'mu_tau_all_df' corresponds to phase=0. I use the value of the absolute magnitude at phase=0 to subtract to the absolute magnitudes at different phases in order to obtain a normalized template.

# Construct the path to the location to save the outputs:
DirSaveOutputs <- DirTemplate1
cat('\n Directory to save the data: \n')
cat(DirSaveOutputs)

# Writting the results to a data file:
if (NormalizedTemp == TRUE){
  MyPathAndNameFile1 <- file.path(DirSaveOutputs, 'Template_phase_mu_stdError_FromR_Norma.dat')
  MyPathAndNameFile2 <- file.path(DirSaveOutputs, 'Template_phase_mu_tau_FromR_Norma.dat')
} else {
  MyPathAndNameFile1 <- file.path(DirSaveOutputs, 'Template_phase_mu_stdError_FromR.dat')
  MyPathAndNameFile2 <- file.path(DirSaveOutputs, 'Template_phase_mu_tau_FromR.dat')
}
# For (phase, mu, stdError)
write.table(phase_mu_stdError, file=MyPathAndNameFile1, sep="   ", row.names = FALSE, col.names = FALSE)

# For (phase, mu, tau)
write.table(phase_mu_tau, file=MyPathAndNameFile2, sep="   ", row.names = FALSE, col.names = FALSE)

#--------------

# cat('')
# cat('Number of SNe for this template:', NumSNeTemplate)

# }
# END LOOP OVER SNe LISTS

cat('--- All done ---')


#######################################################
