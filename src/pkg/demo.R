#------------------------------------#
# Demo for Kinect data analysis tool #
#------------------------------------#


#############
# IMPORTANT #
#############

# Install dependencies

# The tool requires 'randomForest' package to be installed.
# If not installed, please run the following command:
# install.packages('randomForest')

# The tool requires 'rChoiceDialogs' package to be installed.
# If not installed, please run the following command:
# install.packages('rChoiceDialogs')




#load analysis package
source('kinect_analysis_tool.R')
# The maim function 'mainProcess()' takes as an input one single file corresponding to standard walk or heel-toe.
# It will automatically search for the corresponding other file, in the same directory.
# The function also allows missing data.
# After processing both files by removing outlier measurments, predictive models are applied to the data and return 3 scores between 0 and 1.
# One score for standard walk / one score for heel-toe / one aggregated score.
# Each score can be viewed as a probability of being impaired.
# In case of missing data, score is equal to 0.5.


# Details of optional parameters
# mainProcess <- function(file.std = file.choose(), file.ht = NULL, rmSteady = TRUE, rmOdd = TRUE, rmGlitch = TRUE, rotated = TRUE, verbose = FALSE)
# rmSteady: boolean for removing steady start/end of records (default: TRUE)
# rmOdd: boolean for removing odd movements at start/end of records (default: TRUE)
# rmGlitch: boolean for removing glitches occuring in records (default: TRUE)
# rotated: boolean for correcting camera angle variation (default: TRUE)
# verbose: if TRUE, display how data looks at each processing step (default: FALSE)

# Demo 1: Single file
mainProcess(file.std = 'data/144/144_1_W.csv') 

# Demo 2: no file provided
mainProcess()

# Demo 3: display data processing
mainProcess(file.std = 'data/144/144_1_W.csv', verbose = TRUE)

# Demo 4: provide heel-toe data
mainProcess(file.ht = 'data/144/144_1_HT.csv')

# Demo 5: missing heel-toe data
mainProcess(file.std = 'data/145/145_1_W.csv')

# Demo 6: more than one individual (select multiple files)
mainProcess(file.std = c('data/144/144_1_W.csv','data/144/144_2_W.csv'))

# Demo 7: output file
mainProcess(output.file='test/one_individual.csv')

# Demo 8: output file for multiple inputs
mainProcess(file.std = c('data/144/144_1_W.csv','data/144/144_2_W.csv'),output.file='test/multi_individuals.csv')

# Demo 8: output file for multiple chosen inputs
mainProcess(output.file='test/multi_individuals.csv')
