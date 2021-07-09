These folders contain the AUC values from the experiments as well as the data used for plotting the averaged progress curves from which the AUC values are calculated. 
They are structured as follows:

 - "results_dynamicBudget" contains the results for the experiments with dynamic budgets:
 	- scatterplots comparing the different handling mechanisms
 	- plots showing the median performance (AUC^norm^) for different values of L and C
 	- the folder "LC-Tables" contains the average performance for each pair of (L,C) for both time measures, using the best handling mechanism for each instance
 	- averaged data used for plotting the progress curves 

 - "results_dynamicNodeValues" contains the results for the experiments with dynamic node values
 	- plots showing the median performance (rescaled AUC^rel^) for different values of L and C
 	- the folder "LC-Tables" contains the average performance for each pair of (L,C) for both time measures
 	- averaged data used for plotting the progress curves 

  - the csv-files starting with "AUC_PC" contain AUC values. More precisely, if the file name contains the word "normalized", the file contains the AUC^norm^ values, otherwise it contains the AUC^rel^ values
  - if the file name contains "500000", the calculated values only consider the progress curves after 500 000 time units. 

  - "bitVector" corresponds to the results with respect to the time measure "subset" (SS)
  - "evaluation" corresponds to the results with respect to the time measure "function evaluation" (FE)