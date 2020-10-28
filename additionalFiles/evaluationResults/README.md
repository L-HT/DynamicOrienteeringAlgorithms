These folders contain the AUC values from the experiments. They are structured as follows:

-improvementHeuristics: results of the experiments with improvement heuristics VNS_DOP, GSR', EA'
  -"bitVector" corresponds to the results with respect to the time measure "subset" (SS)
  -"evaluation" corresponds to the results with respect to the time measure "function evaluation" (FE)
  -if the file name contains the word "normalized", the file contains the AUC^{norm} values, otherwise it contains the AUC^{rel} values
  -if the file name contains "500000", the calculated values only consider the progress curves after 500 000 time units. 

-standaloneAlgorithms: results of the experiments with "standalone algorithms" VNS_DOP, GSR, EA
  -"bitVector" corresponds to the results with respect to the time measure "subset" (SS)
  -"evaluation" corresponds to the results with respect to the time measure "function evaluation" (FE)
  -if the file name contains the word "normalized", the file contains the AUC^{norm} values, otherwise it contains the AUC^{rel} values
  -the files with "EndPerformance" contain the values for the term "Value(P_best)/UB_t" reached at the end of the time horizon (after 1 000 000 time units have passed) for the static OP instances. These values are used for the statistical tests.
