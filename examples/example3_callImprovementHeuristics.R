# example for improvement heuristics (that require an initial solution)
# the examples here also use dynamic OP instances (where changeFiles are used)

library(DynamicOrienteeringAlgorithms)
library(igraph)

# name of the instance
instanceName <- "LGF-Leipzig-100-80000-1"

# name of initial solution (available: random1, random2, greedy, gsr, ea)
initialSolutionName <- "random2"

# name of timeMeasure and choice of L (number of dynamic changes)
# bitVector = subset, evaluation = function evaluations
# high: L=V0, low: L=V0/2, none=static OP

changesName <- "bitVector-high"
# changesName <- "static"

# use paste to concatenate strings

# read instance data
instanceData <- convertLGFtoR(paste("./instances/", instanceName, sep=""))

# set path to initial solution and dynamic changes using "paste"
pathToInitialSolution <- paste("./initialSolutions/", instanceName, "_", initialSolutionName, sep="")
pathToChanges <- paste("./instances_changes/", instanceName,"-", changesName, sep="")

# name of the resulting logfile
problemName <- paste(instanceName, initialSolutionName, changesName, sep="-")

# call VNS_DOP (improvement heuristic)
# an optional file suffix can be attached with the fileSuffix argument
callVNSImprover(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, 
                budget = 80000,
                problemName= problemName,
                runNumber=1,
                pathToInitialSolution = pathToInitialSolution,
                fileSuffix = "testSuffix",
                pathToChanges = pathToChanges)

# check the "output" folder to see the logdata for this run
# (file: LGF-Leipzig-50-80000-1-random2-bitVector-high-vns-3-testSuffix)

# see the "solutions" folder to see the best solution calculated at the end of the run
# (file: LGF-Leipzig-50-80000-1-random2-bitVector-high_testSuffix_vns)
# plotting these solutions can be done similarly to "example2_callStandaloneAlgorithms.R"

##################################

# examples with other algorithms

# call EA (improvement heuristic)
callEa4OpImprover(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, 
                  budget = instanceData$budget,
                  problemName= problemName,
                  runNumber=2,
                  pathToInitialSolution = pathToInitialSolution,
                  fileSuffix = "",
                  pathToChanges = pathToChanges)

# call GSR (improvement heuristic)
callGraspSrImprover(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, instanceData$budget,
                    problemName= problemName,
                    runNumber=3,
                    pathToInitialSolution = pathToInitialSolution,
                    fileSuffix = "",
                    pathToChanges = pathToChanges)
