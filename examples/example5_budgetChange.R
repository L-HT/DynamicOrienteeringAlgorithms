# example for improvement heuristics (that require an initial solution)
# the examples here also use dynamic OP instances (where changeFiles are used)

library(DynamicOrienteeringAlgorithms)
library(igraph)

# LGF-Leipzig-100-80000-1-evaluation-L50-absolute-C10-budget
# name of the instance
instanceName <- "LGF-Leipzig-100-80000-1"

# name of initial solution (available: random1, random2, greedy, gsr, ea)
initialSolutionName <- "random2"

# name of timeMeasure and choice of L (number of dynamic changes)
# bitVector = subset, evaluation = function evaluations
# high: L=V0, low: L=V0/2, none=static OP

changesName <- "evaluation-L50-absolute-C10-budget"
# changesName <- "static"

# use paste to concatenate strings

# read instance data
instanceData <- convertLGFtoR(paste("./instances/", instanceName, sep=""))

# set path to initial solution and dynamic changes using "paste"
pathToInitialSolution <- paste("./initialSolutions/", instanceName, "_", initialSolutionName, sep="")
pathToChanges <- paste("./instances_changes_budget/", instanceName,"-", changesName, sep="")

# name of the resulting logfile
problemName <- paste(instanceName, initialSolutionName, changesName, sep="-")

pathToDistanceMatrix <- paste("./distanceMatrices/", instanceName, sep="")

# handlingMode:
# 0 - no handling
# 1 - restart (H-restart)
# 2 - table (H-table)
# 31 - H-length
# 32 - H-value

# call VNS_DOP (improvement heuristic)
# an optional file suffix can be attached with the fileSuffix argument
callVNSImprover(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, 
                budget = instanceData$budget,
                problemName= problemName,
                runNumber=1,
                pathToInitialSolution = pathToInitialSolution,
                fileSuffix = "testSuffix",
                pathToChanges = pathToChanges, 
                pathToDistanceMatrix = pathToDistanceMatrix,
                budgetChangeHandlingMode = 1, 
                minBudgetToHandle = 30000,
                maxBudgetToHandle = 100000,
                budgetChangeTableSize = 10
                )



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
                  pathToChanges = pathToChanges,
                  pathToDistanceMatrix = pathToDistanceMatrix,
                  budgetChangeHandlingMode = 1,
                  minBudgetToHandle = 30000,
                  maxBudgetToHandle = 100000,
                  budgetChangeTableSize = 30
)
# call GSR (improvement heuristic)
for (i in 1:10){
    # print(i)
    try(
    callGraspSrImprover(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, 
              budget = instanceData$budget,
              problemName= problemName,
              runNumber=i,
              pathToInitialSolution = pathToInitialSolution,
              fileSuffix = "",
              pathToChanges = pathToChanges,
              pathToDistanceMatrix = pathToDistanceMatrix,
              budgetChangeHandlingMode = 31,
              minBudgetToHandle = 30000,
              maxBudgetToHandle = 100000,
              budgetChangeTableSize = 10
            ))
}
