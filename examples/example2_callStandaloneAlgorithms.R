# example for standalone algorithms (that do not require an initial solution)
# the test runs in this file work with static OP instances
# see "example3_callImprovementHeuristics.R" for examples with dynamic OP instances

library(DynamicOrienteeringAlgorithms)
library(igraph)

# read instanceData from file
instanceData <- convertLGFtoR("./instances/LGF-Leipzig-50-80000-1")

# call VNS_DOP (abbreviated as VNS) to calculate a solution with a budget of 15000
# see src/VNS.cpp for the source code
callVNSSolver(instanceData$nodeDf, 
                  instanceData$arcDf, 
                  instanceData$problemDf, 
                  budget = 15000,
                  problemName= "myTestProblem",
                  runNumber = 1234,
                  fileSuffix = ""
                  )

# check the folder "output" to see the log data for this run
# (file "myTestProblem-gsr-1234"'")

# check the folder "initialSolutions" to see the best solution generated during this run
# (file "myTestProblem_gsr")

# if you want to plot this solution:
# read the solution as a long string and split it into an array of numbers
testSolution <- readLines("./initialSolutions/myTestProblem_vns",warn = F)
testSolution <- as.integer(strsplit(testSolution, split = ",")[[1]])

# plot the solution
plotNetworkRoute_R(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, testSolution)

##########################################

# examples on how to call other algorithms
# (they also write log data into the aforementioned folders)
# plotting the solutions can also be done as described above

# simple greedy heuristic
# see NaiveMethod.cpp for the source code
naiveInitialSolution(instanceData$nodeDf, 
                     instanceData$arcDf, 
                     instanceData$problemDf, 
                     15000,
                     problemName="myTestProblem2"
)

# EA4OP (or just EA in the paper)
# see EA4OPSolver.cpp for the code
callEa4OpSolver(instanceData$nodeDf, 
                instanceData$arcDf, 
                instanceData$problemDf, 
                80000,
                problemName="myTestProblem3",
                runNumber = 9999
)

# GRASP-SR (GSR in the paper)
# see GraspSrSolver.cpp for the code
callGraspSrSolver(instanceData$nodeDf, 
                instanceData$arcDf, 
                instanceData$problemDf, 
                80000,
                problemName="myTestProblem4",
                runNumber = 9999
)




