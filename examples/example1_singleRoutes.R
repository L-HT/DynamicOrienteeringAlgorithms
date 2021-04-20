# examples on how single solutions are plotted

library(DynamicOrienteeringAlgorithms)
library(igraph)

# read instanceData from file
instanceData <- convertLGFtoR("./instances/LGF-Leipzig-50-80000-1")

# use print(instanceData$problemDf) to see the IDs of the colored nodes
# the start node has type=2; one way to get the ID of the start node is:

startNode <- instanceData$problemDf[instanceData$problemDf$type==2,"id"]

# plot the graph without any routes
plotNetworkRoute_R(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, c(startNode))

# plot the graph without colored nodes (the warning can be ignored... for now)
# plotNetworkRoute(NULL, instanceData$nodeDf, instanceData$arcDf)

# Here is an example route. 
myRoute <- c(3878, 35, 276, 3878)

# plot the example route
plotNetworkRoute_R(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, myRoute)

# the routes are traversed from "cold" to "hot", or from blue to red

# feel free to test and plot different routes, for example this one:
# (only routes between nodes from problemDf are supported, otherwise an error occurs)
myRoute <- c(6056, 2333, 1835, 3878)

# this function only calculates the length and value of myRoute (without plotting)
callEvaluator(instanceData$nodeDf, instanceData$arcDf, instanceData$problemDf, 0, myRoute)
