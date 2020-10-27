# Functions to plot a grid network with flows

getEdgeIndex <- function(edgeDf,a,b){
	edgeDf[,1]==a & edgeDf[,2]==b
}

convertToColorScale <- function(x, min = 0, max = 1, colorVector = c("green4", "yellow", "red")){
	myPalette <- colorRampPalette(colorVector)(length(x))#heat.colors(length(x))
	res <- myPalette[ findInterval( x, seq( min, max, length.out= (length(x)) ), rightmost.closed= T ) ]
	res
}

# for a route in a road graph, different data is needed:
# 	-current position of agent
#	-positions of cities to visit
#	-planned route
# Thus, nodes need to export:
#	-whether a node is a destination
#	-whether the agent is currently on that node
# Edges need the following information:
#	-does it exist?
#	-is it on the route?

#' @export
plotNetworkRoute <- function(g, nodeDf, edgeDf){

	COLOR_DESTINATION_NODE <- "indianred1"
	COLOR_CURRENT_NODE <- "blue"
	COLOR_PATH <- "red"
	COLOR_NEUTRAL <- "black"
	DEFAULT_EDGE_WIDTH <- 0.5
	VERTEX_SIZE_SMALL <- 1
	VERTEX_SIZE_MEDIUM <- 4
	VERTEX_SIZE_LARGE <- 8

	##################
	# data about vertices

	#print(g)
    if (is.null(g)){
	    g <- igraph::graph_from_data_frame(d=edgeDf, vertices=nodeDf)
    }
    #testGraph <- igraph::remove.vertex.attribute(testGraph, "name")

	vertexColors <- rep(COLOR_NEUTRAL, igraph::vcount(g))
	vertexLabels <- rep("", igraph::vcount(g))
	vertexSizes <- rep(VERTEX_SIZE_SMALL, igraph::vcount(g))


	#vertexLabels[nodeDf$type >= 1 & nodeDf$value != 0] <- nodeDf$value[nodeDf$type >= 1 & nodeDf$value != 0]
	##################
	# data about edges

	edgeColors <- rep("darkgrey", nrow(edgeDf))
	edgeWidths <- rep(DEFAULT_EDGE_WIDTH, nrow(edgeDf))
	edgeWidths[edgeDf$onRouteOf>=1] <- DEFAULT_EDGE_WIDTH + 2.5#edgeDf$onRouteOf*
    if ("onRouteOf" %in% colnames(edgeDf)){
    	maxRouteNumber <- max(edgeDf$onRouteOf)
    	relRouteNumber <- edgeDf$onRouteOf[edgeDf$onRouteOf>=1] / maxRouteNumber
    	edgeColors[edgeDf$onRouteOf>=1] <- convertToColorScale(relRouteNumber, colorVector = c("blue", "red"))#COLOR_PATH
    }

	maxValue <- max(nodeDf$value)
	relValue <- nodeDf$value / maxValue
	relValue <- relValue[relValue != 0]


	#vertexColors <- rep("darkorange", length(relValue))
	vertexColors[nodeDf$value != 0] <- convertToColorScale(relValue, colorVector = c("darkred", "yellow"))
	vertexSizes[nodeDf$value != 0] <- VERTEX_SIZE_MEDIUM

	subNodeDf <- nodeDf[nodeDf$type != 0, ]
	vertexColors[nodeDf$type == 2] <- COLOR_CURRENT_NODE

	vertexSizes[nodeDf$type == 2] <- VERTEX_SIZE_LARGE
	myLayout = as.matrix(nodeDf[,c(3,2)])

	igraph::plot.igraph(g, layout=myLayout, vertex.color=vertexColors,
		vertex.label.color="black", vertex.label=vertexLabels,
		vertex.size = vertexSizes,
		edge.curved=0,
		#vertex.size=10,
		# edge.curved=0,
		edge.color = edgeColors,
		edge.width= edgeWidths,
		vertex.frame.color="NA",
		margin = -0.1,
		edge.arrow.mode="0"
	)
}

#' @export
plotNetworkRouteNodesOnly <- function(g, nodeDf, edgeDf){

    COLOR_DESTINATION_NODE <- "indianred1"
    COLOR_CURRENT_NODE <- "blue"
    COLOR_PATH <- "red"
    COLOR_NEUTRAL <- "black"
    DEFAULT_EDGE_WIDTH <- 0.5
    VERTEX_SIZE_SMALL <- 1
    VERTEX_SIZE_MEDIUM <- 4
    VERTEX_SIZE_LARGE <- 6

    ##################
    # data about vertices

    #print(g)
    if (is.null(g)){
        g <- igraph::graph_from_data_frame(d=edgeDf, vertices=nodeDf)
    }
    #testGraph <- igraph::remove.vertex.attribute(testGraph, "name")

    vertexColors <- rep(COLOR_NEUTRAL, igraph::vcount(g))
    vertexLabels <- rep("", igraph::vcount(g))
    vertexSizes <- rep(VERTEX_SIZE_SMALL, igraph::vcount(g))


    #vertexLabels[nodeDf$type >= 1 & nodeDf$value != 0] <- nodeDf$value[nodeDf$type >= 1 & nodeDf$value != 0]
    ##################
    # data about edges

    edgeColors <- rep("darkgrey", nrow(edgeDf))
    edgeWidths <- rep(DEFAULT_EDGE_WIDTH, nrow(edgeDf))
    edgeWidths[edgeDf$onRouteOf>=1] <- DEFAULT_EDGE_WIDTH + 1.5#edgeDf$onRouteOf*
    if ("onRouteOf" %in% colnames(edgeDf)){
        maxRouteNumber <- max(edgeDf$onRouteOf)
        relRouteNumber <- edgeDf$onRouteOf[edgeDf$onRouteOf>=1] / maxRouteNumber
        edgeColors[edgeDf$onRouteOf>=1] <- convertToColorScale(relRouteNumber, colorVector = c("blue", "red"))#COLOR_PATH
    }

    maxValue <- max(nodeDf$value)
    relValue <- nodeDf$value / maxValue
    relValue <- relValue[relValue != 0]


    #vertexColors <- rep("darkorange", length(relValue))
    vertexColors[nodeDf$value != 0] <- convertToColorScale(relValue, colorVector = c("darkred", "yellow"))
    vertexSizes[nodeDf$value != 0] <- VERTEX_SIZE_MEDIUM

    subNodeDf <- nodeDf[nodeDf$type != 0, ]
    vertexColors[nodeDf$type == 2] <- COLOR_CURRENT_NODE

    vertexSizes[nodeDf$type == 2] <- VERTEX_SIZE_LARGE
    myLayout = as.matrix(nodeDf[,c(3,2)])

    igraph::plot.igraph(g, layout=myLayout, vertex.color=vertexColors,
                        vertex.label.color="black", vertex.label=vertexLabels,
                        vertex.size = vertexSizes,
                        edge.curved=0,
                        #vertex.size=10,
                        # edge.curved=0,
                        edge.color = "white",
                        edge.width= 0,
                        vertex.frame.color="NA",
                        margin = -0.1,
                        edge.arrow.mode="0"
    )
}

