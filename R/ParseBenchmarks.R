# functions for parsing some commonly used benchmarks
# only cyclic paths are considered here -> the end node is removed so that every path ends with the starting node

filePath <- "/home/htl/Documents/networks/Benchmarks/Chao/p1.2.a.txt"

#' @export
parseChao1 <- function(filePath){
    nodes <- read.table(filePath, header = F, sep = ";", quote = "", skip = 3, stringsAsFactors = F)
    nodeDf <- data.frame(id = 1:nrow(nodes), x = nodes[,1], y = nodes[,2])
    nodeDf <- nodeDf[-nrow(nodeDf),]
    problemDf <- data.frame(id = 1:nrow(nodes), type = c(2, rep(1, nrow(nodes)-1)), value = nodes[,3])
    problemDf <- problemDf[-nrow(problemDf),]
    return(list(
        nodeDf = nodeDf,
        arcDf = calculateArcDfEuclidian(nodeDf),
        problemDf = problemDf
    ))
}

filePath <- "/home/htl/Documents/networks/Benchmarks/Tsiligrides3/tsiligirides_problem_3_budget_025.txt"

#' @export
parseTsiligrides <- function(filePath){
    nodes <- read.table(filePath, header = F, sep = "\t", quote = "", skip = 1, stringsAsFactors = F)
    nodes <- nodes[-2,]
    nodeDf <- data.frame(id = 1:nrow(nodes), x = nodes[,1], y = nodes[,2])
    problemDf <- data.frame(id = 1:nrow(nodes), type = c(2, rep(1, nrow(nodes)-1)), value = nodes[,3])

    return(list(
        nodeDf = nodeDf,
        arcDf = calculateArcDfEuclidian(nodeDf),
        problemDf = problemDf
    ))
}

#' @export
calculateArcDfEuclidian <- function(nodeDf){
    result <- data.frame(from=rep(nodeDf$id, each=nrow(nodeDf)), to=rep(nodeDf$id, nrow(nodeDf)), length=1)
    result <- result[result$from != result$to,]
    result$length <- apply(
        result,
        MARGIN=1,
        FUN=function(row){
            row <- as.numeric(row)
            xDiff <- nodeDf$x[as.numeric(row[1])] - nodeDf$x[as.numeric(row[2])]
            yDiff <- nodeDf$y[as.numeric(row[1])] - nodeDf$y[as.numeric(row[2])]
            sqrt(xDiff*xDiff + yDiff*yDiff)
        }
    )
    return(result)
}

#' @export
plotEuclidianInstance <- function(nodeDf, arcDf, problemDf){
    nodeDf$type = 0
    nodeDf$type[problemDf$id] = problemDf$type
    nodeDf$value = 0
    nodeDf$value[problemDf$id] = problemDf$value
    g <- igraph::graph_from_data_frame(arcDf, vertices=nodeDf)
    #arcDf <- arcDf[1,]
    plotNetworkRouteNodesOnly(g, nodeDf, arcDf)
}

#' @export
plotEuclidianInstanceFromList <- function(myList){
    nodeDf <- myList$nodeDf
    arcDf <- myList$arcDf
    problemDf <- myList$problemDf
    plotEuclidianInstance(nodeDf, arcDf, problemDf)
}


#' @export
parseOPLib <- function(filePath){
    header <- read.table(filePath, header = F, sep = ":", quote = "", nrows = 6, fill = TRUE, stringsAsFactors = F)
    budget <- as.numeric(header[grepl("COST_LIMIT", header$V1), 2])

    nodes <- read.table(filePath, header = F, sep = " ", quote = "", skip = 6, fill = TRUE, stringsAsFactors = F)
    invalidPositions <- which(nodes$V1 %in% c("NODE_COORD_SECTION", "NODE_SCORE_SECTION", "DEPOT_SECTION", "EOF"))

    nodeDf <- nodes[(invalidPositions[1]+1):(invalidPositions[2]-1),]
    problemDf <- nodes[(invalidPositions[2]+1):(invalidPositions[3]-1),c(1,2)]

    # problemDf[,1] <- as.numeric(problemDf[,1])
    problemDf <- cbind(problemDf, "type" = c(2,rep(1, nrow(problemDf)-1)))
    problemDf <- problemDf[,c(1,3,2)]

    for (i in 1:ncol(nodeDf)){
        nodeDf[,i] <- as.numeric(nodeDf[,i])
    }
    for (i in 1:ncol(problemDf)){
        problemDf[,i] <- as.numeric(problemDf[,i])
    }

    colnames(nodeDf) <- c("id", "x", "y")
    colnames(problemDf) <- c("id", "type", "value")

    # remove the value of the start node
    problemDf[1, "value"] <- 0
    print(paste("file:", filePath))
    print(paste(nrow(nodeDf), "nodes"))
    # if (nrow(nodeDf) > 1500){
    #     return(NULL)
    # }
    return(list(
        nodeDf = nodeDf,
        arcDf = calculateArcDfEuclidian(nodeDf),
        problemDf = problemDf,
        budget = budget
    ))
}

filePath <- "/home/htl/Documents/networks/Benchmarks/OPLib-master/instances/gen4/brazil58-gen4-45.oplib"
# upper: brazil58-gen4-45.oplib
# lower: pa561-gen4-90.oplib

#' @export
parseOPLibNonEuclidian <- function(filePath, type = "upperMatrix", withDiagonal = FALSE){
    header <- read.table(filePath, header = F, sep = ":", quote = "", nrows = 7, fill = TRUE, stringsAsFactors = F)
    if (gsub(" ","", header[grepl("EDGE_WEIGHT_TYPE", header$V1), 2]) != "EXPLICIT") {
        stop("distance matrix not found")
    }
    budget <- as.numeric(header[grepl("COST_LIMIT", header$V1), 2])

    nodes <- read.table(filePath, header = F, sep = c(" ", "\t"), quote = "", skip = 7, fill = TRUE, stringsAsFactors = F)

    invalidPositions <- which(nodes$V1 %in% c("EDGE_WEIGHT_SECTION", "DISPLAY_DATA_SECTION", "NODE_SCORE_SECTION", "DEPOT_SECTION", "EOF"))
    #invalidPositions <- which(!substring(nodes$V1, 1,1) %in% 0:9)

    distMatrix <- nodes[(invalidPositions[1]+1):(invalidPositions[2]-1),]#ncol(nodes):1]
    problemDf <- nodes[(invalidPositions[2]+1):(invalidPositions[3]-1),c(1,2)]



    if (type == "upperMatrix"){
        for (i in 1:(ncol(distMatrix)-1)){
            for (j in ncol(distMatrix):(i+1)){
                distMatrix[i,j] <- distMatrix[i,j-i]
            }
            distMatrix[i,i] <- 0

            if (i > 1){
                for (j in (i-1):1){
                    distMatrix[i,j] <- distMatrix[j,i]
                }
            }
            # print(distMatrix[i,])
        }
        lastRow <- c(distMatrix[,ncol(distMatrix)], 0)
        distMatrix <- rbind(distMatrix, lastRow)
    }
    if (type == "lowerMatrix"){
        invalidPositions <- which(nodes$V1 %in% c("EDGE_WEIGHT_SECTION", "DISPLAY_DATA_SECTION", "NODE_SCORE_SECTION", "DEPOT_SECTION", "EOF"))

        startPos <- invalidPositions[1]+1

        firstRow <- c(0, distMatrix[,1])
        distData <- nodes[(invalidPositions[1]+1):(invalidPositions[2]-1),]#ncol(nodes):1]
        distData <- as.matrix(distData)
        distDataVector <- as.vector(t(distData))
        distDataVector <- distDataVector[distDataVector != ""]


        if (withDiagonal){
            # you get this formula by solving length(distDataVector) = S = 1+2+...+n for n
            n <- (sqrt(1+8*length(distDataVector))-1) / 2
            distMatrix <- matrix(0, nrow=n, ncol=n)
            counter <- 1
            for (i in 1:n){
                for (j in 1:i){
                    distMatrix[i,j] <- as.numeric(distDataVector[counter])
                    distMatrix[j,i] <- distMatrix[i,j]
                    counter <- counter + 1
                }
            }
        } else {
            # you get this formula by solving length(distDataVector) = S = 1+2+...+(n-1) for n
            n <- (sqrt(1+8*length(distDataVector))+1) / 2
            distMatrix <- matrix(0, nrow=n, ncol=n)
            counter <- 1
            for (i in 2:n){
                for (j in 1:(i-1)){
                    distMatrix[i,j] <- as.numeric(distDataVector[counter])
                    distMatrix[j,i] <- distMatrix[i,j]
                    counter <- counter + 1
                }
            }
        }

        # for (i in 2:(ncol(distMatrix))){
        #     if (i < ncol(distMatrix)){
        #         for (j in (i+1):ncol(distMatrix)){
        #             distMatrix[i,j] <- distMatrix[j,i]
        #         }
        #     }
        #     distMatrix[i,i] <- 0
        #     print(distMatrix[i,])
        # }
        # row.names(distMatrix) <- 1:nrow(distMatrix)
    }




    # problemDf[,1] <- as.numeric(problemDf[,1])
    problemDf <- cbind(problemDf, "type" = c(2,rep(1, nrow(problemDf)-1)))
    problemDf <- problemDf[,c(1,3,2)]

    nodeDf <- data.frame("id" = 1:nrow(problemDf), "x" = 0, "y" = 0)


    for (i in 1:ncol(nodeDf)){
        nodeDf[,i] <- as.numeric(nodeDf[,i])
    }
    for (i in 1:ncol(problemDf)){
        problemDf[,i] <- as.numeric(problemDf[,i])
    }

    colnames(nodeDf) <- c("id", "x", "y")
    colnames(problemDf) <- c("id", "type", "value")

    # remove the value of the start node
    problemDf[1, "value"] <- 0

    print(paste("file:", filePath))
    print(paste(nrow(nodeDf), "nodes"))
    # if (nrow(nodeDf) > 1500){
    #     return(NULL)
    # }

    arcDf <- data.frame(from=rep(nodeDf$id, each=nrow(nodeDf)), to=rep(nodeDf$id, nrow(nodeDf)), length=1)
    arcDf <- arcDf[arcDf$from != arcDf$to,]
    arcDf$length <- apply(
        arcDf,
        MARGIN=1,
        FUN=function(row){
            row <- as.numeric(row);
            as.numeric(distMatrix[row[1], row[2]])
        }
    )

    return(list(
        nodeDf = nodeDf,
        arcDf = arcDf,
        problemDf = problemDf,
        budget = budget
    ))
}

# lower with weird format: hk48-gen4-80.oplib
#' @export
parseOPLibLowerWeird <- function(filePath, type = "upperMatrix", withDiagonal = FALSE){
    header <- read.table(filePath, header = F, sep = ":", quote = "", nrows = 7, fill = TRUE, stringsAsFactors = F)
    if (gsub(" ","", header[grepl("EDGE_WEIGHT_TYPE", header$V1), 2]) != "EXPLICIT") {
        stop("distance matrix not found")
    }
    budget <- as.numeric(header[grepl("COST_LIMIT", header$V1), 2])

    nodes <- read.table(filePath, header = F, sep = c(" ", "\t"), quote = "", skip = 7, fill = TRUE, stringsAsFactors = F)

    invalidPositions <- which(nodes$V1 %in% c("EDGE_WEIGHT_SECTION", "DISPLAY_DATA_SECTION", "NODE_SCORE_SECTION", "DEPOT_SECTION", "EOF"))
    #invalidPositions <- which(!substring(nodes$V1, 1,1) %in% 0:9)

    distMatrix <- nodes[(invalidPositions[1]+1):(invalidPositions[2]-1),]#ncol(nodes):1]

    startPos <- invalidPositions[1]+1

    firstRow <- c(0, distMatrix[,1])
    distData <- nodes[(invalidPositions[1]+1):(invalidPositions[2]-1),]#ncol(nodes):1]
    distData <- as.matrix(distData)
    distDataVector <- as.vector(t(distData))
    distDataVector <- distDataVector[distDataVector != "" & !is.na(distDataVector)]

    if (withDiagonal){
        # you get this formula by solving length(distDataVector) = S = 1+2+...+n for n
        n <- (sqrt(1+8*length(distDataVector))-1) / 2
        distMatrix <- matrix(0, nrow=n, ncol=n)
        counter <- 1
        for (i in 1:n){
            for (j in 1:i){
                distMatrix[i,j] <- as.numeric(distDataVector[counter])
                distMatrix[j,i] <- distMatrix[i,j]
                counter <- counter + 1
            }
        }
    } else {
        # you get this formula by solving length(distDataVector) = S = 1+2+...+(n-1) for n
        n <- (sqrt(1+8*length(distDataVector))+1) / 2
        distMatrix <- matrix(0, nrow=n, ncol=n)
        counter <- 1
        for (i in 2:n){
            for (j in 1:(i-1)){
                distMatrix[i,j] <- as.numeric(distDataVector[counter])
                distMatrix[j,i] <- distMatrix[i,j]
                counter <- counter + 1
            }
        }
    }
    row.names(distMatrix) <- 1:nrow(distMatrix)

    # problemDf[,1] <- as.numeric(problemDf[,1])
    problemDf <- nodes[(invalidPositions[3]+1):(invalidPositions[4]-1),c(1,2)]
    problemDf <- cbind(problemDf, "type" = c(2,rep(1, nrow(problemDf)-1)))
    problemDf <- problemDf[,c(1,3,2)]

    nodeDf <- data.frame("id" = 1:nrow(problemDf), "x" = 0, "y" = 0)


    for (i in 1:ncol(nodeDf)){
        nodeDf[,i] <- as.numeric(nodeDf[,i])
    }
    for (i in 1:ncol(problemDf)){
        problemDf[,i] <- as.numeric(problemDf[,i])
    }

    colnames(nodeDf) <- c("id", "x", "y")
    colnames(problemDf) <- c("id", "type", "value")

    # remove the value of the start node
    problemDf[1, "value"] <- 0

    print(paste("file:", filePath))
    print(paste(nrow(nodeDf), "nodes"))
    # if (nrow(nodeDf) > 1500){
    #     return(NULL)
    # }

    arcDf <- data.frame(from=rep(nodeDf$id, each=nrow(nodeDf)), to=rep(nodeDf$id, nrow(nodeDf)), length=1)
    arcDf <- arcDf[arcDf$from != arcDf$to,]
    arcDf$length <- apply(
        arcDf,
        MARGIN=1,
        FUN=function(row){
            row <- as.numeric(row);
            as.numeric(distMatrix[row[1], row[2]])
        }
    )

    return(list(
        nodeDf = nodeDf,
        arcDf = arcDf,
        problemDf = problemDf,
        budget = budget
    ))
}



