#' @export
parseChangeFileName <- function(s){
    subStrings <- strsplit(s, "-")[[1]]
    minIndex <- min(which(subStrings %in% c("evaluation", "bitVector", "distanceEvaluation")))

    # instance with dynamic
    if (length(subStrings) - minIndex + 1 == 2) {
        instanceName <- paste(subStrings[1:(minIndex-1)], collapse="-")
        timeMeasure <- subStrings[minIndex + 0]
        dynamicLevel <- subStrings[minIndex + 1]
    } else {
        if (length(subStrings) - minIndex + 1 == 1) {
            instanceName <- paste(subStrings[1:(minIndex-1)], collapse="-")
            timeMeasure <- NA
            dynamicLevel <- "none"
        } else {
            stop("Format of output file not correct. Output file could not be parsed")
        }
    }

    return(
        list(
            "instanceName" = instanceName,
            "timeMeasure" = timeMeasure,
            "dynamicLevel" = dynamicLevel
        )
    )
}

#' @export
parseOutputFileName <- function(s){
    initialSolvers <- c("greedy", "ea", "gsr", "random1", "random2")
    subStrings <- strsplit(s, "-")[[1]]
    minIndex <- min(which(subStrings %in% initialSolvers))

    # instance with dynamic
    if (any(subStrings == "none")) {
        if (length(subStrings) - minIndex + 1 == 5) {
            instanceName <- paste(subStrings[1:(minIndex-1)], collapse="-")
            initialSolution <- subStrings[minIndex]
            timeMeasure <- NA
            dynamicLevel <- "none"
            algorithm <- subStrings[minIndex + 2]
            runNumber <- subStrings[minIndex + 3]
            fileSuffix <- subStrings[minIndex + 4]
        } else {
            if (length(subStrings) - minIndex + 1 == 4) {
                instanceName <- paste(subStrings[1:(minIndex-1)], collapse="-")
                initialSolution <- subStrings[minIndex]
                timeMeasure <- NA
                dynamicLevel <- "none"
                algorithm <- subStrings[minIndex + 2]
                runNumber <- subStrings[minIndex + 3]
                fileSuffix <- ""
            } else {
                stop("Format of output file not correct. Output file could not be parsed")
            }
        }
    } else {
        if (length(subStrings) - minIndex + 1 == 6) {
            instanceName <- paste(subStrings[1:(minIndex-1)], collapse="-")
            initialSolution <- subStrings[minIndex]
            timeMeasure <- subStrings[minIndex + 1]
            dynamicLevel <- subStrings[minIndex + 2]
            algorithm <- subStrings[minIndex + 3]
            runNumber <- subStrings[minIndex + 4]
            fileSuffix <- subStrings[minIndex + 5]
        } else {
            if (length(subStrings) - minIndex + 1 == 5) {
                instanceName <- paste(subStrings[1:(minIndex-1)], collapse="-")
                initialSolution <- subStrings[minIndex]
                timeMeasure <- subStrings[minIndex + 1]
                dynamicLevel <- subStrings[minIndex + 2]
                algorithm <- subStrings[minIndex + 3]
                runNumber <- subStrings[minIndex + 4]
                fileSuffix <- ""
            } else {
                stop("Format of output file not correct. Output file could not be parsed")
            }
        }
    }

    return(
        list(
            "instanceName" = instanceName,
            "initialSolution" = initialSolution,
            "timeMeasure" = timeMeasure,
            "dynamicLevel" = dynamicLevel,
            "algorithm" = algorithm,
            "runNumber" = as.integer(runNumber),
            "fileSuffix" = fileSuffix
        )
    )
}

#' @export
parseInstanceName <- function(fileName){
    parts <- as.list(strsplit(fileName, "-")[[1]])
    if (parts[[1]] == "LGF") {
        parts <- parts[-1]
    }


    if (grepl(pattern = "oplib", fileName)){

        # names(parts) <- c("instanceName", "numberOfNodes", "budget", "instanceNumber", "additionalInfo")
        numberOfNodes <- stringr::str_extract(parts[[1]], "\\d+")

        result <- list()
        result["instanceName"] <- parts[[1]]
        result["numberOfNodes"] <- as.integer(numberOfNodes)
        result["budget"] <- as.numeric(parts[[4]])
        result["additionalInfo"] <- paste(parts[[2]], parts[[3]], sep="_")

        return(result)
    } else {
        names(parts) <- c("instanceName", "numberOfNodes", "budget", "instanceNumber")
        parts$numberOfNodes <- as.integer(parts$numberOfNodes)
        parts$budget <- as.numeric(parts$budget)
        parts$instanceNumber <- as.numeric(parts$instanceNumber)
        return(parts)
    }

    # return(parts)
}
