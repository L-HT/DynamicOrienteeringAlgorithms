
#' Export to Lemon Graph Format
#' @export
exportRoadGraphToLGF <- function(filePath, g){

    # calculate ids as integers
    myIds = as.integer(substr(igraph::V(g)$id, 2, nchar(igraph::V(g)$id)))

    # create nodeDf
    nodeDf = data.frame("id" = myIds, "x" = igraph::V(g)$x, "y" = igraph::V(g)$y)

    # create arcDf
    arcDf = data.frame("from" = as.integer(igraph::head_of(g,igraph::E(g)))-1, "to" = as.integer(igraph::tail_of(g,igraph::E(g)))-1, "length" = igraph::E(g)$length)

    saveDigraphAsLGF(filePath, nodeDf, arcDf);
}


