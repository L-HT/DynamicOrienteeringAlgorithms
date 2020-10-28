# Preprocessing in order to obtain graph data in GraphML format from OSM data

# JSON des OSM-Parsers auslesen (dieses Mal mit osm-Dateien aus dem Internet)
# nodeDf erstellen mit:
#   -Position des Knotens
# edgeDf mit
#   -Knoten i
#   -Knoten j
#   -Bogenlänge
# die Bögen sind als Adjazenzliste gespeichert?

myNames = c("Leipzig", "Berlin")
for (myName in myNames){

    # commandName = paste("java -jar ../osmparser-0.13.jar -f ../", myName,
    #                     ".osm -i highway -e  addr:street addr:city cycleway barrier amenity natural",
    #                     "public_transport tactile_paving bus tram railway service crossing man_made network stop shop ",
    #                     "bridge path bicycle_road oneway traffic_sign bicycle surface=paving_stones",
    #                     " -o ../output",
    #                     myName, sep="")

    # commandName = paste("java -jar ../osmparser-0.13.jar -f ../",
    # myName, ".osm -i highway=primary highway=secondary highway=tertiary highway=residential -e highway=unclassified ",
    #     "-e ",
    #     "highway=motorway_link highway=trunk_link highway=primary_link highway=secondary_link highway=tertiary_link ",
    #     "highway=living_street highway=service highway=pedestrian highway=footway railway",
    #     " -o ../output", myName, sep="")

	# von osm zu json
	
    commandName = paste("java -jar ./osmparse/osmparser-0.13.jar -f ./osmparse/",
                        myName, ".osm -i highway=primary highway=secondary highway=tertiary highway=residential",
                        " -e asd",
                        " -o ./osmparse/output", myName, sep="")

    # commandName = paste("java -jar ../osmparser-0.13.jar -f ../",
    #                     myName, ".osm -i highway=footway  -o ../output", myName, sep="")

    system(commandName)

    library(rjson)

    json_data <- fromJSON(file=paste("./osmparse/output", myName, ".json", sep=""))

    minLong = -Inf
    maxLong = Inf
    minLat = -Inf
    maxLat = Inf
    
    # Grenzen festlegen (damit leichter reproduzierbar) 
    if (myName == "Koethen"){
        # Köthen (small)
        minLong = 11.93699
        maxLong = 11.99998
        minLat = 51.72771
        maxLat = 51.76380
    }

    if (myName == "Leipzig"){
        # Leipzig (mid)
        minLong = 12.3042
        maxLong = 12.4707
        minLat = 51.2937
        maxLat = 51.3888
    }
    # 
    if (myName == "Berlin"){
        # Berlin (big)
        minLong = 12.9714
        maxLong = 13.8016
        minLat = 52.2879
        maxLat = 52.7167
    }


    longLatFit = function(long, lat){

        long < maxLong & long > minLong & lat < maxLat & lat > minLat
    }

    geodeticDf = data.frame(id=0:(length(json_data)-1), latitude=1:length(json_data), longitude=0, numberOfArcs = 0)

    # nodeDf füllen
    # Position mit Mercartor-Projektion schätzen

    numberOfArcs = 0
    
    for (i in 1:length(json_data)){
        if (i %% round(length(json_data)/100) == 0){
            cat(".")
        }
        if (i %% round(length(json_data)/10) == 0){
            print(";")
        }
        geodeticDf$latitude[i] = json_data[[i]]$la
        geodeticDf$longitude[i] = json_data[[i]]$lo
        geodeticDf$numberOfArcs[i] = length(json_data[[i]]$e)
    }

    print("")
    
    geodeticDf = geodeticDf[longLatFit(geodeticDf$longitude, geodeticDf$latitude),]
    nodeDf = geodeticDf
    colnames(nodeDf)[c(2,3)] = c("x","y")
    
    # von Longitude/Latitude zu XY durch Mercator-Projektion
    nodeDf[,c(2,3)] = geosphere::mercator(geodeticDf[,c(2,3)])

    nodeDf$x = nodeDf$x - min(nodeDf$x)
    nodeDf$y = nodeDf$y - min(nodeDf$y)

    edgeDf = data.frame(edgeID=1:sum(nodeDf$numberOfArcs), id1=0, id2=0, length=0)
    counter = 1
    bigCounter = 1
    # edgeDf füllen
    for (i in 1:length(json_data)){
        if (i %% round(length(json_data)/100) == 0){
            cat("*")
        }
        if (i %% round(length(json_data)/10) == 0){
            print("#")
        }
        #Schleife über die inneren Elemente
        for (j in 1:length(json_data[[i]]$e)){
            targetIndex = json_data[[i]]$e[[j]]$i

            if (longLatFit(json_data[[i]]$lo, json_data[[i]]$la) &
                longLatFit(json_data[[targetIndex+1]]$lo, json_data[[targetIndex+1]]$la)
                ){
                edgeDf[counter, c(2)] = as.integer(i-1)
                edgeDf[counter, c(3)] = as.integer(targetIndex)
                edgeDf[counter, c(4)] = json_data[[i]]$e[[j]]$w

                counter = counter + 1
            }
        }
        bigCounter = bigCounter + 1
    }
    print("")
    
    edgeDf <- edgeDf[1:(counter-1),]
    
    # Umwandlung cm zu m
    edgeDf$length = edgeDf$length / 100

    edgeDf$id1 <- as.integer(edgeDf$id1)
    edgeDf$id2 <- as.integer(edgeDf$id2)
    # Redundanzen raus
    #edgeDf = edgeDf[edgeDf$id1 < edgeDf$id2,]

	#########
	# Konvertierung zu graphml
	#########
	
    library(igraph)
    myGraph = graph_from_data_frame(d = edgeDf[,-1], directed=T, vertices = nodeDf)
    E(myGraph)$from = edgeDf[,2]
    E(myGraph)$to = edgeDf[,3]
    myLayout = as.matrix(nodeDf[,c(3,2)])
    #plot.igraph(x = myGraph, layout = myLayout, vertex.size=0.1, curved=0.6, vertex.label="" )

	# unerreichbare Komponenten entfernen
    g = myGraph
    V(g)$comp = components(g)$membership
    mySubgraph = induced_subgraph(g,V(g)$comp==1)
    mySubLayout = myLayout[V(g)$comp==1,]

    # hier wird alte ID gelöscht, sodass nur noch die von igraph zugewiesene ID übrig bleibt (die ab 1 kontinuierlich hochzählt)
    mySubgraph = delete_vertex_attr(mySubgraph, "name")
    mySubgraph = delete_vertex_attr(mySubgraph, "comp")
    mySubgraph = delete_edge_attr(mySubgraph, "from")
    mySubgraph = delete_edge_attr(mySubgraph, "to")

    write_graph(mySubgraph,file = paste("./osmparse/", myName,".graph",sep=""), format="graphml")

	# Plot des Graphen zur Veranschaulichung
	
    edgeDf = edgeDf[edgeDf$id1 < edgeDf$id2,]
    myGraph = graph_from_data_frame(d = edgeDf[,-1], directed=F, vertices = nodeDf)
    myLayout = as.matrix(nodeDf[,c(3,2)])
    g = myGraph
    V(g)$comp = components(g)$membership
    mySubgraph = induced_subgraph(g,V(g)$comp==1)
    mySubLayout = myLayout[V(g)$comp==1,]
    #delete_vertex_attr(mySubgraph, "name")
    #delete_vertex_attr(mySubgraph, "comp")

    png(filename = paste("./osmparse/",myName,"Plot.png", sep=""))
    plot(mySubgraph, layout=mySubLayout, vertex.size=0.1, vertex.label="")
    dev.off()

}
