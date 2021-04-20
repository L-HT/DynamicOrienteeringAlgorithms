#ifndef IO_H
#define IO_H

#include <Rcpp.h>
#include <string>
#include <sstream>
#include <lemon/grid_graph.h>
#include <lemon/list_graph.h>
#include <stdlib.h>
#include <lemon/lgf_writer.h>
#include <lemon/lgf_reader.h>

#include "GraphToList.h"
#include "HelperFunctions.h"
#include "ProblemData.h"

#include "DynamicChange.h"

using namespace Rcpp;
using namespace lemon;

/*
* Instanz einlesen und damit eine Permutation simulieren
*/

ProblemData problemDataFromLGF(std::string filePath, MyGraph& graph, MyGraph::ArcMap<ArcData>& arcMap,
                               MyGraph::NodeMap<NodeData>& nodeMap, MyGraph::Node& startNode,
                               std::vector<MyGraph::Node>& destinations);


void readGridGraphFile(std::string filePath, GridGraph& graph,
                       GridGraph::NodeMap<NodeRouteData>& nodeMap,
                       GridGraph::ArcMap<ArcRouteData>& arcMap,
                       GridGraph::NodeMap<bool>& forbiddenNodes,
                       GridGraph::ArcMap<int>& forbiddenEdges,
                       GridGraph::Node& startNode);

// My first functor ever, yay!
struct nodeArcDataToStringConverter{
    nodeArcDataToStringConverter(){
    }

    std::string operator()(NodeData nd){
        std::stringstream ss;
        ss << nd; //nd.id_ << "," << nd.coordinates_.x << "," << nd.coordinates_.y << "," << nd.type_ << "," << nd.value_;
        return ss.str();
    }

    std::string operator()(ArcData ad){
        std::stringstream ss;
        ss << ad; //ad.length_;
        return ss.str();
    }
};

struct stringToNodeDataConverter{
    stringToNodeDataConverter(){
    }

    NodeData operator()(std::string s){

        NodeData result;
        std::string delimiter = ",";
        std::size_t startSearchPos = 0;
        std::size_t findPos = 0;
        std::string token;

        findPos = s.find(delimiter, startSearchPos);
        token = s.substr(startSearchPos, findPos);
        result.id_ = std::stoi(token);
        startSearchPos = findPos + 1;

        findPos = s.find(delimiter, startSearchPos);
        token = s.substr(startSearchPos, findPos);
        result.coordinates_.x = std::stod(token);
        startSearchPos = findPos + 1;

        findPos = s.find(delimiter, startSearchPos);
        token = s.substr(startSearchPos, findPos);
        result.coordinates_.y = std::stod(token);
        startSearchPos = findPos + 1;

        findPos = s.find(delimiter, startSearchPos);
        token = s.substr(startSearchPos, findPos);
        result.type_ = std::stoi(token);
        startSearchPos = findPos + 1;

        findPos = s.find(delimiter, startSearchPos);
        token = s.substr(startSearchPos, findPos);
        result.value_ = std::stod(token);

        // Rcpp::Rcout << result;
        return result;
    }
};

struct stringToArcDataConverter{
    stringToArcDataConverter(){
    }

    ArcData operator()(std::string s){
        std::string delimiter = ",";
        ArcData result;
        std::size_t startSearchPos = 0;
        std::size_t findPos = 0;
        std::string token;

        findPos = s.find(delimiter, startSearchPos);
        token = s.substr(startSearchPos, findPos);
        result.fromID_ = std::stoi(token);
        startSearchPos = findPos + 1;

        findPos = s.find(delimiter, startSearchPos);
        token = s.substr(startSearchPos, findPos);
        result.toID_ = std::stoi(token);
        startSearchPos = findPos + 1;

        findPos = s.find(delimiter, startSearchPos);
        token = s.substr(startSearchPos, findPos);
        result.length_ = std::stod(token);

        // Rcpp::Rcout << result;
        return result;
    }
};

struct stringToChangeConverter{
    stringToChangeConverter(){
    }

    // Zeitart:
    // changeID
    // absoluteTime
    // Evaluation
    // testedBitVector
    // * Art der Änderung: relativ, absolut
    // * whatIsAffected: node/arc/budget
    // * affectedID1: id steht drin, budget: -1
    // * affectedID2: id steht drin, budget: -1
    // * change: ein double-Wert

    Change operator()(std::string s){
        char delimiter = ',';
        Change result;

        std::stringstream ss(s);
        std::vector<std::string> splitResult;
        while(ss.good()){
            std::string part;
            std::getline(ss, part, delimiter);
            splitResult.push_back(part);
        }

        result.id_ = std::stoi(splitResult[0]);
        result.atTime_ = std::stod(splitResult[1]);
        result.atEvaluation_ = std::stoi(splitResult[2]);
        result.atTestedBitVector_ = std::stoi(splitResult[3]);
        result.atDistanceEvaluation_ = std::stoi(splitResult[4]);
        result.isRelativeChange_ = (splitResult[5] == "relative");
        result.nodeID1_ = std::stoi(splitResult[7]);
        result.nodeID2_ = std::stoi(splitResult[8]);
        result.change_ = std::stod(splitResult[9]);

        // std::size_t startSearchPos = 0;
        // std::size_t findPos = 0;
        // std::string token;
        //
        // findPos = s.find(delimiter, startSearchPos);
        // token = s.substr(startSearchPos, findPos);
        // Rcpp::Rcout << "token: " << token << std::endl;
        // result.id_ = std::stoi(token);
        // startSearchPos = findPos + 1;
        //
        //
        // findPos = s.find(delimiter, startSearchPos);
        // token = s.substr(startSearchPos, findPos);
        // Rcpp::Rcout << "token: " << token << std::endl;
        // result.atTime_ = std::stod(token);
        // startSearchPos = findPos + 1;
        //
        //
        // findPos = s.find(delimiter, startSearchPos);
        // token = s.substr(startSearchPos, findPos);
        // Rcpp::Rcout << "token: " << token << std::endl;
        // result.atEvaluation_ = std::stod(token);
        // startSearchPos = findPos + 1;
        //
        //
        // findPos = s.find(delimiter, startSearchPos);
        // token = s.substr(startSearchPos, findPos);
        // Rcpp::Rcout << "token: " << token << std::endl;
        // result.isRelativeChange_ = (token == "relative");
        // startSearchPos = findPos + 1;
        //
        // findPos = s.find(delimiter, startSearchPos);
        // token = s.substr(startSearchPos, findPos);
        // //result.isRelativeChange_ = (token == "TRUE");
        // startSearchPos = findPos + 1;
        //
        // findPos = s.find(delimiter, startSearchPos);
        // token = s.substr(startSearchPos, findPos);
        // result.nodeID1_ = std::stoi(token);
        // startSearchPos = findPos + 1;
        //
        // findPos = s.find(delimiter, startSearchPos);
        // token = s.substr(startSearchPos, findPos);
        // result.nodeID2_ = std::stoi(token);
        // startSearchPos = findPos + 1;
        //
        // findPos = s.find(delimiter, startSearchPos);
        // token = s.substr(startSearchPos, findPos);
        // result.change_ = std::stod(token);


        return result;
    }
};

void saveDigraphAsLGF(std::string filePath,
               Rcpp::DataFrame& nodeDf,
               Rcpp::DataFrame& edgeDf);


void saveDigraphWithProblemDataAsLGF(std::string filePath, Rcpp::DataFrame& nodeDf, Rcpp::DataFrame& arcDf, Rcpp::DataFrame& problemDf, double budget);


// void saveDigraphWithProblemDataAsLGF(std::string filePath,
//                                      Rcpp::DataFrame& nodeDf,
//                                      Rcpp::DataFrame& arcDf,
//                                      Rcpp::DataFrame& problemDf);


void readMyGraphIntoVariables(const Rcpp::DataFrame& nodeDf,
                                  const Rcpp::DataFrame& arcDf,
                                  MyGraph& graph,
                                  MyGraph::NodeMap<NodeData>& nodeMap,
                                  MyGraph::ArcMap<ArcData>& arcMap);

void readDistanceMatrixIntoVariables(const Rcpp::DataFrame& nodeDf,
                                     const Rcpp::DataFrame& arcDf,
                                     const Rcpp::DataFrame& problemDf,
                                     MyGraph& graph,
                                     MyGraph::NodeMap<NodeData>& nodeMap,
                                     MyGraph::ArcMap<ArcData>& arcMap,
                                     const Rcpp::NumericMatrix& distanceMatrix);

void readProblemDataIntoMaps(const Rcpp::DataFrame& problemDf,
                             MyGraph& graph,
                             MyGraph::NodeMap<NodeData>& nodeMap,
                             MyGraph::ArcMap<ArcData>& arcMap);

// MyGraph::Node getNodeFromInternalID(const MyGraph& graph, const MyGraph::NodeMap<NodeData>& nodeMap, const int& id);

#endif
