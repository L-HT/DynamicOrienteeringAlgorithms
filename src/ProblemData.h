#ifndef PROBLEMDATA_H
#define PROBLEMDATA_H

#include <Rcpp.h>
#include <list_graph.h>
#include <dim2.h>

#include "DynamicChange.h"
//#include "HelperFunctions.h"


using namespace lemon;

typedef ListDigraph MyGraph;

// struct Change;

struct NodeData{
    int id_;

    double value_;
    int type_;
    dim2::Point<double> coordinates_;

    NodeData(int id, double x, double y, double value, int type)
        : id_(id), value_(value), type_(type){
        coordinates_.x = x;
        coordinates_.y = y;
    }
    NodeData(){
    }
    friend std::ostream& operator<<(std::ostream& os, const NodeData& nd){
        os << nd.id_ << "," << nd.coordinates_.x << "," << nd.coordinates_.y << "," << nd.type_ << "," << nd.value_;
        return os;
    }
    // friend std::istream& operator>>(std::istream& is, NodeData& nd){
    //     // is >> nd.id_;
    //     is >> nd.coordinates_.x;
    //
    //     is >> nd.coordinates_.y;
    //     is >> nd.type_;
    //     is >> nd.value_;
    //     return is;
    // }
};

struct NodeDataOld{
    dim2::Point<int> coordinates_;
    int divergence_;
    bool isTrueSupply_;

    NodeDataOld()
        : divergence_(0), isTrueSupply_(false){
    }
};

struct NodeRouteData{
    dim2::Point<int> coordinates_;

    // 0: normal node, 1: node with value, 2: depot
    int type_;
    int value_;

    NodeRouteData()
        : type_(0), value_(0){
    }
};

struct ArcData{
    int fromID_;
    int toID_;
    double length_;

    ArcData(int from, int to, double length)
        : fromID_(from), toID_(to), length_(length){
    }
    ArcData(){

    }
    ArcData (const ArcData& other){
        fromID_ = other.fromID_;
        toID_ = other.toID_;
        length_ = other.length_;
    }
    ArcData& operator=(const ArcData& other){
        if (this != &other){
            fromID_ = other.fromID_;
            toID_ = other.toID_;
            length_ = other.length_;
        }
        return *this;
    }
    friend std::ostream& operator<<(std::ostream& os, const ArcData& ad){
        os << ad.fromID_ << "," << ad.toID_ << "," << ad.length_;
        return os;
    }

    // friend std::istream& operator>>(std::istream& is, ArcData& ad){
    //
    //     is >> ad.fromID_;
    //     is >> ad.toID_;
    //     is >> ad.length_;
    //
    //     return is;
    // }
};

struct ArcDataOld{
    int capacity_;
    int flow_;

    ArcDataOld()
        : capacity_(0), flow_(0){
    }
};

struct ArcRouteData{

    int row1_, col1_, row2_, col2_;
    bool exists_;
    int isOnRouteOf_;

    ArcRouteData()
        : row1_(-1), col1_(-1), row2_(-1), col2_(-1), exists_(true), isOnRouteOf_(0){
    }
};


struct ProblemData{

    MyGraph& graph_;
    MyGraph::NodeMap<NodeData>& nodeMap_;
    MyGraph::ArcMap<ArcData>& arcMap_;
    MyGraph::Node& startNode_;
    std::vector<MyGraph::Node> destinations_;
    int budget_;
    std::string problemName_;
    std::string terminationCriterion_;
    std::vector<Change> dynamicChangesByEvaluation_;
    std::vector<Change> dynamicChangesByTime_;
    std::vector<Change> dynamicChangesByTestedBitVectors_;
    std::vector<Change> dynamicChangesByDistanceEvaluation_;

    ProblemData(MyGraph& graph,
                MyGraph::NodeMap<NodeData>& nodeMap,
                MyGraph::ArcMap<ArcData>& arcMap,
                MyGraph::Node& startNode,
                std::vector<MyGraph::Node> destinations,
                int budget)
        :   graph_(graph), nodeMap_(nodeMap),
            arcMap_(arcMap),
            startNode_(startNode),
            destinations_(destinations), budget_(budget), problemName_(""),
            terminationCriterion_("absoluteTime"){

    }
    ProblemData(MyGraph& graph,
                MyGraph::NodeMap<NodeData>& nodeMap,
                MyGraph::ArcMap<ArcData>& arcMap,
                MyGraph::Node& startNode,
                std::vector<MyGraph::Node> destinations,
                int budget,
                std::string problemName,
                std::string terminationCriterion_ = "absoluteTime")
        :   graph_(graph), nodeMap_(nodeMap),
            arcMap_(arcMap),
            startNode_(startNode),
            destinations_(destinations),
            budget_(budget),
            problemName_(problemName),
            terminationCriterion_("absoluteTime"){

    }
    ProblemData(const ProblemData& pd)
        :   graph_(pd.graph_), nodeMap_(pd.nodeMap_),
            arcMap_(pd.arcMap_),
            startNode_(pd.startNode_), destinations_(pd.destinations_), budget_(pd.budget_),
            problemName_(pd.problemName_),
            terminationCriterion_(pd.terminationCriterion_){

    }

};

void setStartNode(MyGraph::Node& startNode, const MyGraph& graph, const MyGraph::NodeMap<NodeData>& nodeMap, const Rcpp::DataFrame& problemDf);
void fillDestinations(std::vector<MyGraph::Node>& destinations, const MyGraph& graph, const MyGraph::NodeMap<NodeData>& nodeMap, const Rcpp::DataFrame& problemDf);

void readDynamicChangesFromFile(ProblemData& problemData, std::string path);

#endif
