#include <Rcpp.h>
#include <cstdio>
#include <cstdlib>
#include <sstream>
#include <map>

#include "LINKERN_support_CCrandstate_cc.h"
// #include "hallowelt.h"
#include "LINKERN_support_data.h"
#include "LINKERN_linkern.h"

#include "ProblemData.h"
#include "HelperFunctions.h"
#include "Solver.h"
#include "VectorConversion.h"

/*
 * To-Do: später bei der Distanzmatrix checken, dass auch die richtigen Zeilen und Spalten genommen werden
 * (wegen des Offsets +1)
 */


// note: Lin-Kernighan only supports integer distances
/*
 * problem: dist is the distance matrix with -all- destinations, but only the
 * distances for the given solution are needed
 */
int** buildSubMatrixInC(Rcpp::NumericMatrix dist, const int* incycle, std::size_t n){

    if (dist.nrow() != dist.ncol()){
        //Rcpp::Rcerr << "buildMatrixInC: dist must be a quadratic matrix!" << std::endl;
        Rcpp::stop("buildMatrixInC: dist must be a quadratic matrix!");
    }

    int** result = (int**) std::malloc(sizeof(int*) * n);

    for (std::size_t i = 0; i < n; i++){
        result[i] = (int*) std::malloc(sizeof(int) * n);
        result[i][i] = 0;
    }

    for (std::size_t i = 0; i < n-1; i++){
        for (std::size_t j = i+1; j < n; j++){
            result[i][j] = std::round(dist(incycle[i],incycle[j]));
            result[j][i] = result[i][j]; // it is assumed that the distances are symmetric
        }
    }
    return result;
}

// // convert solution to vector with indices referring to destinations vector (containing the starting node as first element)
// std::vector<int> nodeVectorToIndexVectorWithLimit(const ProblemData& pd, const std::vector<MyGraph::Node>& solution, std::size_t upperLimit){
//     std::vector<int> result(solution.size(), 0);
//     int outerCounter = 0;
//
//     //for (MyGraph::Node n : solution){
//     for (std::size_t i = 0; i < upperLimit; i++){ // -1 because last node is also starting node
//         MyGraph::Node n = solution[i]; // maybe a pointer here?
//         bool nodeFound = false;
//         int innerCounter = 0;
//         do{
//             // Rcpp::Rcout << "innerCounter: " << innerCounter << " / " << pd.destinations_.size() << std::endl;
//             if (n == pd.destinations_[innerCounter]){
//                 // + 1 because index 0 is used for the starting node (which is not contained in pd.destinations_))
//                 result[outerCounter++] = innerCounter + 1;
//                 // Rcpp::Rcout << "write: " << innerCounter + 1 << std::endl;
//                 nodeFound = true;
//             } else {
//                 innerCounter++;
//                 if (innerCounter >= (int) pd.destinations_.size()){
//                     // if (pd.graph_.id(n) == pd.graph_.id(pd.startNode_)){
//                     if (pd.nodeMap_[n].id_ == pd.nodeMap_[pd.startNode_].id_){
//                         result[outerCounter++] = 0;
//                         // Rcpp::Rcout << "write: " << 0 << std::endl;
//                         nodeFound = true;
//                     } else {
//                         std::stringstream ss;
//                         ss << "Cannot find node ID " << pd.nodeMap_[n].id_ << "!" << std::endl;
//                         for (MyGraph::Node n2 : pd.destinations_){
//                            Rcpp::Rcout << pd.nodeMap_[n2].id_ << " ";
//                         }
//                         Rcpp::Rcout << std::endl;
//                         Rcpp::stop(ss.str());
//                     }
//                 }
//             }
//         } while (!nodeFound);
//     }
//
//     return result;
// }
//
// std::vector<int> nodeVectorToIndexVectorWithStartNode(const ProblemData& pd, const std::vector<MyGraph::Node>& solution){
//     return nodeVectorToIndexVectorWithLimit(pd, solution, solution.size() - 1);
// }
// std::vector<int> nodeVectorToIndexVector(const ProblemData& pd, const std::vector<MyGraph::Node>& solution){
//     return nodeVectorToIndexVectorWithLimit(pd, solution, solution.size());
// }
//
// std::vector<MyGraph::Node> indexVectorToNodeVector(const ProblemData& pd, const int* indices,
//                                                    const std::size_t numberOfElements,
//                                                    const std::map<int, int>& indicesMap){
//
//     std::vector<MyGraph::Node> result;
//
//     // shift the solution such that the depot ist the first node in the solution
//     std::size_t startPos = 0;
//     for (std::size_t i = 0; i < numberOfElements; i++){
//         if (indices[i] == 0){
//             startPos = i;
//         }
//     }
//
//     // Rcpp::Rcout << "startPos: " << startPos << ", " << numberOfElements << " elements." << std::endl;
//
//     result.push_back(pd.startNode_);
//     for (std::size_t i = startPos+1; i < numberOfElements; i++){
//         int correspondingIndexInDestinations = indicesMap.at(indices[i]) - 1;
//
//         result.push_back(pd.destinations_[correspondingIndexInDestinations]);
//         // Rcpp::Rcout << pd.graph_.id(pd.destinations_[correspondingIndexInDestinations]) << " - " << pd.graph_.id(result.back()) << std::endl;
//     }
//     for (std::size_t i = 0; i < startPos; i++){
//         int correspondingIndexInDestinations = indicesMap.at(indices[i]) - 1;
//         result.push_back(pd.destinations_[correspondingIndexInDestinations]);
//     }
//     result.push_back(pd.startNode_);
//     return result;
// }


// it is assumed that initialSolution contains the start node as first and last element
std::vector<MyGraph::Node> callRepeatedLinKernighan(const ProblemData& pd, const std::vector<MyGraph::Node>& initialSolution,
                              Rcpp::NumericMatrix dist, Solver& solver){

    std::vector<MyGraph::Node> result;
    AdditionalLogData& logData = solver.additionalLogData_;

    std::map<int, int> indicesMap;

    // printf("Start Lin-Kernighan heuristic!\n");
    // waitForInput("(press Enter)", true);
    // prepare RNG
    CCrandstate myRandstate;
    CCutil_sprand(getRandomNumber(1,100000), &myRandstate);

    compass_data mycd;
    compass_init_data (&mycd);
    compass_data_set_norm (&mycd, CC_MATRIXNORM);

    int N = initialSolution.size()-1;//pd.destinations_.size();
    mycd.n = N;

    // set some default values
    mycd.rhdat.dist_00=0;
    mycd.rhdat.dist_01=0;
    mycd.rhdat.dist_02=0;
    mycd.rhdat.dist_12=0;
    mycd.rhdat.dist_22=0;
    mycd.rhdat.p=0;
    mycd.rhdat.rhlength=0;
    mycd.rhdat.space=NULL;

    int ncount = N;
    int ecount = N*(N-1)/2.0; // N choose 2 (N ueber 2)

    // waitForInput("(start writing elist)", true);
    int* elist = (int*) std::malloc(sizeof(int) * 2* ecount);

    int counter = 0;
    for (int i = 0; i < N-1; i++){
        for (int j = i+1; j < N; j++){
            elist[counter++] = i;
            elist[counter++] = j;
        }
    }

    // for (int i=0; i < 2* ecount; i++){
    //     printf("%d ", elist[i]);
    // }

    int stallcount = 100000000;
    int repeatcount = 0;
    // waitForInput("(getIndexVector...)", true);

    std::vector<int> tempIndexVector = nodeVectorToIndexVectorWithStartNode(pd, initialSolution);
    int* incycle = tempIndexVector.data();// = (int*) std::malloc(sizeof(int) *  tempIndexVector.size()); //tempIndexVector.data();
    // for (std::size_t i = 0; i < tempIndexVector.size(); i++){
    //     incycle[i] = tempIndexVector[i];
    // }

    // for (int i=0; i < ncount; i++){
    //     Rcpp::Rcout << pd.graph_.id(pd.destinations_[incycle[i]-1]) << " :: " << pd.graph_.id(initialSolution[i]) << std::endl;
    // }
    // Rcpp::Rcout << std::endl;
    // printNodeIdsOfVector(pd, initialSolution);

    // waitForInput("(start buildMatrixInC)", true);
    mycd.adj = buildSubMatrixInC(dist, incycle, tempIndexVector.size());


    int* outcycle = (int*) std::malloc(sizeof(int) * ncount);
    for (int i=0; i < ncount; i++){
        indicesMap.insert(std::pair<int,int>(i, incycle[i]));
        // Rcpp::Rcout << "pair: " << i << ", " << incycle[i] << std::endl;
        incycle[i] = i;
        outcycle[i] = i;
    }
    // Rcpp::Rcout << std::endl;
    // waitForInput("input cycle printed", true);

    // print distance matrix
    // for (int i = 0; i < N; i++){
    //     for (int j = 0; j < N; j++){
    //         Rcpp::Rcout << std::setw(2) << mycd.adj[i][j] << " ";
    //     }
    //     Rcpp::Rcout << std::endl;
    // }

    double val = -1.0;

    int silent = 1;
    double time_bound = -1.0;
    double length_bound = -1.0;
    char* saveit_name = (char *) NULL;

    int kicktype = CC_LK_WALK_KICK;

    // waitForInput("(start linkern...)", true);
    CClinkern_tour(ncount, &mycd, ecount,
                   elist, stallcount, repeatcount, incycle,
                   outcycle, &val,
                   silent, time_bound, length_bound,
                   saveit_name, kicktype, &myRandstate, &logData);

    std::vector<MyGraph::Node> resultVector = indexVectorToNodeVector(pd, outcycle, ncount, indicesMap);
    result.assign(resultVector.begin(), resultVector.end());


    std::free(elist);
    std::free(outcycle);

    for (int i = 0; i < N; i++){
        std::free(mycd.adj[i]);
    }
    std::free(mycd.adj);

    // forced logging of the new solution
    // solver.evaluateSolution(solver.problemData_, result, solver.targetCriterion_, true, true);
    return result;
}

// adds the starting nodes at both ends of the vector and removes them at the end
std::vector<MyGraph::Node> callRepeatedLinKernighanWithoutStartNodes(const ProblemData& pd, const std::vector<MyGraph::Node>& initialSolution,
                                                                      Rcpp::NumericMatrix dist, Solver& solver){
    std::vector<MyGraph::Node> tempSolution;
    tempSolution.push_back(pd.startNode_);
    tempSolution.insert(tempSolution.end(), initialSolution.begin(), initialSolution.end());
    tempSolution.push_back(pd.startNode_);

    std::vector<MyGraph::Node> result = callRepeatedLinKernighan(pd, tempSolution, dist, solver);
    result.erase(result.end()-1);
    result.erase(result.begin());

    return result;
}

// //' @export
// //  [[Rcpp::export]]
// int asd(){
//     printf("Hallo Weelt!\n");
//
//     CCrandstate testRandstate;
//     //halloWelt();
//
//     //int k = blub();
//     //printf("blub: %d \n", k);
//     CCutil_sprand(1000, &testRandstate);
//
//     compass_data mycd;
//
//     compass_init_data (&mycd);
//     compass_data_set_norm (&mycd, CC_MATRIXNORM);
//
//     int N = 5;
//     mycd.n = N;
//
//     mycd.adj = (int**) std::malloc (sizeof(int*) * N);
//     for (int i = 0; i < N; i++){
//         mycd.adj[i] = (int*) std::malloc (sizeof(int) * N);
//     }
//
//     for (int i = 0; i < N; i++){
//         for (int j = 0; j < N; j++){
//             mycd.adj[i][j] = (j-i)*(j-i);
//         }
//     }
//     for (int i = 0; i < N; i++){
//         mycd.adj[i][i] = 0;
//     }
//
//     mycd.rhdat.dist_00=0;
//     mycd.rhdat.dist_01=0;
//     mycd.rhdat.dist_02=0;
//     mycd.rhdat.dist_12=0;
//     mycd.rhdat.dist_22=0;
//     mycd.rhdat.p=0;
//     mycd.rhdat.rhlength=0;
//     mycd.rhdat.space=NULL;
//
//     //mycd.adj[0][0] = 0;
//     //mycd.adj[0][1] = 1;
//     //mycd.adj[1][0] = 1;
//     //mycd.adj[1][1] = 0;
//
//
//     int ncount = N;
//     int ecount = 10;
//     int* elist = (int*) std::malloc(sizeof(int) * 2* ecount);
//
//     elist[0] = 0; elist[1] = 1;
//     elist[2] = 0; elist[3] = 4;
//     elist[4] = 0; elist[5] = 3;
//     elist[6] = 0; elist[7] = 2;
//     elist[8] = 1; elist[9] = 3;
//     elist[10] = 1; elist[11] = 4;
//     elist[12] = 1; elist[13] = 2;
//     elist[14] = 2; elist[15] = 3;
//     elist[16] = 2; elist[17] = 4;
//     elist[18] = 3; elist[19] = 4;
//     //for (int i = 0; i < N; i++){
//     //for (int j = 0; j < N; j++){
//     //elist[2*N*i + 2*j] = i;
//     //elist[2*N*i + 2*j + 1] = j;
//     //}
//     //}
//
//     for (int i=0; i < 2* ecount; i++){
//         printf("%d ", elist[i]);
//     }
//
//     int stallcount = 100000000;
//     int repeatcount = 0;
//
//     int* incycle = (int*) std::malloc(sizeof(int) * ncount);
//     int* outcycle = (int*) std::malloc(sizeof(int) * ncount);
//
//     for (int i=0; i < ncount; i++){
//         incycle[i] = (ncount-1)-i;
//         outcycle[i] = i;
//     }
//     double val = -1.0;
//
//     int silent = 0;
//     double time_bound = -1.0;
//     double length_bound = -1.0;
//     char* saveit_name = (char *) NULL;
//
//     int kicktype = CC_LK_WALK_KICK;
//
//     CClinkern_tour(ncount, &mycd, ecount,
//                    elist, stallcount, repeatcount, incycle,
//                    outcycle, &val,
//                    silent, time_bound, length_bound,
//                    saveit_name, kicktype, &testRandstate);
//
//     //CClinkern_tour (int ncount, compass_data *dat, int ecount,
//     //   int *elist, int stallcount, int repeatcount, int *incycle,
//     //    int *outcycle, double *val,
//     //    int silent, double time_bound, double length_bound,
//     //    char *saveit_name, int kicktype, CCrandstate *rstate);
//
//     std::free(elist);
//     std::free(incycle);
//     std::free(outcycle);
//     std::free(mycd.adj);
//
//     printf("Ende Welt!\n");
//     return 0;
// }
