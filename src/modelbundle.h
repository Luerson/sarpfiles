#ifndef MODELBUNDLE_H_INCLUDED
#define MODELBUNDLE_H_INCLUDED

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <vector>
#include <algorithm>	
#include <iterator>
#include <math.h>
#include <cmath>
#include <limits>
#include <float.h>
#include <iomanip>
#include <ctime>
#include <ilcplex/ilocplex.h>
#include <stdlib.h>
#include <iostream>
#include <locale.h>
#include <sys/time.h>
#include <ctime>
#include <unistd.h>
#include "SarpADS.h"
#include "mipbundle.h"
#include "bundleData.h"

using namespace std;

//generate bundles and organize them into clusters
void makeBundles (instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat, clSt *cStat, vector< vector<bParcelStruct> > &clsParcel, probStat* problem);
//calculate bundles profits (for bundles of size 1 and 2, the profit is the one from the first node)
void bundleProfit(instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat);
void initVecs (instanceStat *inst, vector< vector<bParcelStruct> > &clsParcel, bundleStat *bStat, probStat* problem);
void initArcs (instanceStat *inst, bundleStat *bStat, clSt *cStat);
void feasibleBundleArcs (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat, clSt *cStat, int p, probStat* problem);
void feasibleClusterArcs (instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat, clSt *cStat, int p, probStat* problem);
void makeParcelBundles(instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat, probStat* problem);
//obtain start and end times of each bundle 
void makeStartTimes (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat, probStat* problem);
//obtain first and last elements of each bundle
void makeBundleReference (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat);
void makeSmallerProblem(instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, int p, vector< vector<bParcelStruct> > &clsParcel, probStat* problem, int Q);
bool compareCosts(const bParcelStruct &a, const bParcelStruct &b);
void makeParcelSets (instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, vector< vector<int> > &parcelSets);
void nodeSolution (instanceStat *inst, double **mdist, bundleStat *bStat, vector<nodeStat> &nodeVec, solStats *sStat);
void bundleMethod(nodeStat *node, instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, probStat* problem, solStats *sStat);
void stillTimeBundle(instanceStat *inst, double **mdist, bundleStat *bStat, vector<nodeStat> &nodeVec, solStats *sStat);


#endif