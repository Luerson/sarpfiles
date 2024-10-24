#ifndef MODELNODE_H_INCLUDED
#define MODELNODE_H_INCLUDED

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
#include "mipnode.h"
#include "SarpADS.h"

using namespace std;

void initArcs (instanceStat *inst, nodeArcsStruct *nas);
void feasibleArcs (instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);
void viewSol (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, solStats *sStat);
void nodeMethod (nodeStat *node, instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, probStat* problem, solStats *sStat);
void fipnodeMethod (nodeStat *node, instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, probStat* problem, solStats *sStat);

void output(instanceStat *inst, vector<nodeStat> &nodeVec, solStats *sStat, probStat* problem);
void fipArcs(instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist, int stage);
void fipMethod(nodeStat *node, instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, probStat*problem, solStats *sStat);

void addArcsToParcelPickup (int i, instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);
void addArcsToCustomerPickup (int i, instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist); 
void addArcsToDummyFromDepot (int i, instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);
void addArcsToCustomerDelivery (int i, instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);

void removeArcsToTheSameNode (instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);
void removeDeliveryToItsPickup (instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);

void fillInfoDepotToDummy (instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);
void fillInfoToDummy (instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);
void fillInfoFromDepot (instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);
void fillInfoRequests (instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist);

void obligueDirectCustomer (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, const double **mdist);
void limitParcelCapacity (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist);
void limitCustomerCapacity (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist);

#endif