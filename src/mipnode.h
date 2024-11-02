#ifndef MIPNODE_H_INCLUDED
#define MIPNODE_H_INCLUDED

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
#include "modelnode.h"
#include "functions.h"
#include "SarpADS.h"

void mipnode(instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, probStat* problem, nodeArcsStruct *nas, solStats *sStat);
void mipnodefip(instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, probStat* problem, nodeArcsStruct *nas, solStats *sStat, fipStats *fipStat);
void printResults(instanceStat *inst, double **mdist, solStats *sStat, vector<nodeStat> &nodeVec);
void fippass(instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, probStat* problem, nodeArcsStruct *nas, solStats *sStat);
void fipmip(instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, probStat* problem, nodeArcsStruct *nas, solStats *sStat, fipStats *fipStat);
void arcBundle(instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, probStat *problem, nodeArcsStruct *nas, solStats *sStat);

/* Function constraints */
void allCustomersVisited (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x);
void sameRoutePDParcel (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x);
void sameRoutePDCustomer (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x);
void flowConservation (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x);
void conversionConstraints (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x, IloBoolVarArray &y);
void startDepot (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x);
void dummyDepot (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x);
void loadConstraints (const int bigW, const vector< int > loadVector, const nodeArcsStruct *nas, double **mdist, IloModel &model, const IloEnv &env, const IloArray <IloArray <IloBoolVarArray> > &x, const IloNumVarArray &w);
void limitCustomerDetour (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x);

/* aux Functions */
vector< int > parcelLoads(const vector<nodeStat> &nodeVec);
vector< int > customerLoads(const vector<nodeStat> &nodeVec);

void tieServiceTimeToVisit (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloBoolVarArray &y, IloNumVarArray &b);
void arcTimeOrder (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloArray <IloArray <IloBoolVarArray> > &x, IloNumVarArray &b);
void orderPD (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloNumVarArray &b);
void earlyAndLate (const instanceStat *inst, nodeArcsStruct *nas, const probStat* problem, const vector<nodeStat> &nodeVec, double **mdist, IloModel &model, IloEnv &env, IloBoolVarArray &y, IloNumVarArray &b);

#endif