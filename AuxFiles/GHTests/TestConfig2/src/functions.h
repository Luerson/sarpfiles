#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

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
#include "readdata.h"

using namespace std;

struct solStats{
	double solprofit;
	double tParcel;
	double tPass;
	double tBoth;
	double tNone;

	double tStillP;
	double tStillG;
	double tStill;

	double dParcel;
	double dPass;
	double dBoth;
	double dNone;

	double time;

	int servedParcels;
	
    vector< vector<int> > solOrder;
	vector< vector<int> > solInNode;
	vector< vector< pair<int, int> > > solvec;

	vector<double> solBegin;
	vector<double> solLoad;

	bool feasible;

	double pProfit;
	double costs;

};

void solStatIni(solStats *sStat);
void mipSolStats (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, solStats *sStat);
void printStats(instanceStat *inst, solStats *sStat);
void distScale(instanceStat *inst, int *instV, vector <vector <double> > &tempData, double *curAvg, int *scale);

#endif