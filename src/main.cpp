#include <ilcplex/ilocplex.h>     
#include <cstdio>
#include <iostream>
#include <iomanip> 
#include <vector>
#include <stack>
#include <algorithm>
#include <list>
#include <utility>
#include <cmath>
#include <ctime>
#include <sys/time.h>
#include "readdata.h"
#include "functions.h"
#include "modelnode.h"
#include "modelbundle.h"
#include "modelbundle2.h"
#include "modeltwostage.h"
#include "mipnode.h"
#include "hbundle.h"
#include "sarpILS.h"
#include "sarpConstruction.h"
#include "cpuTime.h"
#include "Statistics.h"
#include "SarpADS.h"
#include "osarp.h"

using namespace std;

int main (int argc, char *argv[]) {
	double **distMatrix;
	
	int trialK = 1;
	double trialMulti = 1.5;

	nodeStat node;
	instanceStat inst;
	probStat problem;
	solStats sStat;
	vector<nodeStat> nodeVec;
	
	sStat.feasible = false;
	
	nodeVec.clear();

	readData(argc, argv, &node, &inst, nodeVec, &distMatrix, &problem); 

	solveselect(&node, &inst, distMatrix, nodeVec, &problem, &sStat);

	return 0;
}
