#ifndef READDATA_H_INCLUDED
#define READDATA_H_INCLUDED

#include <iostream>
#include <fstream>
#include <stdlib.h>
// #include <string.h>
#include <cstring>
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

using namespace std;

void readData (int argc, char** argv, nodeStat *node, instanceStat *inst, vector<nodeStat> &nodeVec, double ***Mdist, probStat* problem);
void calcDistCsarp(double **dist, int full, int V, const vector<double> &vxs, const vector<double> &vys, const vector<double> &vxf, const vector<double> &vyf, string instType);
void calcDistGhsarp(double **dist, int full, int V, const vector<double> &vxs, const vector<double> &vys, const vector<double> &vxf, const vector<double> &vyf, string instType);
void calcDistSfsarp(double **dist, int full, int V, const vector<double> &vxs, const vector<double> &vys, const vector<double> &vxf, const vector<double> &vyf, string instType);
void tightWindowDETOUR1(double **dist, int n, int m, vector<double> &ve, vector<double> &vl, double kmPerMin, string instModel);

#endif
