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

using namespace std;

struct nodeStat{

	double xs;
	double ys;
	char label;
	double load;
	double e;
	double l;
	double xf;
	double yf;
    double delta;
    double profit;
    double trip;
    int index;
};

struct instanceStat{

    int n;
    int m;
    int K;
    double T = 24;
    int V;
    double maxTime = 8;
    int dummy;
    double service;

	double minpas = 3.24; //initial fare for passengers
	double paskm = 1.03; //fare per km for passengers

	// double minpas = 5; //initial fare for passengers
	// double paskm = 1.03; //fare per km for passengers

   	double minpar = 2.74; //initial fare for parcels
   	// double minpar = 4;
	double parkm = 0.83; //fare per km for parcels

	double costkm = 0.46;
	// double vmed;
	// double vmed = 19.3;
	// double vmed = 19.3;
	double vmed = 41;
	bool preInst;
	string InstName;
	string instType;

	vector<int> instParam;
	double totalCustomProfit;

	double realAvg = 16;
	double realStddv = 1.03;

	bool min;
};

struct probStat{

	bool sim;
	bool seq;
	string scen;
	string model;
};

void readData (int argc, char** argv, nodeStat *node, instanceStat *inst, vector<nodeStat> &nodeVec, double ***Mdist, probStat* problem, int trialK, double trialMulti);

#endif
