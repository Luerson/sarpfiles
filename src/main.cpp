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
#include "functions.h"
#include "mip.h"

using namespace std;

//x(hours) = Euc_dist/vmed
int main (int argc, char *argv[]) {

	double **distMatrix;
	
	nodeStat node;
	instanceStat inst;
	probStat problem;
	std::vector<nodeStat> nodeVec;
	vector< vector<bool> > arcs;
	// int vertices;
	pair<int, int> fArc;
	vector< pair<int,int> > arcVec;
	vector< pair<int,int> > nodummyarcVec;

	vector< pair<int,int> > auxVec;

	readData(argc, argv, &node, &inst, nodeVec, &distMatrix, &problem);
	
	vector< vector< pair<int,int> > > arcPlus;
	vector< vector< pair<int,int> > > arcMinus;
	vector< pair<int,int> > arcNN;

	cout << "Distance Matrix: " << endl;

	for (int i = 0; i < inst.V + inst.dummy; i++){
		for (int j = 0; j < inst.V + inst.dummy; j++){
			cout << setw(5) << distMatrix[i][j] << " ";
		}
		cout << endl;
	}
	// getchar();
	vector<bool> arcRow;
	
	for (int i = 0; i < inst.V + inst.dummy; i++){
		for (int j = 0; j < inst.V + inst.dummy; j++){
			arcRow.push_back(false);
		}
		arcs.push_back(arcRow);
	}
	// cout << "Feasible arcs: " << endl;
	
	// for (int i = 0; i < inst.V + 1; i++){
	// 	auxVec.push_back(fArc);
	// }

	for (int i = 0; i < inst.V + inst.dummy; i++){
		arcMinus.push_back(auxVec);
		arcPlus.push_back(auxVec);
	}

	// for (int i = 0; i < inst.n; i++){
	// 	arcNN.push_back(auxVec);
	// }

	feasibleArcs(&inst, nodeVec, arcs, fArc, arcPlus, arcMinus, &problem, arcNN);
	

	// cout << "\nArcsPlus: " << endl;
	// for (int i = 0; i < arcPlus.size(); i++){
	// 	cout << i << ": " << endl;
	// 	for (int j = 0; j < arcPlus[i].size(); j++){
	// 		cout << "(" << arcPlus[i][j].first << "," << arcPlus[i][j].second << ") ";
	// 	}
	// 	cout << endl;
	// }
	// getchar();
	// cout << "\nArcsMinus: " << endl;
	// for (int i = 0; i < arcMinus.size(); i++){
	// 	cout << i << ": " << endl;
	// 	for (int j = 0; j < arcMinus[i].size(); j++){
	// 		cout << "(" << arcMinus[i][j].first << ", " << arcMinus[i][j].second << ") ";
	// 	}
	// 	cout << endl;
	// }
	// for (int i = 0; i < inst.V + 1; i++){
	// 	for (int j = 0; j < arcSizePlus[i]; j++){
	// 		auxVec.push_back(fArc);
	// 	}
	// 	arcPlus.push_back(auxVec);
	// 	auxVec.clear();
	// }

	// for (int i = 0; i < nodeVec.size(); i++){
	// 	for (int j = 0; j < arcs[i].size() - inst.dummy; j++){
	// 		if (arcs[i][j]){
	// 			fArc.first = i;
	// 			fArc.first = j;
	// 			nodummyarcVec.push_back(fArc);
	// 		}
	// 	}
	// }


	// cout << "\nSizes: " << endl;
	
	// for (int i = 0; i < inst.V + 1; i++){
	// 	cout << "\ni: " << i << endl;
	// 	cout << "Plus: " << arcSizePlus[i] << endl;
	// 	cout << "Minus: " << arcSizeMinus[i] << endl;
	// }
	
	//getchar();
	
	for (int i = 0; i < nodeVec.size(); i++){
		for (int j = 0; j < nodeVec.size(); j++){
			if(arcs[i][j]){
				fArc.first = i;
				fArc.second = j;
				arcVec.push_back(fArc);
				if (j < inst.V){
					nodummyarcVec.push_back(fArc);
				}
			}
		}
	}

	// cout << "\nArcs arriving at f:" << endl;

	// for (int i = 0; i < arcMinus[inst.V].size(); i++){
	// 	cout << arcMinus[inst.V][i].first << " - " << arcMinus[inst.V][i].second << endl;
	// }
	// getchar();

	// cout << "\nNo dummy Arc Vector:" << endl;
	// for (int i = 0; i < nodummyarcVec.size(); i++){
	// 	cout << "(" << nodummyarcVec[i].first << ", " << nodummyarcVec[i].second << ") ";
	// }
	// cout << endl;

	// getchar();

	cout << "\nArc Vector:" << endl;
	for (int i = 0; i < arcVec.size(); i++){
		cout << "(" << arcVec[i].first << ", " << arcVec[i].second << ") ";
	}
	cout << endl;

	getchar();

	cout << "\nArc NN: " << endl;
	for (int i = 0; i < arcNN.size(); i++){
		cout << "(" << arcNN[i].first << ", " << arcNN[i].second << ") ";

	}
	getchar();

	for (int i = 0; i < inst.V + inst.dummy; i++){
		if (i == 0){
			cout << setw(3) << " ";
		}
		cout << setw(3) << std::right << i;
	}
	cout << endl;
	for (int i = 0; i < inst.V + inst.dummy; i++){
		cout << setw(3) << std::right << i;
		for (int j = 0; j < inst.V + inst.dummy; j++){
			cout << setw(3) << std:: right << arcs[i][j];
		}
		cout << endl;
	}

	// cout << "Instance stats: " << inst.K  << " " << inst.delta << " " <<
	// inst.n << " " << inst.m << " " << inst.T << " " << inst.V << endl;
	// cout << "\nDeltas: ";
	// for (int i =0; i < nodeVec.size(); i++){
	// 	cout << nodeVec[i].load << " ";
	// }
	// cout << endl;
	// getchar();
	

	mip(&inst, nodeVec, arcs, distMatrix, arcVec, nodummyarcVec, arcPlus, arcMinus, arcNN);

	for ( int i = 0; i < inst.V + inst.dummy; i++ ) {
		delete[] distMatrix[i];
	}

	delete[] distMatrix;


	return 0;
}