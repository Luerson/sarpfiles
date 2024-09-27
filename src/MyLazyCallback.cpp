
#include "MyLazyCallback.h"
#include <map>
#include <list>

/********************************************** Class' Constructor **********************************************/
MyLazyCallback::MyLazyCallback(IloEnv env, const IloArray <IloArray <IloBoolVarArray> >& x_ref, bundleStat *bStat, double **mdist, instanceStat *inst, probStat* problem, vector<nodeStat> &nodeVec, int nodes, int veic, int parcels, int customers) : IloCplex::LazyConstraintCallbackI(env), x(x_ref), x_vars(env), bStat(bStat), mdist(mdist), inst(inst), problem(problem), nodeVec(nodeVec), n(nodes), v(veic), nParcels(parcels), nCustomers(customers)
{
	int num = 0;

	int l = 0;
	for (int i = 0; i < bStat->bundleVec.size(); i++) {
        for (int j = 0; j < bStat->bundleVec.size(); j++) {
			if (!bStat->bArcs[i][j]) {
				continue;
			}
			
			for (int k1 = 0; k1 < bStat->arcV[i][j].size(); k1++)
			{
				int k = bStat->arcV[i][j][k1];
				x_vars.add(x[i][j][k]);
			}
        }
    }
}
/*****************************************************************************************************************/

/************************** Return a callback copy. This method is a CPLEX requirement ***************************/
IloCplex::CallbackI* MyLazyCallback::duplicateCallback() const 
{ 
   return new (getEnv()) MyLazyCallback(getEnv(), x, bStat, mdist, inst, problem, nodeVec, n, v, nParcels, nCustomers); 
}
/*****************************************************************************************************************/

/************************************ Callback's code that is runned by CPLEX ************************************/
void MyLazyCallback::main()
{	
	// cout << "LAZY CALLBACK" << endl;

	cout << "teste 1" << endl;

	IloNumArray x_vals(getEnv());
    getValues(x_vals, x_vars);

	const int bundlesPerCustomer = 1 + 3*this->nParcels;

	int relevantNodes = bStat->bundleVec.size() - (2*inst->K + inst->m); 

	// Matrix to store every route in the current solution
	vector<vector<int>> routes;
	map<int, int> arcDirection;
	map<pair<int, int>, int> arcChosen;
	vector<bool> used(bStat->bundleVec.size(), false);
	vector<vector<int>> subtours;

    // Imprimir os valores das variáveis relaxadas
    int firstIndex 	= 0;
	int lastIndex 	= bundlesPerCustomer*this->nCustomers + this->nParcels + 2*this->v - 1;
	int l = 0;

    for (int i = 0; i < bStat->bundleVec.size(); i++) {
        for (int j = 0; j < bStat->bundleVec.size(); j++) {

			if (!bStat->bArcs[i][j]) {
				continue;
			}

			for (int k1 = 0; k1 < bStat->arcV[i][j].size(); k1++)
			{
				int k = bStat->arcV[i][j][k1];

				if (x_vals[l++] > EPSILON)
				{
					arcDirection[i] = j;
				}
			}
        }
    }

    // Liberar a memória
    x_vals.end();

	for (int i = 0; i < this->v; i++)
	{
		int startDepotIndex = bundlesPerCustomer*this->nCustomers + this->nParcels + i;	// Start depot for the current vehicle
		int finalDepotIndex = startDepotIndex + this->v;								// Final depot for the current vehicle

		routes.emplace_back(); // Creating a new route

		routes[i].push_back(startDepotIndex);
		used[startDepotIndex] = true;

		if (arcDirection.find(startDepotIndex) != arcDirection.end())
		{
			int prevNode = startDepotIndex;
			for (int nextNode = arcDirection[startDepotIndex]; nextNode != finalDepotIndex; nextNode = arcDirection[nextNode])
			{
				routes[i].push_back(nextNode);
				prevNode = nextNode;
				used[nextNode] = true;
			}
		}
		routes[i].push_back(finalDepotIndex);
		used[finalDepotIndex] = true;
	}
	
	for (int i = 0; i < bStat->bundleVec.size(); i++)
	{
		if (!used[i] && arcDirection.find(i) != arcDirection.end()) {
			int firstElement = i;
			int aux = arcDirection[firstElement];
			vector<int> subtour;

			used[i] = true;
			subtour.push_back(i);

			while (aux != firstElement) {
				used[aux] = true;
				subtour.push_back(aux);
				aux = arcDirection[aux];
			}

			subtours.push_back(subtour);
		}
	}

	// cout << "depois" << endl;

	// for (int k = 0; k < subtours.size(); k++) {
	// 	for (int i = 0; i < subtours[k].size(); i++) {
	// 		cout << subtours[k][i] << " ";
	// 	}
	// 	cout << endl;
	// }
	// cout << "aqui 1" << endl;
	// getchar();
	// cout << "FIM ROTAS" << endl;

	/********** Printing the routes **********/
	// cout << endl;
	// for (int i = 0; i < routes.size(); i++)
	// {
	// 	 cout << "Route " << i << ": ";
	// 	for (int j = 0; j < routes[i].size()-1; j++)
	// 	{
	// 		 cout << routes[i][j] << " -> ";
	// 	}
	// 	 cout << routes[i].back() << endl;
	// }

	vector<bool> pickupBundles(bStat->bundleVec.size(), false);
	vector<bool> deliveryBundles(bStat->bundleVec.size(), false);
	vector<bool> fullParcelBundle(bStat->bundleVec.size(), false);
	vector<vector<int>> ToPrevent;

	for (int k = 0; k < routes.size(); k++) {
		for (int i = 0; i < routes[k].size(); i++) {
			int h = routes[k][i];

			int u = bStat->bundleVec[h][0];
			int v = bStat->bundleVec[h][1];

			if (bStat->bundleVec[h].size() == 2) {
				// cout << u << " " << v << endl;
				// getchar();

				if (u < inst->n && v >= inst->n + inst->m && v < inst->n + 2*inst->m) {
					deliveryBundles[h] = true;
				}
				
				else

				if (v < inst->n && u >= inst->n && u < inst->n + inst->m) {
					pickupBundles[h] = true;
				}

				else

				{
					fullParcelBundle[h] = true;
				}
			} else {
				if (u >= inst->n && u < inst->n + inst->m) {
					fullParcelBundle[h] = true;
				}
			}
		}
	}

	for (int k = 0; k < routes.size(); k++) {
		vector<int> sequence;
		bool restricted = false;

		for (int i = 0; i < routes[k].size(); i++) {
			int h = routes[k][i];

			if (deliveryBundles[h]) {
				restricted = false;
				sequence.clear();
			} else if ((pickupBundles[h] || fullParcelBundle[h]) && restricted) {
				sequence.push_back(h);
				ToPrevent.push_back(sequence);
				sequence.clear();
				restricted = false;
				i--;
			} else if (pickupBundles[h] == true && !restricted){
				sequence.push_back(h);
				restricted = true;
			} else if (restricted) {
				sequence.push_back(h);
			}
		}
	}

	for (int k = 0; k < routes.size(); k++) {
		vector<int> sequence;
		int lastRelevant = routes[k][0];
		float time = 0;

		sequence.push_back(lastRelevant);

		for (int i = 1; i < routes[k].size(); i++) {
			int u = routes[k][i-1];
			int v = routes[k][i];

			int z1 = bStat->lastElement[u];
			int z2 = bStat->firstElement[v];

			bool isRelevant = ((i == routes[k].size() - 1) || (v < relevantNodes));

			sequence.push_back(v);

			time += mdist[z1][z2]/inst->vmed;

			if (isRelevant) {

				float timeStart = (bStat->bundleStart[lastRelevant] > 9 ? bStat->bundleEnd[lastRelevant] : 9);
				float timeEnd = (bStat->bundleStart[v] > 9 ? bStat->bundleStart[v] : 19);

				if (time > timeEnd - timeStart) {
					ToPrevent.push_back(sequence);

					// if (ToPrevent.back()[0] == 64 && ToPrevent.back().back() == 48) {
					// 	cout << timeEnd << " " << timeStart << endl;
					// 	double testTime = 0;

					// 	for (int i1 = 0; i1 < (int)ToPrevent.back().size() - 1; i1++) {
					// 		for (int j1 = 0; j1 < bStat->bundleVec[ToPrevent.back()[i1]].size() - 1; j1++) {
					// 			int element = bStat->bundleVec[ToPrevent.back()[i1]][j1];
					// 			// int element2 = bStat->bundleVec[ToPrevent.back()[i1]][j1 + 1];
					// 			cout << element << " ";
					// 			// testTime += nodeVec[element].delta + mdist[element][element2];
					// 		}
					// 		int u = bStat->lastElement[ToPrevent.back()[i1]];
					// 		int v = bStat->firstElement[ToPrevent.back()[i1 + 1]];

					// 		testTime += bStat->bundleServVec[ToPrevent.back()[i1]] + mdist[u][v];
					// 	}

					// 	cout << "total time = " << testTime << endl;
					// 	getchar();
					// }
				}

				// cout << lastRelevant << " " << bStat->bundleEnd[lastRelevant] << " vs " << v << " " << bStat->bundleStart[v] << " = " << timeEnd - timeStart << endl;
				// getchar();
				
				// for (int z = 0; !ToPrevent.empty() && z < ToPrevent.back().size() - 1; z++) {
				// 	// for (int j1 = 0; j1 < bStat->bundleVec[ToPrevent.back()[i1]].size() - 1; j1++) {
				// 	// 	int element = bStat->bundleVec[ToPrevent.back()[i1]][j1];
				// 	// 	// int element2 = bStat->bundleVec[ToPrevent.back()[i1]][j1 + 1];
				// 	// 	cout << element << " ";
				// 	// 	// testTime += nodeVec[element].delta + mdist[element][element2];
				// 	// }
				// 	// int u = bStat->lastElement[ToPrevent.back()[i1]];
				// 	// int v = bStat->firstElement[ToPrevent.back()[i1 + 1]];

				// 	// testTime += bStat->bundleServVec[ToPrevent.back()[i1]] + mdist[u][v];
				// 	cout << z << " ";
				// 	getchar();
				// }
				// cout << endl;

				sequence.clear();
				time = 0;
				sequence.push_back(v);
				lastRelevant = v;
			} else {
				time += bStat->bundleServVec[v];
			}
		}
	}

	// for (int k = 0; k < ToPrevent.size(); k++) {

	// 	for (int i = 0; i < ToPrevent[k].size() - 1; i++) {
	// 		int u = ToPrevent[k][i];
	// 		int v = ToPrevent[k][i+1];
	// 		int h = u;

	// 		for (int j = 0; j < bStat->bundleVec[h].size(); j++) {
	// 			cout << bStat->bundleVec[h][j] << " ";
	// 		}
	// 	}
	// 	cout << endl;
	// }
	// cout << "aqui" << endl;
	// getchar();

	// cout << "RESTRICTED" << endl;

	// for (int k = 0; k < routes.size(); k++) {
	// 	for (int i = 0; i < routes[k].size(); i++) {
	// 		int h = routes[k][i];

	// 		for (int j = 0; j < bStat->bundleVec[h].size(); j++) {
	// 			cout << bStat->bundleVec[h][j] << " ";
	// 		}
	// 	}
	// 	cout << endl;
	// }

	// cout << "FIM ROTAS" << endl;

	bool lazyCut = false;

	for (int k = 0; k < ToPrevent.size(); k++) {
		IloExpr expr(getEnv());

		for (int i = 0; i < ToPrevent[k].size() - 1; i++) {
			int u = ToPrevent[k][i];
			int v = ToPrevent[k][i+1];
			int h = u;

			cout << u << " ";

			for (int k1 = 0; k1 < bStat->arcV[u][v].size(); k1++) {
				int k2 = bStat->arcV[u][v][k1];

				expr += x[u][v][k2];
			}

			// for (int j = 0; j < bStat->bundleVec[h].size(); j++) {
			// 	cout << bStat->bundleVec[h][j] << " ";
			// }
		}

		cout << ToPrevent[k].back() << endl;
		// for (int j = 0; j < ToPrevent.back().size(); j++) {
		// 	cout << bStat->bundleVec.back()[j] << " ";
		// }
		// cout << endl;
		// getchar();

		lazyCut = true;

		add(expr <= (int)ToPrevent[k].size() - 2);
	}

	for (int k = 0; k < subtours.size(); k++) {
		IloExpr expr(getEnv());

		for (int i = 0; i < subtours[k].size(); i++) {
			for (int j = 0; j < subtours[k].size(); j++) {
				
				int u  = subtours[k][i];
				int v = subtours[k][j];
				
				for (int k1 = 0; k1 < bStat->arcV[u][v].size(); k1++) {
					int k2 = bStat->arcV[u][v][k1];

					expr += x[u][v][k2];
				}
			}

			// for (int j = 0; j < bStat->bundleVec[h].size(); j++) {
			// 	cout << bStat->bundleVec[h][j] << " ";
			// }
		}
		// cout << endl;
		lazyCut = true;

		add(expr <= (int)subtours[k].size() - 1);
	}

	// if (lazyCut) {
	// 	cout << "LAZYCUT" << endl;
	// }

	cout << "teste 2" << endl;

	/********** Searching for infeasibility **********/
}
/*****************************************************************************************************************/


// bool MyLazyCallback::isCustomer(int index) {
// 	return (index >= 0) && (index < this->nCustomers);
// }


// bool MyLazyCallback::isDepot(int index) {
// 	return (index >= this->nCustomers + 2*this->nParcels) && (index < this->nCustomers + 2*this->nParcels + 2*this->v);
// }

// bool MyLazyCallback::isParcel(int index) {
// 	return (index >= this->nCustomers) && (index < (this->nCustomers + 2*this->nParcels));
// }

// bool MyLazyCallback::isPickUpParcel(int index) {
// 	return (index >= this->nCustomers) && (index < (this->nCustomers + this->nParcels));
// }

// bool MyLazyCallback::isDeliveryParcel(int index) {
// 	return (index >= this->nCustomers + this->nParcels) && (index < (this->nCustomers + 2*this->nParcels));
// }

// bool MyLazyCallback::finalVerifier(vector<int> &route) {
// 	int contP = 0;
// 	int contD = 0;

// 	for (int i = 0; i < route.size(); i++) {
// 		if (isPickUpParcel(route[i])) {
// 			contP--;
// 		} else if (isDeliveryParcel(route[i])) {
// 			contD--;
// 		} else if (isCustomer(route[i])){
// 			contP = min(1, contP + 1);
// 			contD = min(1, contD + 1);
// 		}

// 		if (contP < -1 || contD < -1) {
// 			return false;
// 		}
// 	}

// 	if (contP == -1 || contD == -1) {
// 		return false;
// 	}

// 	return true;
// }

// #include "MyLazyCallback.h"
// #include <map>

// /********************************************** Class' Constructor **********************************************/
// MyLazyCallback::MyLazyCallback(IloEnv env, const IloArray<IloArray<IloBoolVarArray>>& x_ref, nodeArcsStruct *nas, int nodes, int veic, int parcels, int customers) : IloCplex::LazyConstraintCallbackI(env), x(x_ref), x_vars(env), n(nodes), v(veic), arrConj(nas), nParcels(parcels), nCustomers(customers)
// {
// 	int num = 0;
// 	/********** Filling x_vars **********/
// 	for(int i = 0; i < n; i++) {
// 		for(int j = 0; j < n; j++){
// 			if (!nas->arcs[i][j]) {
// 				continue;
// 			}
// 			for (int k1 = 0; k1 < nas->arcV[i][j].size(); k1++) {
// 				int k = nas->arcV[i][j][k1];
// 				x_vars.add(x[i][j][k]);
// 			}
// 		}
// 	}
// 	/************************************/
// }
// /*****************************************************************************************************************/

// /************************** Return a callback copy. This method is a CPLEX requirement ***************************/
// IloCplex::CallbackI* MyLazyCallback::duplicateCallback() const 
// { 
//    return new (getEnv()) MyLazyCallback(getEnv(), x, arrConj, n, v, nParcels, nCustomers); 
// }
// /*****************************************************************************************************************/

// /************************************ Callback's code that is runned by CPLEX ************************************/
// void MyLazyCallback::main()
// {	
// 	// /********** Getting the relaxed variables values **********/
// 	// IloNumArray x_vals(getEnv(), (0.5*(n)*(n-1)));
// 	// getValues(x_vals, x_vars);
// 	// /**********************************************************/
   
// 	// vector <vector<int> > cutSetPool;
// 	// vector<IloConstraint> cons; 

// 	IloNumArray x_vals(getEnv());
//     getValues(x_vals, x_vars);

// 	// Matrix to store every route in the current solution
// 	vector<vector<int>> routes;
// 	map<int, int> arcDirection;

// 	bool restricted = false;

//     // Imprimir os valores das variáveis relaxadas
//     int index = 0;
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
// 			if (!arrConj->arcs[i][j]) {
// 				// // TODO UNCOMMENT // cout << i << " " << j << endl;
// 				continue;
// 			}

// 			for (int k1 = 0; k1 < arrConj->arcV[i][j].size(); k1++)
// 			{
// 				int k = arrConj->arcV[i][j][k1];
// 				if (x_vals[index] > EPSILON)
// 				{
// 					arcDirection[i] = j;
// 					// // TODO UNCOMMENT // cout << i << " " << j << " " << k << endl;
// 				}
// 				index++;
// 			}
//         }
//     }

//     // Liberar a memória
//     x_vals.end();

// 	// Remember: "this->v" is the numver of available vehicles
// 	// // TODO UNCOMMENT // cout << endl;
// 	// // TODO UNCOMMENT // cout << "ROTAS" << endl;
// 	// // TODO UNCOMMENT // cout << "numCustomers = " << this->nCustomers << endl;
// 	// // TODO UNCOMMENT // cout << "numParcels = " << this->nParcels << endl;
// 	// // TODO UNCOMMENT // cout << "numVehicles = " << this->v << endl;
// 	for (int i = 0; i < this->v; i++)
// 	{
// 		int startDepotIndex = this->nCustomers + 2*this->nParcels + i;	// Start depot for the current vehicle
// 		int finalDepotIndex = startDepotIndex + this->v;				// Final depot for the current vehicle

// 		// // TODO UNCOMMENT // cout << endl;
// 		// // TODO UNCOMMENT // cout << "i = " << i << endl;
// 		// // TODO UNCOMMENT // cout << "startDepotIndex = " << startDepotIndex << endl;
// 		// // TODO UNCOMMENT // cout << "finalDepotIndex = " << finalDepotIndex << endl;
// 		// // TODO UNCOMMENT // cout << "arcDirection[" << startDepotIndex << "] = " << arcDirection[startDepotIndex] << endl;
// 		// // TODO UNCOMMENT // cout << endl;

// 		if (arcDirection.find(startDepotIndex) != arcDirection.end())
// 		{
// 			routes.emplace_back(); // Creating a new route

// 			routes[i].push_back(startDepotIndex);
// 			for (int nextNode = arcDirection[startDepotIndex]; nextNode != finalDepotIndex; nextNode = arcDirection[nextNode])
// 			{
// 				// // TODO UNCOMMENT // cout << "nextNode = " << nextNode << endl;
// 				routes[i].push_back(nextNode);
// 			}
// 			routes[i].push_back(finalDepotIndex);
// 		}
// 	}
// 	// // TODO UNCOMMENT // cout << "FIM ROTAS" << endl;

// 	// Printing the routes
// 	// // TODO UNCOMMENT // cout << endl;
// 	// for (int i = 0; i < routes.size(); i++)
// 	// {
// 	// 	// TODO UNCOMMENT // cout << "Route " << i << ": ";
// 	// 	for (int j = 0; j < routes[i].size()-1; j++)
// 	// 	{
// 	// 		// TODO UNCOMMENT // cout << routes[i][j] << " -> ";
// 	// 	}
// 	// 	// TODO UNCOMMENT // cout << routes[i].back() << endl;
// 	// }
// 	// getchar();


// 	/**************** Creating the Lazy Constraints ****************/
// 	// // TODO UNCOMMENT // cout << "SIZES: " << x.getSize() << " " << x[0].getSize() << endl;

// 	for (int i = 0; i < routes.size(); i++)
// 	{
// 		int counter = 0;

// 		for (int j = 0; j < routes[i].size(); j++)
// 		{
// 			if (isParcel(routes[i][j]))
// 			{
// 				counter++;

// 				if (counter >= 4 && isDeliveryParcel(routes[i][j-3]) && isPickUpParcel(routes[i][j]))
// 				{
// 					vector<pair<int, int>> pairs = makePairs(routes[i], j - 3);
// 					// Prohibiting the sequence for all vehicles
// 					// // TODO UNCOMMENT // cout << endl;
// 					// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// 					IloExpr expr(getEnv());
// 					for (int vehicle = 0; vehicle < this->v; vehicle++)
// 					{
// 						// int k = arrConj->arcV[i][j][vehicle];

// 						// TODO UNCOMMENT // cout << endl;
// 						// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

// 						// // // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// 						// expr += x[ routes[i][j-5] ][ routes[i][j-4] ][ vehicle ];

// 						// // // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// 						// expr += x[ routes[i][j-4] ][ routes[i][j-3] ][ vehicle ];

// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// 						// expr += x[ routes[i][j-3] ][ routes[i][j-2] ][ vehicle ];

// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-2] << "," << routes[i][j-1] << "," << vehicle << ")" << endl;
// 						expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-1] << "," << routes[i][j] << "," << vehicle << ")" << endl;
// 						// expr += x[ routes[i][j-1] ][ routes[i][ j ] ][ vehicle ];

// 						for (int i1 = 0; i1 < pairs.size(); i1++) {
// 							for (int j1 = 0; j1 < pairs.size(); j1++) {
// 								if (i1 != j1 && arrConj->arcs[ pairs[i1].second ][ pairs[j1].first ])  {
// 									expr += x[ pairs[i1].second ][ pairs[j1].first ][ vehicle ];
// 								}
// 							}
// 						}

// 						// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// 						// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// 						// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// 						// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// 						// }
// 						// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// 						// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// 						// }
// 						// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// 						// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// 					}
// 					add(expr <= 2);
// 					// // TODO UNCOMMENT // cout << endl;
// 					// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// 				}
// 			}
// 			else // isCustomer = true
// 			{
// 				counter = 0;
// 			}
// 		}
// 	}

// 	// Prohibiting sequence V-P-D-P-D
// 	for (int i = 0; i < routes.size(); i++)
// 	{
// 		int counter = 0;

// 		for (int j = 0; j < routes[i].size() && j < 4; j++)
// 		{
// 			if (isParcel(routes[i][j]))
// 			{
// 				counter++;

// 				if (counter >= 3)
// 				{
// 					// Prohibiting the sequence for all vehicles
// 					// // TODO UNCOMMENT // cout << endl;
// 					// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// 					for (int vehicle = 0; vehicle < this->v; vehicle++)
// 					{
// 						// int k = arrConj->arcV[i][j][vehicle];
// 						IloExpr expr(getEnv());
// 						// TODO UNCOMMENT // cout << endl;
// 						// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

// 						if (arrConj->arcs[ routes[vehicle][0] ][ routes[i][j-2] ]) {
// 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// 							expr += x[ routes[vehicle][0] ][ routes[i][j-2] ][ vehicle ];

// 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// 							expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// 							expr += x[ routes[i][j-1] ][ routes[i][j] ][ vehicle ];

// 							if (arrConj->arcs[ routes[i][j]+nParcels ][ routes[i][j-2] ] ) {
// 								if (arrConj->arcs[ routes[vehicle][0] ][ routes[i][j] ]) {
// 									expr += x[ routes[i][j]+nParcels ][ routes[i][j-2] ][ vehicle ];

// 									expr += x[ routes[vehicle][0] ][ routes[i][j] ][ vehicle ];
// 								}
// 							}

// 							// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// 							// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// 							// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// 							// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// 							// }
// 							// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// 							// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// 							// }
// 							// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// 							// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// 						}
// 						add(expr <= 2);
// 					}
// 					// // TODO UNCOMMENT // cout << endl;
// 					// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// 				}
// 			}
// 			else // isCustomer = true
// 			{
// 				counter = 0;
// 			}
// 		}
// 	}

// 	for (int i = 0; i < routes.size(); i++)
// 	{
// 		int counter = 0;

// 		for (int j = routes[i].size() - 4; j < routes[i].size(); j++)
// 		{
// 			if (isParcel(routes[i][j]))
// 			{
// 				counter++;

// 				if (counter >= 3)
// 				{
// 					// Prohibiting the sequence for all vehicles
// 					// // TODO UNCOMMENT // cout << endl;
// 					// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// 					for (int vehicle = 0; vehicle < this->v; vehicle++)
// 					{
// 						// int k = arrConj->arcV[i][j][vehicle];
// 						IloExpr expr(getEnv());
// 						// TODO UNCOMMENT // cout << endl;
// 						// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

// 						if (arrConj->arcs[ routes[i][j] ][ routes[vehicle].back() ]) {
// 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// 							expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// 							expr += x[ routes[i][j-1] ][ routes[i][j] ][ vehicle ];

// 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// 							expr += x[ routes[i][j] ][ routes[vehicle].back() ][ vehicle ];

// 							if (arrConj->arcs[ routes[i][j] ][ routes[i][j-2] - nParcels ] ) {
// 								if (arrConj->arcs[ routes[i][j-2] ][ routes[vehicle].back() ]) {
// 									expr += x[ routes[i][j] ][ routes[i][j-2] - nParcels ][ vehicle ];

// 									expr += x[ routes[i][j-2] ][ routes[vehicle].back() ][ vehicle ];
// 								}
// 							}

// 							// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// 							// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// 							// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// 							// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// 							// }
// 							// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// 							// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// 							// }
// 							// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// 							// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// 						}
// 						add(expr <= 2);
// 					}
// 					// // TODO UNCOMMENT // cout << endl;
// 					// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// 				}
// 			}
// 			else // isCustomer = true
// 			{
// 				counter = 0;
// 			}
// 		}
// 	}

// 	for (int i = 0; i < routes.size(); i++)
// 	{
// 		// // TODO UNCOMMENT // cout << "aqui" << endl;
// 		if (!finalVerifier(routes[i])) {
// 			IloExpr expr(getEnv());
			
// 			// // TODO UNCOMMENT // cout << "aqui 2" << endl;
// 			for (int j = 0; j < routes[i].size() - 1; j++) {
// 				expr += x[ routes[i][j] ][ routes[i][j+1] ][ i ];
// 			}

// 			// // TODO UNCOMMENT // cout << "aqui 3" << endl;

// 			int maxArcs = routes[i].size() - 2;
// 			add(expr <= maxArcs);
// 		}
// 	}

// 	/***************************************************************/

// 	// double **x_edge = new double*[n];
 
// 	// for (int i = 0; i < n; i++) {
// 	// 	x_edge[i] = new double[n];
// 	// }

// 	// int l = 0;
// 	// for(int i = 0; i < n; i++) {
// 	// 	for(int j = i+1; j < n; j++) {
// 	// 		x_edge[i][j] = x_vals[l++];
// 	// 	}
// 	// }
	
// 	// x_vals.end();

// 	// cutSetPool = MaxBack(x_edge, n);

// 	// /***************** Creating the constraints ***************/
// 	// for (int c = 0; c < cutSetPool.size(); c++) {
// 	// 	IloExpr p(getEnv());
// 	// 	for(int i = 0; i < cutSetPool[c].size(); i++){
// 	// 		for(int j = 0; j < cutSetPool[c].size(); j++){
// 	// 			if(cutSetPool[c][i] < cutSetPool[c][j]){
// 	// 				p += x[cutSetPool[c][i]][cutSetPool[c][j]];
// 	// 			}
// 	// 		}
// 	// 	}
// 	// 	int RHS = cutSetPool[c].size();
// 	// 	cons.push_back(p <= RHS - 1);
// 	// }
// 	// /**********************************************************/

// 	// /*********** Adding the constraints to the model **********/
// 	// for(int i = 0; i < cons.size(); i++){
// 	// 	add(cons.at(i)).end();
// 	// }
// 	// /**********************************************************/

// 	// /******************* Cleaning the memory ******************/
// 	// for (int i = 0; i < n; i++) {
// 	// 	delete[] x_edge[i];
// 	// }
// 	// delete[] x_edge;
// 	/**********************************************************/
// }
// /*****************************************************************************************************************/


// bool MyLazyCallback::isCustomer(int index) {
// 	return (index >= 0) && (index < this->nCustomers);
// }


// bool MyLazyCallback::isDepot(int index) {
// 	return (index >= this->nCustomers + 2*this->nParcels) && (index < this->nCustomers + 2*this->nParcels + 2*this->v);
// }


// bool MyLazyCallback::isParcel(int index) {
// 	return (index >= this->nCustomers) && (index < (this->nCustomers + 2*this->nParcels));
// }

// bool MyLazyCallback::isPickUpParcel(int index) {
// 	return (index >= this->nCustomers) && (index < (this->nCustomers + this->nParcels));
// }

// bool MyLazyCallback::isDeliveryParcel(int index) {
// 	return (index >= this->nCustomers + this->nParcels) && (index < (this->nCustomers + 2*this->nParcels));
// }

// bool MyLazyCallback::finalVerifier(vector<int> &route) {
// 	int contP = 0;
// 	int contD = 0;

// 	for (int i = 0; i < route.size(); i++) {
// 		if (isPickUpParcel(route[i])) {
// 			contP--;
// 		} else if (isDeliveryParcel(route[i])) {
// 			contD--;
// 		} else if (isCustomer(route[i])){
// 			contP = min(1, contP + 1);
// 			contD = min(1, contD + 1);
// 		}

// 		if (contP < -1 || contD < -1) {
// 			return false;
// 		}
// 	}

// 	if (contP == -1 || contD == -1) {
// 		return false;
// 	}

// 	return true;
// }

// vector<pair<int, int>> MyLazyCallback::makePairs(vector<int> &route, int index) {
// 	vector<pair<int, int>> pares;

// 	pares.push_back(make_pair( route[index]-nParcels , route[index] ));
// 	pares.push_back(make_pair( route[index+1] , route[index+2] ));
// 	pares.push_back(make_pair( route[index+3] , route[index+3]+nParcels ));

// 	return pares;
// }