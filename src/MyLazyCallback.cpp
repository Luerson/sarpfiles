
// #include "MyLazyCallback.h"
// #include <map>
// #include <list>

// /********************************************** Class' Constructor **********************************************/
// MyLazyCallback::MyLazyCallback(IloEnv env, const IloArray <IloArray < IloArray <IloBoolVarArray> > >& x_ref, nodeArcsStruct *arrConj, double **mdist, instanceStat *inst, probStat* problem, vector<nodeStat> &nodeVec, double* p_bestSolVal, int nodes, int veic, int parcels, int customers) : IloCplex::LazyConstraintCallbackI(env), x(x_ref), x_vars(env), mdist(mdist), inst(inst), problem(problem), nodeVec(nodeVec), bestSolVal(p_bestSolVal), n(nodes), v(veic), nas(arrConj), nParcels(parcels), nCustomers(customers)
// {
// 	int num = 0;
// 	/********** Filling x_vars **********/
// 	// for(int i = 0; i < n; i++) {
// 	// 	for(int j = 0; j < n; j++){
// 	// 		if (!nas->arcs[i][j]) {
// 	// 			continue;
// 	// 		}
// 	// 		for (int k1 = 0; k1 < nas->arcV[i][j].size(); k1++) {
// 	// 			int k = nas->arcV[i][j][k1];
// 	// 			x_vars.add(x[i][j][k]);
// 	// 		}
// 	// 	}
// 	// }
// 	/************************************/
// 	int V = nCustomers + 2*nParcels + v;

// 	for (int i = 0; i < V; i++){
// 		if (i >= nCustomers && i < nCustomers + 2*nParcels) {
// 			continue;
// 		}

//         for(int j = 0; j < n; ++j){
// 			if (j >= nCustomers && j < V) {
// 				continue;
// 			}

//             if (nas->arcs[i][j] != true){
//                 continue; // If arc i to j is invalid
//             }

// 			for (int a = 0; a < nas->subsequences[make_pair(i, j)].size(); a++) {
// 				for(int k1 = 0; k1 < nas->arcV[i][j].size(); k1++){
// 					int k = nas->arcV[i][j][k1];
// 					x_vars.add(x[i][j][a][k]);
// 					// // TODO UNCOMMENT //  << "x: [" << i << "][" << j << "][" << k << "]" << endl;
// 				}
// 			}
//         }
//     }
// }
// /*****************************************************************************************************************/

// /************************** Return a callback copy. This method is a CPLEX requirement ***************************/
// IloCplex::CallbackI* MyLazyCallback::duplicateCallback() const 
// { 
//    return new (getEnv()) MyLazyCallback(getEnv(), x, nas, mdist, inst, problem, nodeVec, bestSolVal, n, v, nParcels, nCustomers); 
// }
// /*****************************************************************************************************************/

// /************************************ Callback's code that is runned by CPLEX ************************************/
// void MyLazyCallback::main()
// {	
// 	// cout << "chegou no corte" << endl;
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
// 	map<pair<int, int>, int> arcChosen;

//     // Imprimir os valores das variáveis relaxadas
//     int index = 0;
//     for (int i = 0; i < n; i++) {
//         for (int j = 0; j < n; j++) {
// 			if (!nas->arcs[i][j]) {
// 				// // TODO UNCOMMENT // cout << i << " " << j << endl;
// 				continue;
// 			}

// 			pair<int, int> myPair = make_pair(i, j);

// 			for (int a = 0; a < nas->subsequences[myPair].size(); a++) {
// 				for (int k1 = 0; k1 < nas->arcV[i][j].size(); k1++)
// 				{
// 					int k = nas->arcV[i][j][k1];
// 					if (x_vals[index] > EPSILON)
// 					{
// 						arcDirection[i] = j;
// 						arcChosen[myPair] = a;
// 						// // TODO UNCOMMENT // cout << i << " " << j << " " << k << endl;
// 					}
// 					index++;
// 				}
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

// 		routes.emplace_back(); // Creating a new route

// 		routes[i].push_back(startDepotIndex);
// 		if (arcDirection.find(startDepotIndex) != arcDirection.end())
// 		{
// 			int prevNode = startDepotIndex;
// 			for (int nextNode = arcDirection[startDepotIndex]; nextNode != finalDepotIndex; nextNode = arcDirection[nextNode])
// 			{
// 				// // TODO UNCOMMENT // cout << "nextNode = " << nextNode << endl;
// 				pair<int, int> myPair = make_pair(prevNode, nextNode);

// 				for (int a = 0; a < nas->subsequences[myPair][arcChosen[myPair]].size(); a++) {
// 					routes[i].push_back(nas->subsequences[myPair][arcChosen[myPair]][a]);
// 				}
// 				routes[i].push_back(nextNode);
// 				prevNode = nextNode;
// 			}
// 		}
// 		routes[i].push_back(finalDepotIndex);
// 	}
// 	// // TODO UNCOMMENT // cout << "FIM ROTAS" << endl;

// 	// Printing the routes
// 	// cout << endl;
// 	// for (int i = 0; i < routes.size(); i++)
// 	// {
// 	// 	 cout << "Route " << i << ": ";
// 	// 	for (int j = 0; j < routes[i].size()-1; j++)
// 	// 	{
// 	// 		 cout << routes[i][j] << " -> ";
// 	// 	}
// 	// 	 cout << routes[i].back() << endl;
// 	// }
// 	// getchar();

// 	vector< vector<int> > customersRoutes;

// 	for (int k = 0; k < routes.size(); k++) {
// 		customersRoutes.push_back(vector<int>());
// 		for (int i = 0; i < routes[k].size(); i++) {
// 			if (routes[k][i] < inst->n || routes[k][i] >= inst->n + 2*inst->m) {
// 				customersRoutes[k].push_back(routes[k][i]);
// 			}
// 		}
// 	}

// 	int count = 0;

// 	double solVal = 0;
// 	double value = 0;

// 	for (int k = 0; k < customersRoutes.size(); k++) {
// 		//MIP
// 		//Creating environment and model 
// 		char var[100];
// 		IloEnv env;
// 		IloModel model(env, "nSARPLazy");
// 		int currSP;
// 		long M = 2*inst->T;
// 		long M2 = 2*(inst->n + inst->m + 1);
// 		long W = inst->m + 1;
// 		int Q;

// 		int fDepot = inst->n + 2*inst->m;
// 		int fDummy = inst->n + 2*inst->m + inst->K;
		
// 		int decimalPlaces = 4;
// 		double multiplier = std::pow(10, decimalPlaces);

// 		vector< pair<int, int> > auxPairVec;
// 		vector<int> servedParcels;
// 		pair<int, int> auxPair;

// 		// for (int i = 0; i < routes[k].size(); i++) {
// 		// 	int a = routes[k][i];

// 		// 	if (a >= inst->n && a < inst->n + 2*inst->m) {
// 		// 		servedParcels.push_back(a);
// 		// 		servedParcels.push_back(a + inst->m);
// 		// 	}
// 		// }
		
// 		// int P = 2*inst->m;

// 		// map<int, vector<int>> arcOut;
// 		// map<int, vector<int>> arcIn;
// 		// set<pair<int, int>> possArcs;

// 		// arcOut[customersRoutes[k].back()] = vector<int>();
// 		// arcIn[customersRoutes[k][0]] = vector<int>();

// 		// for (int i = 0; i < customersRoutes[k].size() - 1; i++) {
// 		// 	int cur = customersRoutes[k][i];
// 		// 	int next = customersRoutes[k][i+1];

// 		// 	arcOut[cur].push_back(next);
// 		// 	arcIn[next].push_back(cur);
// 		// 	pair<int, int> myPAir = make_pair(cur, next);
// 		// 	possArcs.insert(myPAir);

// 		// 	for (int j1 = 0; j1 < servedParcels.size(); j1++) {
// 		// 		int j = servedParcels[j1];

// 		// 		pair<int, int> myPAir = make_pair(cur, j);
// 		// 		arcOut[cur].push_back(j);
// 		// 		arcIn[j].push_back(cur);
// 		// 		possArcs.insert(myPAir);
// 		// 	}
// 		// }

// 		// for (int i = 0; i < servedParcels.size(); i++) {
// 		// 	int delivery = servedParcels[i];

// 		// 	for (int j1 = 1; j1 < customersRoutes[k].size(); j1++) {
// 		// 		int j = customersRoutes[k][j1];

// 		// 		pair<int, int> myPair = make_pair(delivery, j);
// 		// 		arcOut[delivery].push_back(j);
// 		// 		arcIn[j].push_back(delivery);
// 		// 		possArcs.insert(myPair);
// 		// 	}

// 		// 	for (int j1 = 0; j1 < servedParcels.size(); j1++) {
// 		// 		int j = servedParcels[j1];

// 		// 		if (j == delivery) {
// 		// 			continue;
// 		// 		}

// 		// 		pair<int, int> myPair = make_pair(delivery, j);
// 		// 		arcOut[delivery].push_back(j);
// 		// 		arcIn[j].push_back(delivery);
// 		// 		possArcs.insert(myPair);
// 		// 	}
// 		// }

// 		for (int i = 0; i < routes[k].size(); i++) {
// 			int a = routes[k][i];

// 			if (a >= inst->n && a < inst->n + 2*inst->m) {
// 				servedParcels.push_back(a);
// 			}
// 		}
		
// 		int P = 2*inst->m;

// 		map<int, vector<int>> arcOut;
// 		map<int, vector<int>> arcIn;
// 		set<pair<int, int>> possArcs;

// 		arcOut[customersRoutes[k].back()] = vector<int>();
// 		arcIn[customersRoutes[k][0]] = vector<int>();

// 		for (int i = 0; i < customersRoutes[k].size() - 1; i++) {
// 			int cur = customersRoutes[k][i];
// 			int next = customersRoutes[k][i+1];

// 			arcOut[cur].push_back(next);
// 			arcIn[next].push_back(cur);
// 			pair<int, int> myPAir = make_pair(cur, next);
// 			possArcs.insert(myPAir);

// 			for (int j1 = 0; j1 < servedParcels.size(); j1++) {
// 				int j = servedParcels[j1];

// 				pair<int, int> myPAir = make_pair(cur, j);
// 				arcOut[cur].push_back(j);
// 				arcIn[j].push_back(cur);
// 				possArcs.insert(myPAir);
// 			}
// 		}

// 		for (int i = 0; i < servedParcels.size(); i++) {
// 			int delivery = servedParcels[i];

// 			for (int j1 = 1; j1 < customersRoutes[k].size(); j1++) {
// 				int j = customersRoutes[k][j1];

// 				pair<int, int> myPair = make_pair(delivery, j);
// 				arcOut[delivery].push_back(j);
// 				arcIn[j].push_back(delivery);
// 				possArcs.insert(myPair);
// 			}

// 			for (int j1 = 0; j1 < servedParcels.size(); j1++) {
// 				int j = servedParcels[j1];

// 				if (j == delivery) {
// 					continue;
// 				}

// 				pair<int, int> myPair = make_pair(delivery, j);
// 				arcOut[delivery].push_back(j);
// 				arcIn[j].push_back(delivery);
// 				possArcs.insert(myPair);
// 			}
// 		}

// 		//Creating variables
// 		IloArray <IloBoolVarArray> u(env, nodeVec.size());

// 		for (int i = 0; i < nodeVec.size(); i++){
// 			if (arcOut.find(i) == arcOut.end()){
// 				continue; // If arc i to j is invalid
// 			}

// 			u[i] = IloBoolVarArray (env, nodeVec.size());
// 			for(int j = 0; j < nodeVec.size(); ++j){

// 				pair<int, int> myPair = make_pair(i, j);

// 				if (possArcs.find(myPair) == possArcs.end()){
// 					continue; // If arc i to j is invalid
// 				}

// 				sprintf(var, "u(%d,%d)", i, j);
// 				u[i][j].setName(var);
// 				model.add(u[i][j]);
// 				// // TODO UNCOMMENT //  << "x: [" << i << "][" << j << "][" << k << "]" << endl;
// 			}
// 		}

// 		// //creates boolean variable (y_i = 1 if node i is visited; 0 cc)
// 		// IloBoolVarArray y(env, adjList.size()); 

// 		// for (const auto& i : adjList){
// 		// 	sprintf(var, "y(%d)", i.first);
// 		// 	y[i.first].setName(var);
// 		// 	model.add(y[i.first]);
// 		// }

// 			// Variable start of service time
// 		IloNumVarArray b(env, nodeVec.size(), 9, inst->T);
// 		for (const auto i : arcOut){
// 			sprintf(var, "b(%d)", i.first);
// 			b[i.first].setName(var);
// 			model.add(b[i.first]);
// 		}

// 		IloExpr objFunction(env);

// 		objFunction += 0;

// 		// for (auto a : possArcs) {
// 		// 	int i = a.first;
// 		// 	int j = a.second;

// 		// 	// cout << i << " " << j << endl;

// 		// 	objFunction += nodeVec[i].profit * u[i][j];
// 		// 	objFunction -= (double)inst->costkm*mdist[i][j] * u[i][j];
// 		// }

// 		model.add(IloMaximize(env, objFunction));

// 		//Creating constraints

// 		//Constraint 1 - All elements must be used

// 		for (auto a : arcOut) {
// 			int i = a.first;

// 			if (i == customersRoutes[k].back()) {
// 				continue;
// 			}

// 			IloExpr exp(env);

// 			for (int j1 = 0; j1 < a.second.size(); j1++) {
// 				int j = a.second[j1];

// 				exp += u[i][j];
// 			}

// 			sprintf (var, "Constraint1_%d", i);
// 			IloRange cons = (exp == 1);
// 			cons.setName(var);
// 			model.add(cons);
// 		}

// 		{
// 			int i = customersRoutes[k].back();
			
// 			IloExpr exp(env);

// 			for (int j1 = 0; j1 < arcIn[i].size(); j1++) {
// 				int j = arcIn[i][j1];

// 				exp += u[j][i];
// 			}

// 			sprintf (var, "Constraint1_%d", i);
// 			IloRange cons = (exp == 1);
// 			cons.setName(var);
// 			model.add(cons);
// 		}

// 		//Constraint 4 - Flow conservation

// 		for (auto a : arcOut) {
// 			IloExpr exp1(env);
// 			IloExpr exp2(env);

// 			if (a.first >= inst->n + 2*inst->m) {
// 				continue;
// 			}

// 			int i = a.first;

// 			for (int b = 0; b < arcOut[i].size(); b++){
// 				int j = arcOut[i][b];
// 				exp1 += u[i][j];
// 			}
// 			//Right side: arc enters i
// 			for (int b = 0; b < arcIn[i].size(); b++){
// 				int j = arcIn[i][b];
// 				exp2 += u[j][i];
// 			}
// 			sprintf (var, "Constraint4_%d_%d", a, k);
// 			IloRange cons = ((exp1-exp2) == 0);
// 			cons.setName(var);
// 			model.add(cons);
// 		}

// 		// // Constraint 5 - route start in depot
// 		// {
// 		// 	IloExpr exp1(env);
// 		// 	IloExpr exp2(env);
// 		// 	int i = routes[k][0];

// 		// 	for (int b = 0; b < arcOut[i].size(); b++){
// 		// 		int j = arcOut[i][b];
// 		// 		exp1 += x[i][j];
// 		// 	}
// 		// 	//Right side: arc enters i
// 		// 	for (int b = 0; b < arcIn[i].size(); b++){
// 		// 		int j = arcOut[i][b];
// 		// 		exp2 += x[j][i];
// 		// 	}
// 		// 	sprintf (var, "Constraint4_%d_%d", a, k);
// 		// 	IloRange cons = ((exp1-exp2) == 0);
// 		// 	cons.setName(var);
// 		// 	model.add(cons);
// 		// }	

// 		//Constraint 8 - service of pickup must come before the delivery

// 		for (int i1 = 0; i1 < servedParcels.size(); i1++){
// 			IloExpr exp(env);
// 			int i = servedParcels[i1];

// 			if (i >= inst->n + inst->m) {
// 				continue;
// 			}

// 			exp = b[i] - b[i + inst->m];

// 			sprintf (var, "Constraint8_%d", i);
// 			IloRange cons = (exp <= 0);
// 			cons.setName(var);
// 			model.add(cons);
// 		}

// 		//Constraints 9 - TW 

// 		for (auto a : possArcs){
// 			IloExpr exp(env);
// 			IloExpr sumX(env);

// 			int i = a.first;
// 			int j = a.second;

// 			double cvalue = mdist[i][j]/inst->vmed;
// 			//cvalue = std::round(cvalue * multiplier) / multiplier;
// 			//cvalue = timeRound(cvalue);
// 			exp = b[i] - b[j] + nodeVec[i].delta + (cvalue) - M * (1 - u[i][j]);
// 			sprintf (var, "Constraint9_%d_%d", i, j);
// 			IloRange cons = (exp <= 0);
// 			cons.setName(var);
// 			model.add(cons);			
// 		}

// 		// Constraint 10 - ensure solution sequence
// 		for (int i = 0; i < routes[k].size() - 1; i++) {
// 			int u = routes[k][i];
// 			int v = routes[k][i + 1];
// 			IloExpr exp(env);

// 			exp += b[v] - b[u] - nodeVec[u].delta - mdist[u][v]/inst->vmed;

// 			sprintf (var, "Constraint10_%d", i);
// 			IloRange cons = (exp >= 0);
// 			cons.setName(var);
// 			model.add(cons);
// 		}

// 		//Constraints 11 and 12 - bound the service beginning time by the earlier and later service times for each node

// 		for (auto a : arcOut){
// 			int i = a.first;

// 			IloExpr exp(env);
// 			exp = b[i];

// 			sprintf (var, "Constraint11_%d", i);
// 			IloRange cons1 = (exp <= nodeVec[i].l);
// 			cons1.setName(var);
// 			model.add(cons1);
			
// 			sprintf (var, "Constraint12_%d", i);
// 			IloRange cons2 = (nodeVec[i].e <= exp);
// 			cons2.setName(var);
// 			model.add(cons2);			
// 		}

// 		//Constraints 13 - maximum driving time
// 		{
// 			int i = customersRoutes[k][0];
// 			IloExpr exp(env);
// 			exp = b[i + inst->K] - b[i];

// 			sprintf (var, "Constraint13_%d", i);
// 			IloRange cons1 = (exp <= inst->maxTime);
// 			cons1.setName(var);
// 			model.add(cons1);
// 		}    

// 		int threads;

// 		threads = 1;

// 		IloCplex nSARPLazy(model);
// 		nSARPLazy.exportModel("nSARPLazy.lp");
// 		nSARPLazy.setOut(env.getNullStream());
// 		// nSARPLazy.setWarning(env.getNullStream());
// 		nSARPLazy.setParam(IloCplex::Threads, threads);
// 		nSARPLazy.setParam(IloCplex::Param::TimeLimit, 7200);
// 		// nSARPLazy.setParam(IloCplex::Param::MIP::Display, 0);
// 		// nSARPLazy.setParam(IloCplex::Param::Read::WarningLimit, 0);

// 		IloNum start;
// 		IloNum time;
// 		start = nSARPLazy.getTime();
// 		nSARPLazy.solve();
// 		time = (nSARPLazy.getTime() - start)/threads;
// 		// cout  << "\nSol status: " << nSARPLazy.getStatus() << endl;
// 		// sStat->feasible = nSARPLazy.isPrimalFeasible();

// 		if (nSARPLazy.isPrimalFeasible()) {
// 			// cout << nSARPLazy.getObjValue() << endl;

// 			solVal += nSARPLazy.getObjValue();
			
// 			for (auto a : possArcs) {
// 				int i = a.first;
// 				int j = a.second;

// 				if (nSARPLazy.getValue(u[i][j]) > 0.5){
// 					auxPair.first = i;
// 					auxPair.second = j;
// 					// cout  << i << " " << j << ": " << nSARPLazy.getValue(u[i][j]) << endl;
// 					value += nodeVec[i].profit - mdist[i][j]*inst->costkm;
// 				}
// 			}	
// 			// getchar();

// 			count++;

// 			// cout << "FEASIBLE" << endl;
// 		} else {
// 			IloExpr expr(getEnv());	// Expression for the lazy constraint

// 			for (int i = 0; i < customersRoutes[k].size() - 1; i++) {
// 				int first = customersRoutes[k][i];
// 				int second = customersRoutes[k][i + 1];

// 				pair<int, int> myPair = make_pair(first, second);

// 				int a = arcChosen[myPair];

// 				expr += x[first][second][a][k];
// 			}

// 			add(expr <= (int)customersRoutes[k].size() - 2);

// 			// cout << "INFEASIBLE" << endl;
// 		}
// 		// getchar();

// 		// cout  << " Tree_Size: " <<  nSARP.getNnodes() + nSARP.getNnodesLeft() + 1 << endl;
// 		// cout  << " Total Time: " << time << endl;

// 		// if (sStat->feasible){

// 		// 	// cout  << " LB: " << nSARP.getObjValue() << endl;
// 		// 	// cout  << " UB: " << nSARP.getBestObjValue() << endl;
// 		// 	sStat->solprofit = nSARP.getObjValue();
// 		// 	sStat->solDual = nSARP.getBestObjValue();
// 		// 	sStat->time = time;

// 		// 	if (((nSARP.getBestObjValue() - nSARP.getObjValue())/nSARP.getBestObjValue()) * 100 < 0.01) {
// 		// 		sStat->status = "Optimal";
// 		// 	} else {
// 		// 		sStat->status = "Feasible";
// 		// 	}

// 		// 	for (int k = 0; k < inst->K; k++){
// 		// 		sStat->solvec.push_back(auxPairVec);
// 		// 	}

// 		// 	for (int i = 0; i < nodeVec.size(); i++){
// 		// 		for(int j = 0; j < nodeVec.size(); j++){                
// 		// 			if (nas->arcs[i][j] == true){
// 		// 				for (int k1 = 0; k1 < nas->arcV[i][j].size(); k1++){
// 		// 					int k = nas->arcV[i][j][k1];
				
// 		// 				}
// 		// 			}
// 		// 		}   
// 		// 	}

// 		// 	for (int i = 0; i < nodeVec.size(); i++){
// 		// 		if (nSARP.getValue(b[i]) > 0){
// 		// 			sStat->solBegin.push_back(nSARP.getValue(b[i]));
// 		// 		}
// 		// 		else {
// 		// 			sStat->solBegin.push_back(0);
// 		// 		}
// 		// 	}

// 		// 	for (int i = 0; i < nodeVec.size(); i++){
// 		// 		if (nSARP.getValue(w[i]) > 0.5){
// 		// 			sStat->solLoad.push_back(nSARP.getValue(w[i]));
// 		// 		}
// 		// 		else {
// 		// 			sStat->solLoad.push_back(0);
// 		// 		}
// 		// 	}


// 		// 	for (int i = 0; i < nodeVec.size(); i++){
// 		// 		//if (nSARP.getValue(u[i])){
// 		// 			sStat->solLoad2.push_back(nSARP.getValue(u[i]));
// 		// 		//}
// 		// 		//else {
// 		// 		//    sStat->solLoad2.push_back(0);
// 		// 		//}
// 		// 	}

// 		// }
// 		// if (problem->scen == "PC"){
// 		// 	nSARP.clearModel();
// 		// 	nSARP.end();
// 		// }

// 		env.end();
// 	}

// 	// cout << "saiu do corte" << endl;

// 	if (count == routes.size()) {
// 		// cout << "ALL ROUTES ARES FEASIBLE" << endl;
// 		// IloExpr expr(getEnv());	// Expression for the lazy constraint

// 		// cout << solVal << endl;
// 		*bestSolVal = max(*bestSolVal, solVal);

// 		// cout << "Cur best = " << *bestSolVal << endl;
// 		// getchar();

// 		// int _size = 0;
		
// 		// for (int k = 0; k < customersRoutes.size(); k++) {
// 		// 	_size += customersRoutes[k].size() - 1;
// 		// 	for (int i = 0; i < customersRoutes[k].size() - 1; i++) {
// 		// 		int first = customersRoutes[k][i];
// 		// 		int second = customersRoutes[k][i + 1];

// 		// 		pair<int, int> myPair = make_pair(first, second);

// 		// 		int a = arcChosen[myPair];

// 		// 		expr += x[first][second][a][k];
// 		// 	}
// 		// }

// 		// // cout << value << endl;
// 		// // getchar();
// 		// add(expr <= _size - 1);
// 	}


// 	/**************** Creating the Lazy Co	// // TODO UNCOMMENT // cout << "aqui 2" << endl;
// 			for (int j = 0; j < routes[i].size() - 1; j++) {
// 				expr += x[ routes[i][j] ][ routes[i][j+1] ][ i ];
// 			}

// 			// // TODO UNCOMMENT // cout << "aqui 3" << endl;

// 			int maxArcs = routes[i].size() - 2;
// 			add(expr <= maxArcs);
// 		}nstraints ****************/
// 	// // TODO UNCOMMENT // cout << "SIZES: " << x.getSize() << " " << x[0].getSize() << endl;

// 	// for (int i = 0; i < routes.size(); i++)
// 	// {
// 	// 	int counter = 0;

// 	// 	for (int j = 0; j < routes[i].size(); j++)
// 	// 	{
// 	// 		if (isParcel(routes[i][j]))
// 	// 		{
// 	// 			counter++;

// 	// 			if (counter >= 4 && isDeliveryParcel(routes[i][j-3]) && isPickUpParcel(routes[i][j]))
// 	// 			{
// 	// 				// Prohibiting the sequence for all vehicles
// 	// 				// // TODO UNCOMMENT // cout << endl;
// 	// 				// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// 	// 				IloExpr expr(getEnv());
// 	// 				for (int vehicle = 0; vehicle < this->v; vehicle++)
// 	// 				{
// 	// 					// int k = arrConj->arcV[i][j][vehicle];

// 	// 					// TODO UNCOMMENT // cout << endl;
// 	// 					// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

// 	// 					// // // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// 	// 					// expr += x[ routes[i][j-5] ][ routes[i][j-4] ][ vehicle ];

// 	// 					// // // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// 	// 					// expr += x[ routes[i][j-4] ][ routes[i][j-3] ][ vehicle ];

// 	// 					// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// 	// 					expr += x[ routes[i][j-3] ][ routes[i][j-2] ][ vehicle ];

// 	// 					// // TODO UNCOMMENT // cout << "x(" << routes[i][j-2] << "," << routes[i][j-1] << "," << vehicle << ")" << endl;
// 	// 					expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// 	// 					// // TODO UNCOMMENT // cout << "x(" << routes[i][j-1] << "," << routes[i][j] << "," << vehicle << ")" << endl;
// 	// 					expr += x[ routes[i][j-1] ][ routes[i][ j ] ][ vehicle ];

// 	// 					// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// 	// 					// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// 	// 					// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// 	// 					// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// 	// 					// }
// 	// 					// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// 	// 					// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// 	// 					// }
// 	// 					// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// 	// 					// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// 	// 				}
// 	// 				add(expr <= 2);
// 	// 				// // TODO UNCOMMENT // cout << endl;
// 	// 				// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// 	// 			}
// 	// 		}
// 	// 		else // isCustomer = true
// 	// 		{
// 	// 			counter = 0;
// 	// 		}
// 	// 	}
// 	// }

// 	// Prohibiting sequence V-P-D-P-D
// 	// for (int i = 0; i < routes.size(); i++)
// 	// {
// 	// 	int counter = 0;

// 	// 	for (int j = 0; j < routes[i].size() && j < 4; j++)
// 	// 	{
// 	// 		if (isParcel(routes[i][j]))
// 	// 		{
// 	// 			counter++;

// 	// 			if (counter >= 3)
// 	// 			{
// 	// 				// Prohibiting the sequence for all vehicles
// 	// 				// // TODO UNCOMMENT // cout << endl;
// 	// 				// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// 	// 				for (int vehicle = 0; vehicle < this->v; vehicle++)
// 	// 				{
// 	// 					// int k = arrConj->arcV[i][j][vehicle];
// 	// 					// // TODO UNCOMMENT // cout << endl;
// 	// 					// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;
// 	// 					if (arrConj->arcs[ vehicle ][ routes[i][j-2] ]) {
// 	// 						IloExpr expr(getEnv()); //!! Isso estava fora do "if". Eu movi para dentro.
// 	// 												//!! Mas por que não colocar fora do laço de repetição?

// 	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// 	// 						expr += x[ routes[vehicle][0] ][ routes[i][j-2] ][ vehicle ];

// 	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// 	// 						expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// 	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// 	// 						expr += x[ routes[i][j-1] ][ routes[i][j] ][ vehicle ];

// 	// 						// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// 	// 						// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// 	// 						// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// 	// 						// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// 	// 						// }
// 	// 						// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// 	// 						// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// 	// 						// }
// 	// 						// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// 	// 						// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// 	// 						add(expr <= 2); //!! Isso estava fora do "if". Eu movi para dentro
// 	// 					}
// 	// 				}
// 	// 				// // TODO UNCOMMENT // cout << endl;
// 	// 				// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// 	// 			}
// 	// 		}
// 	// 		else // isCustomer = true
// 	// 		{
// 	// 			counter = 0;
// 	// 		}
// 	// 	}
// 	// }

// 	// for (int i = 0; i < routes.size(); i++)
// 	// {
// 	// 	int counter = 0;

// 	// 	for (int j = routes[i].size() - 4; j < routes[i].size(); j++)
// 	// 	{
// 	// 		if (isParcel(routes[i][j]))
// 	// 		{
// 	// 			counter++;

// 	// 			if (counter >= 3)
// 	// 			{
// 	// 				// Prohibiting the sequence for all vehicles
// 	// 				// // TODO UNCOMMENT // cout << endl;
// 	// 				// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// 	// 				for (int vehicle = 0; vehicle < this->v; vehicle++)
// 	// 				{
// 	// 					// int k = arrConj->arcV[i][j][vehicle];
// 	// 					// // TODO UNCOMMENT // cout << endl;
// 	// 					// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

// 	// 					if (arrConj->arcs[ routes[i][j] ][ routes[vehicle].back() ]) {
// 	// 						IloExpr expr(getEnv()); //!! Isso estava fora do "if". Eu movi para dentro.
// 	// 												//!! Mas por que não colocar fora do laço de repetição?

// 	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// 	// 						expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// 	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// 	// 						expr += x[ routes[i][j-1] ][ routes[i][j] ][ vehicle ];

// 	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// 	// 						expr += x[ routes[i][j] ][ routes[vehicle].back() ][ vehicle ];

// 	// 						// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// 	// 						// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// 	// 						// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// 	// 						// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// 	// 						// }
// 	// 						// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// 	// 						// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// 	// 						// }
// 	// 						// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// 	// 						// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// 	// 						add(expr <= 2); //!! Isso estava fora do "if". Eu movi para dentro
// 	// 					}
// 	// 				}
// 	// 				// // TODO UNCOMMENT // cout << endl;
// 	// 				// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// 	// 			}
// 	// 		}
// 	// 		else // isCustomer = true
// 	// 		{
// 	// 			counter = 0;
// 	// 		}
// 	// 	}
// 	// }

// 	// bool restricted = false;

// 	// for (int i = 0; i < routes.size(); i++)
// 	// {
// 	// 	// // TODO UNCOMMENT // cout << "aqui" << endl;
// 	// 	for (int j = 1; j < routes[i].size(); j++) {
// 	// 		if (j + 1 < routes[i].size() && isPickUpParcel(j) && isDeliveryParcel(j + 1)) {
// 	// 			int k;
// 	// 			for (k = j+1; k + 3 < routes[i].size();) {
// 	// 				if (isCustomer(k+1) && isPickUpParcel(k+2) && isDeliveryParcel(k+3)) {
// 	// 					k += 3;
// 	// 				} else {
// 	// 					break;
// 	// 				}
// 	// 			}

// 	// 			if (!(isCustomer(j-1) || isCustomer(k + 1))) {
// 	// 				IloExpr expr(getEnv());

// 	// 				for (int vehicle = 0; vehicle < this->v; vehicle++) {
// 	// 					for (int k1 = j; k1 < k; k1++) {
// 	// 						if (arrConj->arcs[k1][k1+1]) expr += x[k1][k1+1][vehicle];
// 	// 					}

// 	// 					expr += x[j-1][j][vehicle];
// 	// 					expr += x[k][k+1][vehicle];

// 	// 					for (int k1 = 0; k1 < this->n; k1++) {
// 	// 						if (isDeliveryParcel(k1) && arrConj->arcs[k1][j]) {
// 	// 							expr += x[k1][j][vehicle];
// 	// 						}

// 	// 						if (isPickUpParcel(k1) && arrConj->arcs[k][k1]) {
// 	// 							expr += x[k][k1][vehicle];
// 	// 						}
// 	// 					}
// 	// 				}

// 	// 				add(expr <= k - j);
// 	// 				restricted = true;
// 	// 			}
				
// 	// 			// j = k+1;
// 	// 		}
// 	// 	}
// 	// }


// 	// Lazy Constraint
// 	// In the case P - D - c - P - D, at least one of the parcels must be connected to
// 	// another customer

// 	// for (int i = 0; i < routes.size(); i++)
// 	// {
// 	// 	int counter = 0;
	
// 	// 	for (int j = 0; j < routes[i].size(); j++)
// 	// 	{
// 	// 		if (isPickUpParcel(routes[i][j]))
// 	// 		{
// 	// 			if ((counter == 0 || counter == 2) && isDeliveryParcel(routes[i][j+1]))	counter++;
// 	// 			else																	counter = 0;

// 	// 			// If counter == 3, then we found a pattern P - D - c - P - D
// 	// 			if (counter == 3)
// 	// 			{
// 	// 				if (!isCustomer(routes[i][j-4]) && !isCustomer(routes[i][j+2]))
// 	// 				{
// 	// 					// TODO UNCOMMENT // cout 	<< "EXPRESSAO P - D - c - P - D: "
// 	// 							<< routes[i][j-3] << " " << routes[i][j-2] << " " << routes[i][j-1] << " "
// 	// 							<< routes[i][j] << " " << routes[i][j+1] << endl;

// 	// 					IloExpr expr(getEnv());

// 	// 					for (int k = 0; k < this->v; k++)
// 	// 					{
// 	// 						expr += x[ routes[i][j-3] ][ routes[i][j-2] ][ k ];
// 	// 						expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ k ];
// 	// 						expr += x[ routes[i][j-1] ][ routes[i][j] ][ k ];
// 	// 						expr += x[ routes[i][j] ][ routes[i][j+1] ][ k ];
// 	// 					}

// 	// 					for (int customer = 0; customer < nCustomers; customer++)
// 	// 					{
// 	// 						for (int k = 0; k < this->v; k++)
// 	// 						{
// 	// 							if (arrConj->arcs[customer][routes[i][j-3]]) expr -= x[ customer ][ routes[i][j-3] ][ k ];
// 	// 							if (arrConj->arcs[routes[i][j+1]][customer]) expr -= x[ routes[i][j+1] ][ customer ][ k ];
// 	// 						}
// 	// 					}

// 	// 					////////////////////////////////
// 	// 					add(expr <= 3);
// 	// 				}
					
// 	// 				////////////////////////
// 	// 				counter = 1;
// 	// 			}
// 	// 		}
// 	// 		else if (isCustomer(routes[i][j]))
// 	// 		{
// 	// 			if (counter == 1) 	counter++;
// 	// 			else				counter = 0;
// 	// 		}
// 	// 	}
// 	// }


// 	// for (int k = 0; k < routes.size(); k++) {
// 	// 	for (int i = 0; i < routes[k].size(); i++) {
// 	// 		cout << routes[k][i] << " ";
// 	// 	}
// 	// 	cout << endl;
// 	// }
// 	// getchar();

// 	// Creating restrictions for every subsequence that violates the formulation
// 	// NO PERMUTATION!

// 	// for (int r = 0; r < routes.size(); r++)
// 	// {
// 	// 	for (int size = 2; size <= (routes[r].size()-2); size++)
// 	// 	{
// 	// 		for (int startIndex = 1; startIndex < (routes[r].size()-1); startIndex++)
// 	// 		{
// 	// 			int lastIndex = startIndex + size - 1;
// 	// 			if (lastIndex >= routes[r].size()-1)
// 	// 			{
// 	// 				break;
// 	// 			}
// 	// 			// getchar();

// 	// 			// If it's a sequence of type P - D - ... - P - D, or simply P - D
// 	// 			if (isPickUpParcel(routes[r][startIndex]) && isDeliveryParcel(routes[r][startIndex+1]) && isPickUpParcel(routes[r][lastIndex-1]) && isDeliveryParcel(routes[r][lastIndex]))
// 	// 			{
// 	// 				int counter	= 0;

// 	// 				for (int i = startIndex; i <= lastIndex; i++)
// 	// 				{
// 	// 					if (isPickUpParcel(routes[r][i]))	counter--;
// 	// 					else if (isCustomer(routes[r][i]))	counter++;
// 	// 				}

// 	// 				if ((counter + isCustomer(routes[r][startIndex-1]) + isCustomer(routes[r][lastIndex+1])) < 0)
// 	// 				{
// 	// 					// Creating the restriction
// 	// 					IloExpr expr(getEnv());

// 	// 					for (int currNode = startIndex; currNode <= lastIndex; currNode++)
// 	// 					{
// 	// 						for (int nextNode = startIndex; nextNode <= lastIndex; nextNode++) {
// 	// 							if (!arrConj->arcs[ routes[r][currNode] ][ routes[r][nextNode] ]) continue;

// 	// 							for (int k1 = 0; k1 < arrConj->arcV[ routes[r][currNode] ][ routes[r][nextNode] ].size(); k1++)
// 	// 							{
// 	// 								int k = arrConj->arcV[ routes[r][currNode] ][ routes[r][nextNode] ][ k1 ];
// 	// 								expr += x[ routes[r][currNode] ][ routes[r][nextNode] ][ k ];
// 	// 							}
// 	// 						}
// 	// 						// cout << currNode << " ";
// 	// 					}

						
// 	// 					// int start;

// 	// 					// if (routes[r][startIndex - 1] < nCustomers) start = 0;
// 	// 					// else start = nCustomers;

// 	// 					// for (int left = start; left < n; left++) {
// 	// 					// 	if (!arrConj->arcs[ left ][ routes[r][startIndex] ]) continue;

// 	// 					// 	for (int k1 = 0; k1 < arrConj->arcV[ left ][ routes[r][startIndex] ].size(); k1++)
// 	// 					// 	{
// 	// 					// 		int k = arrConj->arcV[ left ][ routes[r][startIndex] ][ k1 ];
// 	// 					// 		expr += x[ left ][ routes[r][startIndex] ][ k ];
// 	// 					// 	}
// 	// 					// }

// 	// 					// if (routes[r][lastIndex + 1] < nCustomers) start = 0;
// 	// 					// else start = nCustomers;
						
// 	// 					// for (int right = start; right < n; right++) {
// 	// 					// 	if (!arrConj->arcs[ routes[r][lastIndex] ][ right ]) continue;

// 	// 					// 	for (int k1 = 0; k1 < arrConj->arcV[ routes[r][lastIndex] ][ right ].size(); k1++)
// 	// 					// 	{
// 	// 					// 		int k = arrConj->arcV[ routes[r][lastIndex] ][ right ][ k1 ];
// 	// 					// 		expr += x[ routes[r][lastIndex] ][ right ][ k ];
// 	// 					// 	}
// 	// 					// }

// 	// 					if (arrConj->arcs[ routes[r][startIndex - 1] ][ routes[r][startIndex] ]) {
// 	// 						for (int k1 = 0; k1 < arrConj->arcV[ routes[r][startIndex-1] ][ routes[r][startIndex] ].size(); k1++)
// 	// 						{
// 	// 							int k = arrConj->arcV[ routes[r][startIndex - 1] ][ routes[r][startIndex] ][ k1 ];
// 	// 							expr += x[ routes[r][startIndex - 1] ][ routes[r][startIndex] ][ k ];
// 	// 						}
// 	// 					}

// 	// 					if (arrConj->arcs[ routes[r][lastIndex] ][ routes[r][lastIndex + 1] ]) {
// 	// 						for (int k1 = 0; k1 < arrConj->arcV[ routes[r][lastIndex] ][ routes[r][lastIndex + 1] ].size(); k1++)
// 	// 						{
// 	// 							int k = arrConj->arcV[ routes[r][lastIndex] ][ routes[r][lastIndex + 1] ][ k1 ];
// 	// 							expr += x[ routes[r][lastIndex] ][ routes[r][lastIndex + 1] ][ k ];
// 	// 						}
// 	// 					}
						
// 	// 					// cout << lastIndex << endl;
// 	// 					// getchar();

// 	// 					// for (int c = 0; c < nCustomers; c++)
// 	// 					// {
// 	// 					// 	if (arrConj->arcs[ c ][ routes[r][startIndex] ])
// 	// 					// 	{
// 	// 					// 		for (int k1 = 0; k1 < arrConj->arcV[ c ][ routes[r][startIndex] ].size(); k1++)
// 	// 					// 		{
// 	// 					// 			int k = arrConj->arcV[ c ][ routes[r][startIndex] ][ k1 ];
// 	// 					// 			expr -= x[ c ][ routes[r][startIndex] ][ k ];
// 	// 					// 		}
// 	// 					// 	}

// 	// 					// 	if (arrConj->arcs[ routes[r][lastIndex] ][ c ])
// 	// 					// 	{
// 	// 					// 		for (int k1 = 0; k1 < arrConj->arcV[ routes[r][lastIndex] ][ c ].size(); k1++)
// 	// 					// 		{
// 	// 					// 			int k = arrConj->arcV[ routes[r][lastIndex] ][ c ][ k1 ];
// 	// 					// 			expr -= x[ routes[r][lastIndex] ][ c ][ k ];
// 	// 					// 		}
// 	// 					// 	}
// 	// 					// }

// 	// 					expr -= size;
// 	// 					// expr -= counter;

// 	// 					// Adding the constraint
// 	// 					add(expr <= 0);

// 	// 					// Printing the sebsequence
// 	// 					cout << endl;
// 	// 					cout << "route: ";
// 	// 					for (int i = 0; i < routes[r].size(); i++)
// 	// 					{
// 	// 						cout << routes[r][i] << " ";
// 	// 					}
// 	// 					cout << endl;

// 	// 					cout << endl;
// 	// 					cout << "sequence: ";
// 	// 					for (int i = startIndex; i <= lastIndex; i++)
// 	// 					{
// 	// 						cout << routes[r][i] << " ";
// 	// 					}
// 	// 					cout << endl;
// 	// 					cout << "sequence: ";
// 	// 					for (int i = startIndex; i <= lastIndex; i++)
// 	// 					{
// 	// 						if (isCustomer(routes[r][i])) 		cout << "c";
// 	// 						if (isPickUpParcel(routes[r][i])) 	cout << "P";
// 	// 						if (isDeliveryParcel(routes[r][i])) cout << "D";

// 	// 						if (i != lastIndex)	cout << " - ";
// 	// 					}
// 	// 					cout << endl;
// 	// 					cout << endl;
// 	// 					// getchar();
// 	// 				}
// 	// 			}
// 	// 		}
// 	// 	}
// 	// }






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


// // bool MyLazyCallback::isCustomer(int index) {
// // 	return (index >= 0) && (index < this->nCustomers);
// // }


// // bool MyLazyCallback::isDepot(int index) {
// // 	return (index >= this->nCustomers + 2*this->nParcels) && (index < this->nCustomers + 2*this->nParcels + 2*this->v);
// // }

// // bool MyLazyCallback::isParcel(int index) {
// // 	return (index >= this->nCustomers) && (index < (this->nCustomers + 2*this->nParcels));
// // }

// // bool MyLazyCallback::isPickUpParcel(int index) {
// // 	return (index >= this->nCustomers) && (index < (this->nCustomers + this->nParcels));
// // }

// // bool MyLazyCallback::isDeliveryParcel(int index) {
// // 	return (index >= this->nCustomers + this->nParcels) && (index < (this->nCustomers + 2*this->nParcels));
// // }

// // bool MyLazyCallback::finalVerifier(vector<int> &route) {
// // 	int contP = 0;
// // 	int contD = 0;

// // 	for (int i = 0; i < route.size(); i++) {
// // 		if (isPickUpParcel(route[i])) {
// // 			contP--;
// // 		} else if (isDeliveryParcel(route[i])) {
// // 			contD--;
// // 		} else if (isCustomer(route[i])){
// // 			contP = min(1, contP + 1);
// // 			contD = min(1, contD + 1);
// // 		}

// // 		if (contP < -1 || contD < -1) {
// // 			return false;
// // 		}
// // 	}

// // 	if (contP == -1 || contD == -1) {
// // 		return false;
// // 	}

// // 	return true;
// // }

// // #include "MyLazyCallback.h"
// // #include <map>

// // /********************************************** Class' Constructor **********************************************/
// // MyLazyCallback::MyLazyCallback(IloEnv env, const IloArray<IloArray<IloBoolVarArray>>& x_ref, nodeArcsStruct *nas, int nodes, int veic, int parcels, int customers) : IloCplex::LazyConstraintCallbackI(env), x(x_ref), x_vars(env), n(nodes), v(veic), arrConj(nas), nParcels(parcels), nCustomers(customers)
// // {
// // 	int num = 0;
// // 	/********** Filling x_vars **********/
// // 	for(int i = 0; i < n; i++) {
// // 		for(int j = 0; j < n; j++){
// // 			if (!nas->arcs[i][j]) {
// // 				continue;
// // 			}
// // 			for (int k1 = 0; k1 < nas->arcV[i][j].size(); k1++) {
// // 				int k = nas->arcV[i][j][k1];
// // 				x_vars.add(x[i][j][k]);
// // 			}
// // 		}
// // 	}
// // 	/************************************/
// // }
// // /*****************************************************************************************************************/

// // /************************** Return a callback copy. This method is a CPLEX requirement ***************************/
// // IloCplex::CallbackI* MyLazyCallback::duplicateCallback() const 
// // { 
// //    return new (getEnv()) MyLazyCallback(getEnv(), x, arrConj, n, v, nParcels, nCustomers); 
// // }
// // /*****************************************************************************************************************/

// // /************************************ Callback's code that is runned by CPLEX ************************************/
// // void MyLazyCallback::main()
// // {	
// // 	// /********** Getting the relaxed variables values **********/
// // 	// IloNumArray x_vals(getEnv(), (0.5*(n)*(n-1)));
// // 	// getValues(x_vals, x_vars);
// // 	// /**********************************************************/
   
// // 	// vector <vector<int> > cutSetPool;
// // 	// vector<IloConstraint> cons; 

// // 	IloNumArray x_vals(getEnv());
// //     getValues(x_vals, x_vars);

// // 	// Matrix to store every route in the current solution
// // 	vector<vector<int>> routes;
// // 	map<int, int> arcDirection;

// // 	bool restricted = false;

// //     // Imprimir os valores das variáveis relaxadas
// //     int index = 0;
// //     for (int i = 0; i < n; i++) {
// //         for (int j = 0; j < n; j++) {
// // 			if (!arrConj->arcs[i][j]) {
// // 				// // TODO UNCOMMENT // cout << i << " " << j << endl;
// // 				continue;
// // 			}

// // 			for (int k1 = 0; k1 < arrConj->arcV[i][j].size(); k1++)
// // 			{
// // 				int k = arrConj->arcV[i][j][k1];
// // 				if (x_vals[index] > EPSILON)
// // 				{
// // 					arcDirection[i] = j;
// // 					// // TODO UNCOMMENT // cout << i << " " << j << " " << k << endl;
// // 				}
// // 				index++;
// // 			}
// //         }
// //     }

// //     // Liberar a memória
// //     x_vals.end();

// // 	// Remember: "this->v" is the numver of available vehicles
// // 	// // TODO UNCOMMENT // cout << endl;
// // 	// // TODO UNCOMMENT // cout << "ROTAS" << endl;
// // 	// // TODO UNCOMMENT // cout << "numCustomers = " << this->nCustomers << endl;
// // 	// // TODO UNCOMMENT // cout << "numParcels = " << this->nParcels << endl;
// // 	// // TODO UNCOMMENT // cout << "numVehicles = " << this->v << endl;
// // 	for (int i = 0; i < this->v; i++)
// // 	{
// // 		int startDepotIndex = this->nCustomers + 2*this->nParcels + i;	// Start depot for the current vehicle
// // 		int finalDepotIndex = startDepotIndex + this->v;				// Final depot for the current vehicle

// // 		// // TODO UNCOMMENT // cout << endl;
// // 		// // TODO UNCOMMENT // cout << "i = " << i << endl;
// // 		// // TODO UNCOMMENT // cout << "startDepotIndex = " << startDepotIndex << endl;
// // 		// // TODO UNCOMMENT // cout << "finalDepotIndex = " << finalDepotIndex << endl;
// // 		// // TODO UNCOMMENT // cout << "arcDirection[" << startDepotIndex << "] = " << arcDirection[startDepotIndex] << endl;
// // 		// // TODO UNCOMMENT // cout << endl;

// // 		if (arcDirection.find(startDepotIndex) != arcDirection.end())
// // 		{
// // 			routes.emplace_back(); // Creating a new route

// // 			routes[i].push_back(startDepotIndex);
// // 			for (int nextNode = arcDirection[startDepotIndex]; nextNode != finalDepotIndex; nextNode = arcDirection[nextNode])
// // 			{
// // 				// // TODO UNCOMMENT // cout << "nextNode = " << nextNode << endl;
// // 				routes[i].push_back(nextNode);
// // 			}
// // 			routes[i].push_back(finalDepotIndex);
// // 		}
// // 	}
// // 	// // TODO UNCOMMENT // cout << "FIM ROTAS" << endl;

// // 	// Printing the routes
// // 	// // TODO UNCOMMENT // cout << endl;
// // 	// for (int i = 0; i < routes.size(); i++)
// // 	// {
// // 	// 	// TODO UNCOMMENT // cout << "Route " << i << ": ";
// // 	// 	for (int j = 0; j < routes[i].size()-1; j++)
// // 	// 	{
// // 	// 		// TODO UNCOMMENT // cout << routes[i][j] << " -> ";
// // 	// 	}
// // 	// 	// TODO UNCOMMENT // cout << routes[i].back() << endl;
// // 	// }
// // 	// getchar();


// // 	/**************** Creating the Lazy Constraints ****************/
// // 	// // TODO UNCOMMENT // cout << "SIZES: " << x.getSize() << " " << x[0].getSize() << endl;

// // 	for (int i = 0; i < routes.size(); i++)
// // 	{
// // 		int counter = 0;

// // 		for (int j = 0; j < routes[i].size(); j++)
// // 		{
// // 			if (isParcel(routes[i][j]))
// // 			{
// // 				counter++;

// // 				if (counter >= 4 && isDeliveryParcel(routes[i][j-3]) && isPickUpParcel(routes[i][j]))
// // 				{
// // 					vector<pair<int, int>> pairs = makePairs(routes[i], j - 3);
// // 					// Prohibiting the sequence for all vehicles
// // 					// // TODO UNCOMMENT // cout << endl;
// // 					// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// // 					IloExpr expr(getEnv());
// // 					for (int vehicle = 0; vehicle < this->v; vehicle++)
// // 					{
// // 						// int k = arrConj->arcV[i][j][vehicle];

// // 						// TODO UNCOMMENT // cout << endl;
// // 						// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

// // 						// // // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// // 						// expr += x[ routes[i][j-5] ][ routes[i][j-4] ][ vehicle ];

// // 						// // // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// // 						// expr += x[ routes[i][j-4] ][ routes[i][j-3] ][ vehicle ];

// // 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// // 						// expr += x[ routes[i][j-3] ][ routes[i][j-2] ][ vehicle ];

// // 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-2] << "," << routes[i][j-1] << "," << vehicle << ")" << endl;
// // 						expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// // 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-1] << "," << routes[i][j] << "," << vehicle << ")" << endl;
// // 						// expr += x[ routes[i][j-1] ][ routes[i][ j ] ][ vehicle ];

// // 						for (int i1 = 0; i1 < pairs.size(); i1++) {
// // 							for (int j1 = 0; j1 < pairs.size(); j1++) {
// // 								if (i1 != j1 && arrConj->arcs[ pairs[i1].second ][ pairs[j1].first ])  {
// // 									expr += x[ pairs[i1].second ][ pairs[j1].first ][ vehicle ];
// // 								}
// // 							}
// // 						}

// // 						// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// // 						// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// // 						// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// // 						// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// // 						// }
// // 						// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// // 						// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// // 						// }
// // 						// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// // 						// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// // 					}
// // 					add(expr <= 2);
// // 					// // TODO UNCOMMENT // cout << endl;
// // 					// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// // 				}
// // 			}
// // 			else // isCustomer = true
// // 			{
// // 				counter = 0;
// // 			}
// // 		}
// // 	}

// // 	// Prohibiting sequence V-P-D-P-D
// // 	for (int i = 0; i < routes.size(); i++)
// // 	{
// // 		int counter = 0;

// // 		for (int j = 0; j < routes[i].size() && j < 4; j++)
// // 		{
// // 			if (isParcel(routes[i][j]))
// // 			{
// // 				counter++;

// // 				if (counter >= 3)
// // 				{
// // 					// Prohibiting the sequence for all vehicles
// // 					// // TODO UNCOMMENT // cout << endl;
// // 					// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// // 					for (int vehicle = 0; vehicle < this->v; vehicle++)
// // 					{
// // 						// int k = arrConj->arcV[i][j][vehicle];
// // 						IloExpr expr(getEnv());
// // 						// TODO UNCOMMENT // cout << endl;
// // 						// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

// // 						if (arrConj->arcs[ routes[vehicle][0] ][ routes[i][j-2] ]) {
// // 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// // 							expr += x[ routes[vehicle][0] ][ routes[i][j-2] ][ vehicle ];

// // 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// // 							expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// // 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// // 							expr += x[ routes[i][j-1] ][ routes[i][j] ][ vehicle ];

// // 							if (arrConj->arcs[ routes[i][j]+nParcels ][ routes[i][j-2] ] ) {
// // 								if (arrConj->arcs[ routes[vehicle][0] ][ routes[i][j] ]) {
// // 									expr += x[ routes[i][j]+nParcels ][ routes[i][j-2] ][ vehicle ];

// // 									expr += x[ routes[vehicle][0] ][ routes[i][j] ][ vehicle ];
// // 								}
// // 							}

// // 							// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// // 							// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// // 							// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// // 							// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// // 							// }
// // 							// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// // 							// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// // 							// }
// // 							// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// // 							// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// // 						}
// // 						add(expr <= 2);
// // 					}
// // 					// // TODO UNCOMMENT // cout << endl;
// // 					// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// // 				}
// // 			}
// // 			else // isCustomer = true
// // 			{
// // 				counter = 0;
// // 			}
// // 		}
// // 	}

// // 	for (int i = 0; i < routes.size(); i++)
// // 	{
// // 		int counter = 0;

// // 		for (int j = routes[i].size() - 4; j < routes[i].size(); j++)
// // 		{
// // 			if (isParcel(routes[i][j]))
// // 			{
// // 				counter++;

// // 				if (counter >= 3)
// // 				{
// // 					// Prohibiting the sequence for all vehicles
// // 					// // TODO UNCOMMENT // cout << endl;
// // 					// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
// // 					for (int vehicle = 0; vehicle < this->v; vehicle++)
// // 					{
// // 						// int k = arrConj->arcV[i][j][vehicle];
// // 						IloExpr expr(getEnv());
// // 						// TODO UNCOMMENT // cout << endl;
// // 						// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

// // 						if (arrConj->arcs[ routes[i][j] ][ routes[vehicle].back() ]) {
// // 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
// // 							expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

// // 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
// // 							expr += x[ routes[i][j-1] ][ routes[i][j] ][ vehicle ];

// // 							// // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
// // 							expr += x[ routes[i][j] ][ routes[vehicle].back() ][ vehicle ];

// // 							if (arrConj->arcs[ routes[i][j] ][ routes[i][j-2] - nParcels ] ) {
// // 								if (arrConj->arcs[ routes[i][j-2] ][ routes[vehicle].back() ]) {
// // 									expr += x[ routes[i][j] ][ routes[i][j-2] - nParcels ][ vehicle ];

// // 									expr += x[ routes[i][j-2] ][ routes[vehicle].back() ][ vehicle ];
// // 								}
// // 							}

// // 							// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

// // 							// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
// // 							// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
// // 							// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
// // 							// }
// // 							// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
// // 							// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
// // 							// }
// // 							// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
// // 							// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
// // 						}
// // 						add(expr <= 2);
// // 					}
// // 					// // TODO UNCOMMENT // cout << endl;
// // 					// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
// // 				}
// // 			}
// // 			else // isCustomer = true
// // 			{
// // 				counter = 0;
// // 			}
// // 		}
// // 	}

// // 	for (int i = 0; i < routes.size(); i++)
// // 	{
// // 		// // TODO UNCOMMENT // cout << "aqui" << endl;
// // 		if (!finalVerifier(routes[i])) {
// // 			IloExpr expr(getEnv());
			
// // 			// // TODO UNCOMMENT // cout << "aqui 2" << endl;
// // 			for (int j = 0; j < routes[i].size() - 1; j++) {
// // 				expr += x[ routes[i][j] ][ routes[i][j+1] ][ i ];
// // 			}

// // 			// // TODO UNCOMMENT // cout << "aqui 3" << endl;

// // 			int maxArcs = routes[i].size() - 2;
// // 			add(expr <= maxArcs);
// // 		}
// // 	}

// // 	/***************************************************************/

// // 	// double **x_edge = new double*[n];
 
// // 	// for (int i = 0; i < n; i++) {
// // 	// 	x_edge[i] = new double[n];
// // 	// }

// // 	// int l = 0;
// // 	// for(int i = 0; i < n; i++) {
// // 	// 	for(int j = i+1; j < n; j++) {
// // 	// 		x_edge[i][j] = x_vals[l++];
// // 	// 	}
// // 	// }
	
// // 	// x_vals.end();

// // 	// cutSetPool = MaxBack(x_edge, n);

// // 	// /***************** Creating the constraints ***************/
// // 	// for (int c = 0; c < cutSetPool.size(); c++) {
// // 	// 	IloExpr p(getEnv());
// // 	// 	for(int i = 0; i < cutSetPool[c].size(); i++){
// // 	// 		for(int j = 0; j < cutSetPool[c].size(); j++){
// // 	// 			if(cutSetPool[c][i] < cutSetPool[c][j]){
// // 	// 				p += x[cutSetPool[c][i]][cutSetPool[c][j]];
// // 	// 			}
// // 	// 		}
// // 	// 	}
// // 	// 	int RHS = cutSetPool[c].size();
// // 	// 	cons.push_back(p <= RHS - 1);
// // 	// }
// // 	// /**********************************************************/

// // 	// /*********** Adding the constraints to the model **********/
// // 	// for(int i = 0; i < cons.size(); i++){
// // 	// 	add(cons.at(i)).end();
// // 	// }
// // 	// /**********************************************************/

// // 	// /******************* Cleaning the memory ******************/
// // 	// for (int i = 0; i < n; i++) {
// // 	// 	delete[] x_edge[i];
// // 	// }
// // 	// delete[] x_edge;
// // 	/**********************************************************/
// // }
// // /*****************************************************************************************************************/


// // bool MyLazyCallback::isCustomer(int index) {
// // 	return (index >= 0) && (index < this->nCustomers);
// // }


// // bool MyLazyCallback::isDepot(int index) {
// // 	return (index >= this->nCustomers + 2*this->nParcels) && (index < this->nCustomers + 2*this->nParcels + 2*this->v);
// // }


// // bool MyLazyCallback::isParcel(int index) {
// // 	return (index >= this->nCustomers) && (index < (this->nCustomers + 2*this->nParcels));
// // }

// // bool MyLazyCallback::isPickUpParcel(int index) {
// // 	return (index >= this->nCustomers) && (index < (this->nCustomers + this->nParcels));
// // }

// // bool MyLazyCallback::isDeliveryParcel(int index) {
// // 	return (index >= this->nCustomers + this->nParcels) && (index < (this->nCustomers + 2*this->nParcels));
// // }

// // bool MyLazyCallback::finalVerifier(vector<int> &route) {
// // 	int contP = 0;
// // 	int contD = 0;

// // 	for (int i = 0; i < route.size(); i++) {
// // 		if (isPickUpParcel(route[i])) {
// // 			contP--;
// // 		} else if (isDeliveryParcel(route[i])) {
// // 			contD--;
// // 		} else if (isCustomer(route[i])){
// // 			contP = min(1, contP + 1);
// // 			contD = min(1, contD + 1);
// // 		}

// // 		if (contP < -1 || contD < -1) {
// // 			return false;
// // 		}
// // 	}

// // 	if (contP == -1 || contD == -1) {
// // 		return false;
// // 	}

// // 	return true;
// // }

// // vector<pair<int, int>> MyLazyCallback::makePairs(vector<int> &route, int index) {
// // 	vector<pair<int, int>> pares;

// // 	pares.push_back(make_pair( route[index]-nParcels , route[index] ));
// // 	pares.push_back(make_pair( route[index+1] , route[index+2] ));
// // 	pares.push_back(make_pair( route[index+3] , route[index+3]+nParcels ));

// // 	return pares;
// // }