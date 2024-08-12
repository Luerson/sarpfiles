
#include "MyLazyCallback.h"
#include <map>
#include <list>

/********************************************** Class' Constructor **********************************************/
MyLazyCallback::MyLazyCallback(IloEnv env, const IloArray <IloArray < IloArray <IloBoolVarArray> > >& x_ref, nodeArcsStruct *arrConj, int nodes, int veic, int parcels, int customers) : IloCplex::LazyConstraintCallbackI(env), x(x_ref), x_vars(env), n(nodes), v(veic), nas(arrConj), nParcels(parcels), nCustomers(customers)
{
	int num = 0;
	/********** Filling x_vars **********/
	// for(int i = 0; i < n; i++) {
	// 	for(int j = 0; j < n; j++){
	// 		if (!nas->arcs[i][j]) {
	// 			continue;
	// 		}
	// 		for (int k1 = 0; k1 < nas->arcV[i][j].size(); k1++) {
	// 			int k = nas->arcV[i][j][k1];
	// 			x_vars.add(x[i][j][k]);
	// 		}
	// 	}
	// }
	/************************************/
	int V = nCustomers + 2*nParcels + v;

	for (int i = 0; i < V; i++){
		if (i >= nCustomers && i < nCustomers + 2*nParcels) {
			continue;
		}

        for(int j = 0; j < n; ++j){
			if (j >= nCustomers && j < V) {
				continue;
			}

            if (nas->arcs[i][j] != true){
                continue; // If arc i to j is invalid
            }

			for (int a = 0; a < nas->subsequences[make_pair(i, j)].size(); a++) {
				for(int k1 = 0; k1 < nas->arcV[i][j].size(); k1++){
					int k = nas->arcV[i][j][k1];
					x_vars.add(x[i][j][a][k]);
					// // TODO UNCOMMENT //  << "x: [" << i << "][" << j << "][" << k << "]" << endl;
				}
			}
        }
    }
}
/*****************************************************************************************************************/

/************************** Return a callback copy. This method is a CPLEX requirement ***************************/
IloCplex::CallbackI* MyLazyCallback::duplicateCallback() const 
{ 
   return new (getEnv()) MyLazyCallback(getEnv(), x, nas, n, v, nParcels, nCustomers); 
}
/*****************************************************************************************************************/

/************************************ Callback's code that is runned by CPLEX ************************************/
void MyLazyCallback::main()
{	
	// /********** Getting the relaxed variables values **********/
	// IloNumArray x_vals(getEnv(), (0.5*(n)*(n-1)));
	// getValues(x_vals, x_vars);
	// /**********************************************************/
   
	// vector <vector<int> > cutSetPool;
	// vector<IloConstraint> cons; 

	IloNumArray x_vals(getEnv());
    getValues(x_vals, x_vars);

	// Matrix to store every route in the current solution
	vector<vector<int>> routes;
	map<int, int> arcDirection;

    // Imprimir os valores das variáveis relaxadas
    int index = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
			if (!nas->arcs[i][j]) {
				// // TODO UNCOMMENT // cout << i << " " << j << endl;
				continue;
			}

			pair<int, int> myPair = make_pair(i, j);

			for (int a = 0; a < nas->subsequences[])
			for (int k1 = 0; k1 < nas->arcV[i][j].size(); k1++)
			{
				int k = nas->arcV[i][j][k1];
				if (x_vals[index] > EPSILON)
				{
					arcDirection[i] = j;
					// // TODO UNCOMMENT // cout << i << " " << j << " " << k << endl;
				}
				index++;
			}
        }
    }

    // Liberar a memória
    x_vals.end();

	// Remember: "this->v" is the numver of available vehicles
	// // TODO UNCOMMENT // cout << endl;
	// // TODO UNCOMMENT // cout << "ROTAS" << endl;
	// // TODO UNCOMMENT // cout << "numCustomers = " << this->nCustomers << endl;
	// // TODO UNCOMMENT // cout << "numParcels = " << this->nParcels << endl;
	// // TODO UNCOMMENT // cout << "numVehicles = " << this->v << endl;
	for (int i = 0; i < this->v; i++)
	{
		int startDepotIndex = this->nCustomers + 2*this->nParcels + i;	// Start depot for the current vehicle
		int finalDepotIndex = startDepotIndex + this->v;				// Final depot for the current vehicle

		// // TODO UNCOMMENT // cout << endl;
		// // TODO UNCOMMENT // cout << "i = " << i << endl;
		// // TODO UNCOMMENT // cout << "startDepotIndex = " << startDepotIndex << endl;
		// // TODO UNCOMMENT // cout << "finalDepotIndex = " << finalDepotIndex << endl;
		// // TODO UNCOMMENT // cout << "arcDirection[" << startDepotIndex << "] = " << arcDirection[startDepotIndex] << endl;
		// // TODO UNCOMMENT // cout << endl;

		routes.emplace_back(); // Creating a new route

		routes[i].push_back(startDepotIndex);
		if (arcDirection.find(startDepotIndex) != arcDirection.end())
		{
			for (int nextNode = arcDirection[startDepotIndex]; nextNode != finalDepotIndex; nextNode = arcDirection[nextNode])
			{
				// // TODO UNCOMMENT // cout << "nextNode = " << nextNode << endl;
				routes[i].push_back(nextNode);
			}
		}
		routes[i].push_back(finalDepotIndex);
	}
	// // TODO UNCOMMENT // cout << "FIM ROTAS" << endl;

	// Printing the routes
	// // TODO UNCOMMENT // cout << endl;
	// for (int i = 0; i < routes.size(); i++)
	// {
	// 	// TODO UNCOMMENT // cout << "Route " << i << ": ";
	// 	for (int j = 0; j < routes[i].size()-1; j++)
	// 	{
	// 		// TODO UNCOMMENT // cout << routes[i][j] << " -> ";
	// 	}
	// 	// TODO UNCOMMENT // cout << routes[i].back() << endl;
	// }
	// getchar();


	/**************** Creating the Lazy Co	// // TODO UNCOMMENT // cout << "aqui 2" << endl;
			for (int j = 0; j < routes[i].size() - 1; j++) {
				expr += x[ routes[i][j] ][ routes[i][j+1] ][ i ];
			}

			// // TODO UNCOMMENT // cout << "aqui 3" << endl;

			int maxArcs = routes[i].size() - 2;
			add(expr <= maxArcs);
		}nstraints ****************/
	// // TODO UNCOMMENT // cout << "SIZES: " << x.getSize() << " " << x[0].getSize() << endl;

	// for (int i = 0; i < routes.size(); i++)
	// {
	// 	int counter = 0;

	// 	for (int j = 0; j < routes[i].size(); j++)
	// 	{
	// 		if (isParcel(routes[i][j]))
	// 		{
	// 			counter++;

	// 			if (counter >= 4 && isDeliveryParcel(routes[i][j-3]) && isPickUpParcel(routes[i][j]))
	// 			{
	// 				// Prohibiting the sequence for all vehicles
	// 				// // TODO UNCOMMENT // cout << endl;
	// 				// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
	// 				IloExpr expr(getEnv());
	// 				for (int vehicle = 0; vehicle < this->v; vehicle++)
	// 				{
	// 					// int k = arrConj->arcV[i][j][vehicle];

	// 					// TODO UNCOMMENT // cout << endl;
	// 					// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

	// 					// // // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
	// 					// expr += x[ routes[i][j-5] ][ routes[i][j-4] ][ vehicle ];

	// 					// // // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
	// 					// expr += x[ routes[i][j-4] ][ routes[i][j-3] ][ vehicle ];

	// 					// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
	// 					expr += x[ routes[i][j-3] ][ routes[i][j-2] ][ vehicle ];

	// 					// // TODO UNCOMMENT // cout << "x(" << routes[i][j-2] << "," << routes[i][j-1] << "," << vehicle << ")" << endl;
	// 					expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

	// 					// // TODO UNCOMMENT // cout << "x(" << routes[i][j-1] << "," << routes[i][j] << "," << vehicle << ")" << endl;
	// 					expr += x[ routes[i][j-1] ][ routes[i][ j ] ][ vehicle ];

	// 					// // TODO UNCOMMENT // cout << "END EXPRESSION!!" << endl;

	// 					// // TODO UNCOMMENT // cout << "Termos da expressão:" << endl;
	// 					// for (IloExpr::LinearIterator it(expr); it.ok(); ++it) {
	// 					// 	// TODO UNCOMMENT // cout << it.getVar() << " * " << it.getCoef() << endl;
	// 					// }
	// 					// for (IloExpr::QuadIterator it(expr); it.ok(); ++it) {
	// 					// 	// TODO UNCOMMENT // cout << it.getVar1() << " * " << it.getVar2() << " * " << it.getCoef() << endl;
	// 					// }
	// 					// // TODO UNCOMMENT // cout << "LAZY CONS" << endl;
	// 					// // TODO UNCOMMENT // cout << "END LAZY CONS" << endl;
	// 				}
	// 				add(expr <= 2);
	// 				// // TODO UNCOMMENT // cout << endl;
	// 				// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
	// 			}
	// 		}
	// 		else // isCustomer = true
	// 		{
	// 			counter = 0;
	// 		}
	// 	}
	// }

	// Prohibiting sequence V-P-D-P-D
	// for (int i = 0; i < routes.size(); i++)
	// {
	// 	int counter = 0;

	// 	for (int j = 0; j < routes[i].size() && j < 4; j++)
	// 	{
	// 		if (isParcel(routes[i][j]))
	// 		{
	// 			counter++;

	// 			if (counter >= 3)
	// 			{
	// 				// Prohibiting the sequence for all vehicles
	// 				// // TODO UNCOMMENT // cout << endl;
	// 				// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
	// 				for (int vehicle = 0; vehicle < this->v; vehicle++)
	// 				{
	// 					// int k = arrConj->arcV[i][j][vehicle];
	// 					// // TODO UNCOMMENT // cout << endl;
	// 					// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;
	// 					if (arrConj->arcs[ vehicle ][ routes[i][j-2] ]) {
	// 						IloExpr expr(getEnv()); //!! Isso estava fora do "if". Eu movi para dentro.
	// 												//!! Mas por que não colocar fora do laço de repetição?

	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
	// 						expr += x[ routes[vehicle][0] ][ routes[i][j-2] ][ vehicle ];

	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
	// 						expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
	// 						expr += x[ routes[i][j-1] ][ routes[i][j] ][ vehicle ];

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
	// 						add(expr <= 2); //!! Isso estava fora do "if". Eu movi para dentro
	// 					}
	// 				}
	// 				// // TODO UNCOMMENT // cout << endl;
	// 				// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
	// 			}
	// 		}
	// 		else // isCustomer = true
	// 		{
	// 			counter = 0;
	// 		}
	// 	}
	// }

	// for (int i = 0; i < routes.size(); i++)
	// {
	// 	int counter = 0;

	// 	for (int j = routes[i].size() - 4; j < routes[i].size(); j++)
	// 	{
	// 		if (isParcel(routes[i][j]))
	// 		{
	// 			counter++;

	// 			if (counter >= 3)
	// 			{
	// 				// Prohibiting the sequence for all vehicles
	// 				// // TODO UNCOMMENT // cout << endl;
	// 				// // TODO UNCOMMENT // cout << "VEHICLES!" << endl;
	// 				for (int vehicle = 0; vehicle < this->v; vehicle++)
	// 				{
	// 					// int k = arrConj->arcV[i][j][vehicle];
	// 					// // TODO UNCOMMENT // cout << endl;
	// 					// // TODO UNCOMMENT // cout << "EXPRESSION!!" << endl;

	// 					if (arrConj->arcs[ routes[i][j] ][ routes[vehicle].back() ]) {
	// 						IloExpr expr(getEnv()); //!! Isso estava fora do "if". Eu movi para dentro.
	// 												//!! Mas por que não colocar fora do laço de repetição?

	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-4] << "," << routes[i][j-3] << "," << vehicle << ")" << endl;
	// 						expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ vehicle ];

	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-3] << "," << routes[i][j-2] << "," << vehicle << ")" << endl;
	// 						expr += x[ routes[i][j-1] ][ routes[i][j] ][ vehicle ];

	// 						// // TODO UNCOMMENT // cout << "x(" << routes[i][j-5] << "," << routes[i][j-4] << "," << vehicle << ")" << endl;
	// 						expr += x[ routes[i][j] ][ routes[vehicle].back() ][ vehicle ];

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
	// 						add(expr <= 2); //!! Isso estava fora do "if". Eu movi para dentro
	// 					}
	// 				}
	// 				// // TODO UNCOMMENT // cout << endl;
	// 				// // TODO UNCOMMENT // cout << "END VEHICLES!" << endl;
	// 			}
	// 		}
	// 		else // isCustomer = true
	// 		{
	// 			counter = 0;
	// 		}
	// 	}
	// }

	// bool restricted = false;

	// for (int i = 0; i < routes.size(); i++)
	// {
	// 	// // TODO UNCOMMENT // cout << "aqui" << endl;
	// 	for (int j = 1; j < routes[i].size(); j++) {
	// 		if (j + 1 < routes[i].size() && isPickUpParcel(j) && isDeliveryParcel(j + 1)) {
	// 			int k;
	// 			for (k = j+1; k + 3 < routes[i].size();) {
	// 				if (isCustomer(k+1) && isPickUpParcel(k+2) && isDeliveryParcel(k+3)) {
	// 					k += 3;
	// 				} else {
	// 					break;
	// 				}
	// 			}

	// 			if (!(isCustomer(j-1) || isCustomer(k + 1))) {
	// 				IloExpr expr(getEnv());

	// 				for (int vehicle = 0; vehicle < this->v; vehicle++) {
	// 					for (int k1 = j; k1 < k; k1++) {
	// 						if (arrConj->arcs[k1][k1+1]) expr += x[k1][k1+1][vehicle];
	// 					}

	// 					expr += x[j-1][j][vehicle];
	// 					expr += x[k][k+1][vehicle];

	// 					for (int k1 = 0; k1 < this->n; k1++) {
	// 						if (isDeliveryParcel(k1) && arrConj->arcs[k1][j]) {
	// 							expr += x[k1][j][vehicle];
	// 						}

	// 						if (isPickUpParcel(k1) && arrConj->arcs[k][k1]) {
	// 							expr += x[k][k1][vehicle];
	// 						}
	// 					}
	// 				}

	// 				add(expr <= k - j);
	// 				restricted = true;
	// 			}
				
	// 			// j = k+1;
	// 		}
	// 	}
	// }


	// Lazy Constraint
	// In the case P - D - c - P - D, at least one of the parcels must be connected to
	// another customer

	// for (int i = 0; i < routes.size(); i++)
	// {
	// 	int counter = 0;
	
	// 	for (int j = 0; j < routes[i].size(); j++)
	// 	{
	// 		if (isPickUpParcel(routes[i][j]))
	// 		{
	// 			if ((counter == 0 || counter == 2) && isDeliveryParcel(routes[i][j+1]))	counter++;
	// 			else																	counter = 0;

	// 			// If counter == 3, then we found a pattern P - D - c - P - D
	// 			if (counter == 3)
	// 			{
	// 				if (!isCustomer(routes[i][j-4]) && !isCustomer(routes[i][j+2]))
	// 				{
	// 					// TODO UNCOMMENT // cout 	<< "EXPRESSAO P - D - c - P - D: "
	// 							<< routes[i][j-3] << " " << routes[i][j-2] << " " << routes[i][j-1] << " "
	// 							<< routes[i][j] << " " << routes[i][j+1] << endl;

	// 					IloExpr expr(getEnv());

	// 					for (int k = 0; k < this->v; k++)
	// 					{
	// 						expr += x[ routes[i][j-3] ][ routes[i][j-2] ][ k ];
	// 						expr += x[ routes[i][j-2] ][ routes[i][j-1] ][ k ];
	// 						expr += x[ routes[i][j-1] ][ routes[i][j] ][ k ];
	// 						expr += x[ routes[i][j] ][ routes[i][j+1] ][ k ];
	// 					}

	// 					for (int customer = 0; customer < nCustomers; customer++)
	// 					{
	// 						for (int k = 0; k < this->v; k++)
	// 						{
	// 							if (arrConj->arcs[customer][routes[i][j-3]]) expr -= x[ customer ][ routes[i][j-3] ][ k ];
	// 							if (arrConj->arcs[routes[i][j+1]][customer]) expr -= x[ routes[i][j+1] ][ customer ][ k ];
	// 						}
	// 					}

	// 					////////////////////////////////
	// 					add(expr <= 3);
	// 				}
					
	// 				////////////////////////
	// 				counter = 1;
	// 			}
	// 		}
	// 		else if (isCustomer(routes[i][j]))
	// 		{
	// 			if (counter == 1) 	counter++;
	// 			else				counter = 0;
	// 		}
	// 	}
	// }


	// for (int k = 0; k < routes.size(); k++) {
	// 	for (int i = 0; i < routes[k].size(); i++) {
	// 		cout << routes[k][i] << " ";
	// 	}
	// 	cout << endl;
	// }
	// getchar();

	// Creating restrictions for every subsequence that violates the formulation
	// NO PERMUTATION!

	// for (int r = 0; r < routes.size(); r++)
	// {
	// 	for (int size = 2; size <= (routes[r].size()-2); size++)
	// 	{
	// 		for (int startIndex = 1; startIndex < (routes[r].size()-1); startIndex++)
	// 		{
	// 			int lastIndex = startIndex + size - 1;
	// 			if (lastIndex >= routes[r].size()-1)
	// 			{
	// 				break;
	// 			}
	// 			// getchar();

	// 			// If it's a sequence of type P - D - ... - P - D, or simply P - D
	// 			if (isPickUpParcel(routes[r][startIndex]) && isDeliveryParcel(routes[r][startIndex+1]) && isPickUpParcel(routes[r][lastIndex-1]) && isDeliveryParcel(routes[r][lastIndex]))
	// 			{
	// 				int counter	= 0;

	// 				for (int i = startIndex; i <= lastIndex; i++)
	// 				{
	// 					if (isPickUpParcel(routes[r][i]))	counter--;
	// 					else if (isCustomer(routes[r][i]))	counter++;
	// 				}

	// 				if ((counter + isCustomer(routes[r][startIndex-1]) + isCustomer(routes[r][lastIndex+1])) < 0)
	// 				{
	// 					// Creating the restriction
	// 					IloExpr expr(getEnv());

	// 					for (int currNode = startIndex; currNode <= lastIndex; currNode++)
	// 					{
	// 						for (int nextNode = startIndex; nextNode <= lastIndex; nextNode++) {
	// 							if (!arrConj->arcs[ routes[r][currNode] ][ routes[r][nextNode] ]) continue;

	// 							for (int k1 = 0; k1 < arrConj->arcV[ routes[r][currNode] ][ routes[r][nextNode] ].size(); k1++)
	// 							{
	// 								int k = arrConj->arcV[ routes[r][currNode] ][ routes[r][nextNode] ][ k1 ];
	// 								expr += x[ routes[r][currNode] ][ routes[r][nextNode] ][ k ];
	// 							}
	// 						}
	// 						// cout << currNode << " ";
	// 					}

						
	// 					// int start;

	// 					// if (routes[r][startIndex - 1] < nCustomers) start = 0;
	// 					// else start = nCustomers;

	// 					// for (int left = start; left < n; left++) {
	// 					// 	if (!arrConj->arcs[ left ][ routes[r][startIndex] ]) continue;

	// 					// 	for (int k1 = 0; k1 < arrConj->arcV[ left ][ routes[r][startIndex] ].size(); k1++)
	// 					// 	{
	// 					// 		int k = arrConj->arcV[ left ][ routes[r][startIndex] ][ k1 ];
	// 					// 		expr += x[ left ][ routes[r][startIndex] ][ k ];
	// 					// 	}
	// 					// }

	// 					// if (routes[r][lastIndex + 1] < nCustomers) start = 0;
	// 					// else start = nCustomers;
						
	// 					// for (int right = start; right < n; right++) {
	// 					// 	if (!arrConj->arcs[ routes[r][lastIndex] ][ right ]) continue;

	// 					// 	for (int k1 = 0; k1 < arrConj->arcV[ routes[r][lastIndex] ][ right ].size(); k1++)
	// 					// 	{
	// 					// 		int k = arrConj->arcV[ routes[r][lastIndex] ][ right ][ k1 ];
	// 					// 		expr += x[ routes[r][lastIndex] ][ right ][ k ];
	// 					// 	}
	// 					// }

	// 					if (arrConj->arcs[ routes[r][startIndex - 1] ][ routes[r][startIndex] ]) {
	// 						for (int k1 = 0; k1 < arrConj->arcV[ routes[r][startIndex-1] ][ routes[r][startIndex] ].size(); k1++)
	// 						{
	// 							int k = arrConj->arcV[ routes[r][startIndex - 1] ][ routes[r][startIndex] ][ k1 ];
	// 							expr += x[ routes[r][startIndex - 1] ][ routes[r][startIndex] ][ k ];
	// 						}
	// 					}

	// 					if (arrConj->arcs[ routes[r][lastIndex] ][ routes[r][lastIndex + 1] ]) {
	// 						for (int k1 = 0; k1 < arrConj->arcV[ routes[r][lastIndex] ][ routes[r][lastIndex + 1] ].size(); k1++)
	// 						{
	// 							int k = arrConj->arcV[ routes[r][lastIndex] ][ routes[r][lastIndex + 1] ][ k1 ];
	// 							expr += x[ routes[r][lastIndex] ][ routes[r][lastIndex + 1] ][ k ];
	// 						}
	// 					}
						
	// 					// cout << lastIndex << endl;
	// 					// getchar();

	// 					// for (int c = 0; c < nCustomers; c++)
	// 					// {
	// 					// 	if (arrConj->arcs[ c ][ routes[r][startIndex] ])
	// 					// 	{
	// 					// 		for (int k1 = 0; k1 < arrConj->arcV[ c ][ routes[r][startIndex] ].size(); k1++)
	// 					// 		{
	// 					// 			int k = arrConj->arcV[ c ][ routes[r][startIndex] ][ k1 ];
	// 					// 			expr -= x[ c ][ routes[r][startIndex] ][ k ];
	// 					// 		}
	// 					// 	}

	// 					// 	if (arrConj->arcs[ routes[r][lastIndex] ][ c ])
	// 					// 	{
	// 					// 		for (int k1 = 0; k1 < arrConj->arcV[ routes[r][lastIndex] ][ c ].size(); k1++)
	// 					// 		{
	// 					// 			int k = arrConj->arcV[ routes[r][lastIndex] ][ c ][ k1 ];
	// 					// 			expr -= x[ routes[r][lastIndex] ][ c ][ k ];
	// 					// 		}
	// 					// 	}
	// 					// }

	// 					expr -= size;
	// 					// expr -= counter;

	// 					// Adding the constraint
	// 					add(expr <= 0);

	// 					// Printing the sebsequence
	// 					cout << endl;
	// 					cout << "route: ";
	// 					for (int i = 0; i < routes[r].size(); i++)
	// 					{
	// 						cout << routes[r][i] << " ";
	// 					}
	// 					cout << endl;

	// 					cout << endl;
	// 					cout << "sequence: ";
	// 					for (int i = startIndex; i <= lastIndex; i++)
	// 					{
	// 						cout << routes[r][i] << " ";
	// 					}
	// 					cout << endl;
	// 					cout << "sequence: ";
	// 					for (int i = startIndex; i <= lastIndex; i++)
	// 					{
	// 						if (isCustomer(routes[r][i])) 		cout << "c";
	// 						if (isPickUpParcel(routes[r][i])) 	cout << "P";
	// 						if (isDeliveryParcel(routes[r][i])) cout << "D";

	// 						if (i != lastIndex)	cout << " - ";
	// 					}
	// 					cout << endl;
	// 					cout << endl;
	// 					// getchar();
	// 				}
	// 			}
	// 		}
	// 	}
	// }






	/***************************************************************/

	// double **x_edge = new double*[n];
 
	// for (int i = 0; i < n; i++) {
	// 	x_edge[i] = new double[n];
	// }

	// int l = 0;
	// for(int i = 0; i < n; i++) {
	// 	for(int j = i+1; j < n; j++) {
	// 		x_edge[i][j] = x_vals[l++];
	// 	}
	// }
	
	// x_vals.end();

	// cutSetPool = MaxBack(x_edge, n);

	// /***************** Creating the constraints ***************/
	// for (int c = 0; c < cutSetPool.size(); c++) {
	// 	IloExpr p(getEnv());
	// 	for(int i = 0; i < cutSetPool[c].size(); i++){
	// 		for(int j = 0; j < cutSetPool[c].size(); j++){
	// 			if(cutSetPool[c][i] < cutSetPool[c][j]){
	// 				p += x[cutSetPool[c][i]][cutSetPool[c][j]];
	// 			}
	// 		}
	// 	}
	// 	int RHS = cutSetPool[c].size();
	// 	cons.push_back(p <= RHS - 1);
	// }
	// /**********************************************************/

	// /*********** Adding the constraints to the model **********/
	// for(int i = 0; i < cons.size(); i++){
	// 	add(cons.at(i)).end();
	// }
	// /**********************************************************/

	// /******************* Cleaning the memory ******************/
	// for (int i = 0; i < n; i++) {
	// 	delete[] x_edge[i];
	// }
	// delete[] x_edge;
	/**********************************************************/
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