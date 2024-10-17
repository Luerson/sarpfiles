
#include "MylazyNode.h"
#include <map>
#include <list>

/********************************************** Class' Constructor **********************************************/
MylazyNode::MylazyNode(IloEnv env, const IloArray<IloArray<IloBoolVarArray>>& x_ref, nodeArcsStruct *nas, instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, int nodes, int veic, int parcels, int customers) : IloCplex::LazyConstraintCallbackI(env), x(x_ref), x_vars(env), n(nodes), v(veic), arrConj(nas), inst(inst), nodeVec(nodeVec), mdist(mdist), nParcels(parcels), nCustomers(customers)
{
	int num = 0;
	/********** Filling x_vars **********/
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++){
			if (!nas->arcs[i][j]) {
				continue;
			}
			for (int k1 = 0; k1 < nas->arcV[i][j].size(); k1++) {
				int k = nas->arcV[i][j][k1];
				x_vars.add(x[i][j][k]);
			}
		}
	}
	/************************************/
}
/*****************************************************************************************************************/

/************************** Return a callback copy. This method is a CPLEX requirement ***************************/
IloCplex::CallbackI* MylazyNode::duplicateCallback() const 
{ 
   return new (getEnv()) MylazyNode(getEnv(), x, arrConj, inst, nodeVec, mdist, n, v, nParcels, nCustomers); 
}
/*****************************************************************************************************************/

/************************************ Callback's code that is runned by CPLEX ************************************/
void MylazyNode::main()
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
	vector<bool> used(n, false);
	vector<vector<int>> subtours;

    // Imprimir os valores das variáveis relaxadas
    int index = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
			if (!arrConj->arcs[i][j]) {
				continue;
			}

			for (int k1 = 0; k1 < arrConj->arcV[i][j].size(); k1++)
			{
				int k = arrConj->arcV[i][j][k1];
				if (x_vals[index] > 0.5)
				{	
					arcDirection[i] = j;
				}
				index++;
			}
        }
    }

    x_vals.end();

	// Getting the routes made by the solver
	for (int i = 0; i < this->v; i++)
	{
		int startDepotIndex = this->nCustomers + 2*this->nParcels + i;	// Start depot for the current vehicle
		int finalDepotIndex = startDepotIndex + this->v;				// Final depot for the current vehicle

		routes.emplace_back(); // Creating a new route

		routes[i].push_back(startDepotIndex);
		used[startDepotIndex] = true;
		if (arcDirection.find(startDepotIndex) != arcDirection.end())
		{
			for (int nextNode = arcDirection[startDepotIndex]; nextNode != finalDepotIndex; nextNode = arcDirection[nextNode])
			{
				routes[i].push_back(nextNode);
				used[nextNode] = true;
			}
		}
		routes[i].push_back(finalDepotIndex);
		used[finalDepotIndex] = true;
	}

	/********** Subtour constraints **********/
	// Getting the subtours that have no connection to the depot
	for (int i = 0; i < n; i++)
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

	// Creating subtour constraints
	for (int k = 0; k < subtours.size(); k++) {
		IloExpr expr(getEnv());

		for (int i = 0; i < subtours[k].size(); i++) {
			for (int j = 0; j < subtours[k].size(); j++) {
				
				int u  = subtours[k][i];
				int v = subtours[k][j];
				
				for (int k1 = 0; k1 < arrConj->arcV[u][v].size(); k1++) {
					int k2 = arrConj->arcV[u][v][k1];

					expr += x[u][v][k2];
				}
			}
		}

		add(expr <= (int)subtours[k].size() - 1);
	}

	/********** duration constraints **********/

	vector<vector<int>> ToPrevent;

	/********** The time limit between customers shouldn't be exceeded **********/
	for (int k = 0; k < routes.size(); k++) {
		vector<int> sequence;
		int lastRelevant = routes[k][0];
		float time = 0;

		sequence.push_back(lastRelevant);
		time += nodeVec[lastRelevant].delta;

		for (int i = 1; i < routes[k].size(); i++) {
			int u = routes[k][i-1];
			int v = routes[k][i];

			bool isRelevant = ((v == routes[k].back()) || (v < nCustomers));

			sequence.push_back(v);

			time += mdist[u][v]/inst->vmed;

			if (isRelevant) {

				float timeStart = nodeVec[lastRelevant].e;
				float timeEnd = nodeVec[v].l;

				if (time > timeEnd - timeStart) {
					ToPrevent.push_back(sequence);
				}

				sequence.clear();
				time = 0;
				sequence.push_back(v);
				lastRelevant = v;
			}
			time += nodeVec[v].delta;
		}
	}

	for (int k = 0; k < ToPrevent.size(); k++)
	{
		// cout << "DURATION CONSTRAINT " << s << endl;
		IloExpr expr(getEnv());

		const int numEdges = ToPrevent[k].size() - 1;

		/********** The actual expressions/lazy constraints **********/
		for (int i = 0; i < ToPrevent[k].size() - 1; i++) {
			int u = ToPrevent[k][i];
			int v = ToPrevent[k][i+1];

			for (int k1 = 0; k1 < arrConj->arcV[u][v].size(); k1++) {
				int k2 = arrConj->arcV[u][v][k1];

				expr += x[u][v][k2];
			}
		}

		add(expr <= numEdges - 1);
	}
	ToPrevent.clear();

	/********** precedence constraints **********/

	/********** Adding sequences of type "P - d - ... - P" to ToPrevent **********/
	for (int k = 0; k < routes.size(); k++) {
		vector<int> sequence;
		bool restricted = false;

		int currPickup = -1;

		for (int i = 0; i < routes[k].size(); i++) {
			int h = routes[k][i];

			bool isDelivery = (h >= inst->n + inst->m && h < inst->n + 2*inst->m);
			bool isPickup = (h >= inst->n && h < inst->n + inst->m);
			bool isDepot = (h >= inst->n + 2*inst->m);

			if ((currPickup != -1) && (h >= inst->n) && (h != currPickup + inst->m)) {
				sequence.push_back(h);
				ToPrevent.push_back(sequence);

				sequence.clear();
				currPickup = -1;
			}

			if (h == currPickup + inst->m) {
				currPickup = -1;
				sequence.clear();
			}
			
			if (isPickup) {
				currPickup = h;
				sequence.clear();
			}

			sequence.push_back(h);
		}
	}

	for (int k = 0; k < ToPrevent.size(); k++)
	{
		// cout << "DURATION CONSTRAINT " << s << endl;
		IloExpr expr(getEnv());

		const int numEdges = ToPrevent[k].size() - 1;

		/********** The actual expressions/lazy constraints **********/
		for (int i = 0; i < ToPrevent[k].size() - 1; i++) {
			int u = ToPrevent[k][i];
			int v = ToPrevent[k][i+1];

			for (int k1 = 0; k1 < arrConj->arcV[u][v].size(); k1++) {
				int k2 = arrConj->arcV[u][v][k1];

				expr += x[u][v][k2];
			}
		}

		add(expr <= numEdges - 1);
	}
	ToPrevent.clear();
}
/*****************************************************************************************************************/


bool MylazyNode::isCustomer(int index) {
	return (index >= 0) && (index < this->nCustomers);
}


bool MylazyNode::isDepot(int index) {
	return (index >= this->nCustomers + 2*this->nParcels) && (index < this->nCustomers + 2*this->nParcels + 2*this->v);
}

bool MylazyNode::isParcel(int index) {
	return (index >= this->nCustomers) && (index < (this->nCustomers + 2*this->nParcels));
}

bool MylazyNode::isPickUpParcel(int index) {
	return (index >= this->nCustomers) && (index < (this->nCustomers + this->nParcels));
}

bool MylazyNode::isDeliveryParcel(int index) {
	return (index >= this->nCustomers + this->nParcels) && (index < (this->nCustomers + 2*this->nParcels));
}

bool MylazyNode::finalVerifier(vector<int> &route) {
	int contP = 0;
	int contD = 0;

	for (int i = 0; i < route.size(); i++) {
		if (isPickUpParcel(route[i])) {
			contP--;
		} else if (isDeliveryParcel(route[i])) {
			contD--;
		} else if (isCustomer(route[i])){
			contP = min(1, contP + 1);
			contD = min(1, contD + 1);
		}

		if (contP < -1 || contD < -1) {
			return false;
		}
	}

	if (contP == -1 || contD == -1) {
		return false;
	}

	return true;
}