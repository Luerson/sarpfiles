/***************************************************/
/* Functions prototypes by Prof. Anand Subramanian */
/***************************************************/

#ifndef MY_LAZY_CALLBACK_H
#define MY_LAZY_CALLBACK_H

#include <ilcplex/ilocplex.h>
#include <set>
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits.h>
#include <stdlib.h>
#include <stack>
#include <algorithm>
#include <exception>
#include "SarpADS.h"

#define EPSILON 1e-6

/**************************** Lazy callback ****************************/ 
class MyLazyCallback : public IloCplex::LazyConstraintCallbackI 
{
   public:
 
        //Class constructor
        MyLazyCallback(IloEnv env, const IloArray <IloArray < IloArray <IloBoolVarArray> > >& x_ref, nodeArcsStruct *nas, double **mdist, instanceStat *inst, probStat* problem, vector<nodeStat> &nodeVec, int nodes, int veic, int parcels, int customers);
        
        //Method that returned a callback's copy. CPLEX requirement
        IloCplex::CallbackI* duplicateCallback() const;
        
        //Callback's code that is runned by CPLEX
        void main();

        /**************** Auxiliary Methods ****************/
        //   bool isCustomer(int);

        //   bool isDepot(int);

        //   bool isParcel(int);

        //   bool isPickUpParcel(int);

        //   bool isDeliveryParcel(int);

        //   bool finalVerifier(vector<int> &);

        //   vector<pair<int, int>> makePairs(vector<int> &, int);

   private:
	
        //Used to allocate x into x_vars, when the object is created
        IloArray <IloArray < IloArray <IloBoolVarArray> > > x;

        //x_vars contains all the x variables values. A more fast way to get the variables values from CPLEX 
        IloNumVarArray x_vars;

        //Problem's dimension
        int n;
        int v;
        int nCustomers;
        int nParcels;
        vector<nodeStat> &nodeVec;
        nodeArcsStruct* nas;
        instanceStat *inst;
        probStat* problem;
        double **mdist;
};
/***********************************************************************/

#endif