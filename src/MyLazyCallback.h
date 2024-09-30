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
#include "bundleData.h"

#define EPSILON 1e-6

/**************************** Lazy callback ****************************/ 
class MyLazyCallback : public IloCplex::LazyConstraintCallbackI 
{
     public:
 
          //Class constructor
          MyLazyCallback(IloEnv env, const IloArray <IloArray <IloBoolVarArray> >& x_ref, bundleStat *bStat, double **mdist, instanceStat *inst, probStat* problem, vector<nodeStat> &nodeVec, int nodes, int veic, int parcels, int customers);
          
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
          IloArray <IloArray <IloBoolVarArray> > x;

          //x_vars contains all the x variables values. A more fast way to get the variables values from CPLEX 
          IloNumVarArray x_vars;

          //Problem's dimension
          int n;
          int v;
          int nCustomers;
          int nParcels;
          vector<nodeStat> &nodeVec;
          bundleStat* bStat;
          instanceStat *inst;
          probStat* problem;
          double **mdist;

          /********** Auxiliary Vectors to easily identify types of bundles **********/
          vector<bool> pickupBundles;   // Bundles of type "P - d"
          vector<bool> deliveryBundles; // Bundles of type "d - D"
          vector<bool> pdBundles;       // pickup-delivery bundles
          vector<bool> customerBundles; // Bundles with only the customer alone
          vector<bool> pcdBundles;      // Bundles of type "P - d - D"
          vector<bool> depotBundles;    // bundles composed by a single depot (either start or ending)

          /********** Auxiliary Vectors to easily iterate through types of bundles **********/
          vector<int> pickupBundles_Iter;    // Bundles of type "P - d"
          vector<int> deliveryBundles_Iter;  // Bundles of type "d - D"
          vector<int> pdBundles_Iter;        // pickup-delivery bundles
          vector<int> customerBundles_Iter;  // Bundles with only the customer alone
          vector<int> pcdBundles_Iter;       // Bundles of type "P - d - D"
          vector<int> depotBundles_Iter;     // bundles composed by a single depot (either start or ending)


          // Temporary variables
          // double *bestSolVal;
};
/***********************************************************************/

#endif