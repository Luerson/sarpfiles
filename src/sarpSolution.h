#ifndef SARPSOLUTION_H
#define SARPSOLUTION_H

#include <algorithm>

#include "sarpRoute.h"
#include "readdata.h"
using namespace std;

class sarpSolution {
   protected:

    vector<sarpRoute> routes;
    double cost;

    vector<nodeStat> unservParc;

    int usedK;
    // vector<vector<vector<bool>>> pairNeighbStatus;

    // void createRoutesVector();
    // void updateCost();
    // void updateRoutes();

   public:
    
    sarpSolution() { this->cost = INT32_MAX; this->usedK = 0; }

    // Solution(Data *data);
    // Solution(Data *data, vector<Route *> routesVec);
    // Solution(const Solution &other);
    ~sarpSolution(){};

    inline double getCost() { return cost; }
    vector<sarpRoute> getRoutes() { return routes; }
    sarpRoute getRoute(int idx) { return routes[idx]; }

    inline int getvehicle() { return this->usedK; }
    
    inline int getRoutesSize() const { return this->routes.size(); }

    void addRoute(sarpRoute* route);
    
    void updateRoutes(sarpRoute *route,  int idr);

    void addCost(double delta);

    // access routes
    vector<sarpRoute>::iterator begin() { return routes.begin(); }
    vector<sarpRoute>::const_iterator begin() const { return routes.begin(); }
    vector<sarpRoute>::iterator end() { return routes.end(); }
    vector<sarpRoute>::const_iterator end() const { return routes.end(); }
    
    void updateCost();
    void updateVehicles();

    void printSol(instanceStat *inst);
    // sarpRoute *getRoute(int idx) { return routes[idx]; }
    // inline int getRoutesSize() const { return this->routes.size(); }s
};

#endif