#ifndef SARPROUTE_H_
#define SARPROUTE_H_

#include <vector>
#include <algorithm>

#include "readdata.h"
#include "Statistics.h"

class sarpBlock {

protected:
    double profit_;
    // double cost_;

    vector<int> block;

    double starttime;
    double endtime;

    //Endpos: Point of ending the loop for building the block (the last request in 
    //the block is in position endPos-1). 
    int iniPos, endPos;

    int lastPass, lastPassPos, firstPass, firstPassPos;


public:

    double profit() const{ return profit_; };
    int getBlockReq(int pos) { return block[pos]; };
    int getLastReq() { return block.back(); };
    int getiniPos () const{ return iniPos; };
    int getendPos () const{ return endPos; };
    
    int getLastPass() { return lastPassPos; }; //returns the position of blocks last passenger
    int getFirstPass() { return firstPassPos; }; //returns a position of blocks first passenger

    double getStart() { return starttime; };
    double getEnd() { return endtime; };

    int getBlockSize() { return block.size(); };

    vector<int> getBlock() { return block; };

    void makeBlock(vector<int> oriRoute, int starting, int ending);
    

    void clearBlock() { block.clear(); };

    void blockProfit(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist);
    void calcBlockTimes(instanceStat *inst, 
                        vector<nodeStat> &nodeVec, double **Mdist);
};

class sarpRoute {
protected:

    // the cost of this route
    double cost_;
    //start time
    double starttime;
    //end time
    double endtime;
    //first and last passenger requests
    int firstPass;
    int lastPass;
    //first and last passenger request positions
    int firstPassPos;
    int lastPassPos;
    //route id
    int id;

    vector <pair <nodeStat, int> > passandpos; //passengers and their positions in the route
    vector < pair <int, int> > pdvec; //first: position of parcel pickup; second: position of parcel delivery in the route.
    vector <pair <int, int> > loadofroute; //first:actual load; second:passengers visited while carrying parcels.


public:
    // the route
    vector<int> nodes_;
    Runtime stats;

    // default constructor: route with only starting and ending depot visits
    sarpRoute(instanceStat *inst, int vehicle);
    sarpRoute(const sarpRoute *rt);

    // access route elements
    double cost() const { return cost_; };
    double startTime() const { return starttime; };
    double endTime() const { return endtime; };
    int fPass() const { return firstPass; };
    int lPass() const { return lastPass; };
    
    void printTotalTime();

    void clearRoute();

    inline int getNodesSize() { return nodes_.size(); };
    inline int getLoadSize() { return loadofroute.size(); };
    inline int getPassPosSize() { return passandpos.size(); };
    inline int getPaPdvecSize() { return pdvec.size(); };

    pair <nodeStat, int> getPassReq(int i) { return passandpos[i]; } ;
    pair <int, int> getPDReq(int i) { return pdvec[i]; };

    int firstpos() const { return nodes_[1]; };
    int lastpos() const { return nodes_[nodes_.size()-2]; };
    int getReq(int pos) const { return nodes_[pos]; };

    int getReqLoad(int pos) const { return loadofroute[pos].first; };

    int getPU(int request) const { return pdvec[request].first; };
    int getDL(int request) const { return pdvec[request].second; };

    int getNextPass(int request);
    int getPrevPass(int request);

    inline int getId() { return id; };

    double getProfit (vector<nodeStat> &nodeVec, int pos) { return nodeVec[nodes_[pos]].profit; };
    
    void updateCost(double delta) { cost_ += delta; }

    void clearPDVec();

    vector<int>::iterator begin() { return nodes_.begin(); };
    vector<int>::iterator end() { return nodes_.end(); };
    vector<int> getNodes() {return nodes_; };

    bool testInsertion(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist, int position, int request);
    bool testInsertionParcel(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist, int pos1, int pos2, int pu, int dl);

    bool testSwap(instanceStat *inst, vector<nodeStat> &nodeVec,
                         double **Mdist,int pos1, int pos2, 
                         pair <int, int> inter1, pair <int, int> inter2);
    bool testRelocate(instanceStat *inst, vector<nodeStat> &nodeVec,
                            double **Mdist,int pos1, int pos2, 
                            pair <int, int> inter1);
    void calcCost(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist);

    //returns the new endtime without the block
    double blockrmvTime(instanceStat *inst, 
                        vector<nodeStat> &nodeVec, 
                        double **Mdist,
                        int iniPos);    

    bool testBlockIns(instanceStat *inst, 
                        vector<nodeStat> &nodeVec, 
                        double **Mdist, double &newEnd,
                        int strPos, int endPos, sarpBlock newBlock);

    //only when the first insertion is a passenger request.
    bool fInsertion(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist, int request);

    // int load() const { return load_; };
    bool fInsertionParcel(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist, int pu, int dl);
    //update starting and end times after insertion or erasing.
    void updateTimes(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist);


    //updates passenger list with their positions.
    void updatePass(instanceStat *inst, vector<nodeStat> &nodeVec);

    void updateLoad(instanceStat *inst, vector<nodeStat> &nodeVec);
    //determine available positions to add new passenger request 
    void availablePos(instanceStat *inst, vector<nodeStat> &nodeVec, int request, probStat* problem, vector<int> &positions);

    void updateParcels(int request, int pos1, int pos2);
    // // evaluate the cost of cheapest insertion of node in this route
    // // return a <position, cost> pair
    // // If no insertion is feasible, the cost is INT_MAX
    pair<int, double> cheapestInsertion(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist, int node, vector<int> &positions);
    
    //testing insertion of parcel pickup and delivery simultaneously
    void cheapestInsertionParcel(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist, int node, int node2, vector<int> &positions, vector<int> &positions2, vector< pair<int, double> > &bestMove, probStat* problem);

    // // pre-condition: the insertion is feasible
    void insert(instanceStat *inst, double **Mdist, int node, int position, double profit);

    void erase(instanceStat *inst, double **Mdist, int position, double profit);

    void insertBlock(instanceStat *inst, double **Mdist, 
                    vector<int> blockNodes, int position, double profit);
    
    void eraseBlock(instanceStat *inst, double **Mdist, 
                        int pos1, int pos2, double profit);

    void printLoad();
    void printPDVec();
    void printPP();

    pair <int, int> getInterval(int req);
    
    bool checkInterval(instanceStat *inst, int pos1, int pos2, pair <int, int> inter1, pair <int, int> inter2);
    bool checkDelivery(instanceStat *inst, int pos1, int pos2, probStat* problem);

    double Swap(instanceStat *inst, double **Mdist, vector<nodeStat> &nodeVec, probStat* problem);
    
    double relocateK(instanceStat *inst, double **Mdist, vector<nodeStat> &nodeVec, probStat* problem, int k);

    double TwoOpt(instanceStat *inst, double **Mdist, vector<nodeStat> &nodeVec, probStat* problem);

    // double ThreeOpt(instanceStat *inst, double **Mdist, vector<nodeStat> &nodeVec, probStat* problem);
    //updates passenger, loads, times, clears PD and updates this vector.
    void updateAll (instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist);

    double rmvVal (instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist, int candidate, bool isparcel);


};

#endif
