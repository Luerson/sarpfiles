#include "sarpSolution.h"


void sarpSolution::addRoute(sarpRoute *route){
    this->routes.push_back(*route);
    this->updateCost();
}

void sarpSolution::updateCost() {
    this->cost = 0;
    for (int i = 0; i < routes.size(); ++i) {
        this->cost += routes[i].cost();
    }
}

void sarpSolution::updateRoutes(sarpRoute *route, int idr) {
    this->routes[idr] = *route;
}

void sarpSolution::updateVehicles(){
    this->usedK++;
}

void sarpSolution::printSol(instanceStat *inst) {

    sarpRoute sroute(inst, usedK);

    cout << "SARP Solution: " << endl;
    for (int i = 0; i < routes.size(); i++){
        sroute = routes[i];
        for(auto j: sroute){
            cout << j << " - ";
        }
        sroute.printTotalTime();
        cout << endl;
    }
    cout << endl;
}
void sarpSolution::printLoadSol(instanceStat *inst) {

    sarpRoute sroute(inst, usedK);

    cout << "Solution Load: " << endl;
    for (int i = 0; i < routes.size(); i++){
        sroute = routes[i];
        cout << endl << "**********Load************" << endl;
        cout << "Route " << i << ": " << endl;
        sroute.printLoad();
        cout << endl << "**********end Load********" << endl;
    }
    cout << endl;
}
void sarpSolution::addtounserved(int candidate) {
    unserved.push_back(candidate);
}   


void sarpSolution::printCosts(){
    cout << "Cost of routes: " << endl;

    for (int i = 0; i < this->routes.size(); i++){
        cout << "Route " << i << ": " << routes[i].cost() << endl;
    }

    cout << "\nCost of routes: " << this->cost << endl;
}

double sarpSolution::relocate (instanceStat *inst, vector<nodeStat> &nodeVec, 
                                double **Mdist, probStat* problem,
                                 int rid1, int rid2, int &currCand, pair <int, int> &currPairPos){

    sarpRoute r1(inst, 0);
    sarpRoute r2(inst, 0);
    int rsize1, rsize2;

    r1 = getRoute(rid1);
    r2 = getRoute(rid2);

    rsize1 = r1.getNodesSize();

    // cout << "Size of route 1: " << rsize1 << endl;
    // getchar();
    rsize2 = r2.getNodesSize(); 

    bool feasible = 0;

    int candidate, candidate2;
    int bestpos1, bestpos2;
    int inCurrent;

    vector<int> inspositions; 
    vector<int> inspositions2;

    pair <int, double> cheapestpair; //position and delta
    
    cheapestpair.first = -1;
    cheapestpair.second = -100000;
    
    vector< pair <int, double> > cheapestMove;
    
    cheapestMove.push_back(cheapestpair);
    cheapestMove.push_back(cheapestpair); //chapest move has 2 cheapest pairs. Only 1 is used for passengers.
    
    double best_cost, compareCost, iniCost, rval1, rval2, rmvVal, addVal;

    rval1 = r1.cost();
    rval2 = r2.cost();

    best_cost = rval1 + rval2;
    iniCost = rval1 + rval2;

    currCand = -1;

    cout << "Routes: " << rid1 << " and " << rid2;

    for (int i = 1; i < rsize1-1; i++){
        candidate = r1.getReq(i);
        compareCost = rval1 + rval2;
        cout << "\ncandidate: " << candidate << endl;

        cout << "Initial cost of both routes: " << compareCost << endl;

        if(nodeVec[candidate].load < 0){
            // cout << "\n Delivery node (skip)";
            // getchar();
            continue;
        }

        inspositions.clear();
        r2.availablePos(inst, nodeVec, candidate, problem, inspositions);

        if (inspositions.size() < 1){
            continue;
        }
        
        // cout << "Available positions for insertion: " << endl;
        // for (int p = 0; p < inspositions.size(); p++){
        //     cout << inspositions[p] << " - ";
        // }
        // cout << endl;
        // getchar();

        if (candidate < inst->n){
            cheapestpair = r2.cheapestInsertion(inst, nodeVec, Mdist, candidate, inspositions);
            // cout << "end of cheapest insertion." << endl;
            // getchar();
            addVal = cheapestpair.second;

            cout << "Passenger add val: " << addVal << endl;

            if (addVal > 0){
                rmvVal = r1.rmvVal(inst, nodeVec, Mdist, i, 0);
                cout << "Calculate rmvval: " << rmvVal << endl;
                compareCost += addVal + rmvVal;           
            }
            else{
                continue;
            }                
            
            cout << "Passenger added at position: " << cheapestpair.first << endl;

            // cout << "Compare Cost: " << compareCost;
        }
        else {
            cheapestMove[0].first = -1;
            cheapestMove[0].second = -100000;
            cheapestMove[1].first = -1;
            cheapestMove[1].second = -100000;
            // cout << "cheapest move 2 BEFORE: " <<  cheapestMove[1].second << endl;
            int node2 = candidate + inst->m;
            inspositions2.clear();
            r2.cheapestInsertionParcel(inst, nodeVec, Mdist, candidate, node2, inspositions, inspositions2, cheapestMove, problem);
            
            // if (inspositions2.size() < 1){
            //     continue;
            // }

            if(cheapestMove[0].first != cheapestMove[1].first){
                addVal = cheapestMove[0].second + cheapestMove[1].second;
            }
            else{
                addVal = nodeVec[candidate].profit + cheapestMove[1].second;
                // cout << "cheapest move 2: " <<  cheapestMove[1].second << endl;
            }  
 
            if (addVal > 0){
                cout << "New add value: " << addVal << endl;
                // getchar();

                rmvVal = r1.rmvVal(inst, nodeVec, Mdist, i, 1);

                cout << "New rmv value: " << rmvVal << endl;
                // getchar();
                
                cout << "Adding at position: " << cheapestMove[0].first << " and " << cheapestMove[1].first << endl;
                // getchar();

                compareCost += addVal + rmvVal;  
            }
            else{
                continue;
            }
            // cout << "Compare Cost: " << compareCost;
        }

        if (compareCost > best_cost){
            // cout << "\nThere was an improvement" << endl;
            // cout << "\nCandidate: " << candidate << endl;
            // getchar();
            // inCurrent = i;
            currCand = i;
            // currCand = &inCurrent;
            cout << "Current candidate in solutions function: " << currCand << endl;
            // getchar();
            if (candidate < inst->n){
                best_cost = compareCost;
                currPairPos.first = cheapestpair.first;
                currPairPos.second = -1;

                // cout << "insert in position: "<< currPairPos.first << endl; 
                // cout << "with an improvement of " << best_cost << endl;
            }
            else{
                best_cost = compareCost;
                currPairPos.first = cheapestMove[0].first;
                currPairPos.second = cheapestMove[1].first;

                // cout << "insert in position: "<< currPairPos.first << endl; 
                // cout << "with an improvement of " << best_cost << endl;
            }

        }
        // else{
        //     currCand = -1;
        // }
    }

    double delta;

    delta = iniCost - best_cost;

    cout << "Best cost: " << best_cost << endl;
    cout << "cost (previous cost): " << iniCost << endl;
    cout << "delta: " << delta << endl;
    // getchar();
    // currCand = &inCurrent;

    cout << "Current candidate passing out of function: " << currCand << endl;
    // getchar();

    return delta;
}

//Needs to be tested with block of only parcels
pair <double, double> sarpSolution::calcBlockTimes(instanceStat *inst, 
                  vector<nodeStat> &nodeVec, 
                  double **Mdist, vector<int> newBlock){

    pair <double, double> timeForBlock;
    timeForBlock.first = 0;
    
    timeForBlock.second = 0;

    int lastPass, lastPassPos;

    lastPassPos = -1;
    lastPass = -1;

    bool newpass;

    //accumulating travel time for total travel time involved.
    for (int i = 0; i < newBlock.size() - 1; i++){
        int req = newBlock[i];
        int nextreq = newBlock[i + 1];
        newpass = 0;

        timeForBlock.second += nodeVec[req].delta +
                   ((Mdist[req][nextreq])/inst->vmed);

        //setting up for calculating block endtime
        if (req < inst->n){
            newpass = 1;
            lastPass = req;
            lastPassPos = i;
            timeForBlock.first = 0;
        }
        if (lastPass > -1){
            timeForBlock.first += nodeVec[req].delta +
                    ((Mdist[req][nextreq])/inst->vmed);
        }
    }

    timeForBlock.second += nodeVec[newBlock.back()].delta;

    if (lastPass > -1){
        timeForBlock.first += nodeVec[lastPass].e;
        timeForBlock.first += nodeVec[newBlock.back()].delta;
    }
    else{
        timeForBlock.first = timeForBlock.second;
    }

    return timeForBlock;
}


pair <double, double> sarpSolution::TwoOpt(instanceStat *inst, double **Mdist,
                    vector<nodeStat> &nodeVec, int rid1, int rid2, 
                    pair <int, int> &currPairPos, probStat* problem){

    pair<double, double> delta, tBlock1, tBlock2;

    delta.first = 0;
    delta.second = 0;

	pair<double, double> bestDelta;

	double newCost = 0;

    sarpRoute sroute1(inst, 0);
    sarpRoute sroute2(inst, 0);
	
    sroute1 = getRoute(rid1);
    sroute2 = getRoute(rid2);

    bestDelta.first = sroute1.cost();
    bestDelta.second = sroute2.cost();

    int pos1 = 0;
	int pos2 = 0;

    int jstart, jend;

    pair <int, int> interpos1, interpos2;
    
    vector<int> tempVector, zeroPos1, zeroPos2, zeroVec1, zeroVec2, bestBlock1, bestBlock2;
    
    int temp;

    bool feasible, improve;

    improve = 0;
    // cout << "Current cost: " << cost_ << endl;
    // getchar();
    int loadSize1, loadSize2;
    double newEndTime1, newEndTime2, bestp1, bestp2;
    
    loadSize1 = sroute1.getLoadSize();
    loadSize2 = sroute2.getLoadSize();

    //-2 because two opt with a single request is already included into 'relocate inter'
    for (int i = 1; i < loadSize1 - 2; i++){
        int currLoad, request;
        currLoad = sroute1.getReqLoad(i);
        request = sroute1.getReq(i);
        if (currLoad == 0){
            if (nodeVec[request].load == 0){
                zeroPos1.push_back(i);
            }
        }
        else if (currLoad == 1){
            if (nodeVec[request].load > 0 ){
                zeroPos1.push_back(i);
            }
        }
    }

    for (int i = 1; i < loadSize2 - 2; i++){
        int currLoad, request;
        currLoad = sroute2.getReqLoad(i);
        request = sroute2.getReq(i);
        if (currLoad == 0){
            if (nodeVec[request].load == 0){
                zeroPos2.push_back(i);
            }
        }
        else if (currLoad == 1){
            if (nodeVec[request].load > 0 ){
                zeroPos2.push_back(i);
            }
        }
    }

    cout << "Zero Pos vec (1): " << endl;

    for (int j = 0; j < zeroPos1.size(); j++){
        cout << zeroPos1[j] << " ";
    }
    cout << endl;
    // getchar();

    cout << "Zero Pos vec (2): " << endl;

    for (int j = 0; j < zeroPos2.size(); j++){
        cout << zeroPos2[j] << " ";
    }
    cout << endl;
    // getchar();

    for (int i = 0; i < zeroPos1.size(); i++){
        pos1 = zeroPos1[i];
        zeroVec1.clear();
        //making zero sum vec (new block) from r1
        for (int z = pos1; z < loadSize1 - 1; z++){
            int req = sroute1.getReq(z);
            zeroVec1.push_back(req);
        }
        
        tBlock1 = calcBlockTimes(inst, nodeVec, Mdist, zeroVec1);

        int iniPos1 = pos1;
        int endPos1 = loadSize1 - 1;
        cout << "r1: inipos: " << iniPos1 << " - endPos: " << endPos1 << endl;
        getchar();

        if (iniPos1 > 1){
            newEndTime1 = sroute1.blockrmvTime(inst, nodeVec, Mdist, iniPos1, endPos1);
        }
        else{
            newEndTime1 = 0;
        }

        cout << "New End time: " << newEndTime1 << endl;
        getchar();

        for (int j = 0; j < zeroPos2.size(); j++){
            pos2 = zeroPos2[j];
            zeroVec2.clear();
            //making zero sum vec (new block) from r2
            for (int z = pos2; z < loadSize2 - 1; z++){
                int req = sroute2.getReq(z);
                zeroVec2.push_back(req);
            }

            tBlock2 = calcBlockTimes(inst, nodeVec, Mdist, zeroVec2);
            int iniPos2 = pos2;
            int endPos2 = loadSize2 - 1;
            cout << "r2: inipos: " << iniPos2 << " - endPos: " << endPos2 << endl;
            getchar();

            if (iniPos2 > 1){
                newEndTime2 = sroute1.blockrmvTime(inst, nodeVec, Mdist, iniPos2, endPos2);
            }
            else{
                newEndTime2 = 0;
            }

            cout << "Times for block 2: " << tBlock2.first << " - " << tBlock2.second << endl;
            getchar();

            newEndTime2 = 0;

            //testing piece of route 2 into route 1

            cout << "Testing feasibility r1" << endl;
            feasible = sroute1.testBlockIns(inst, nodeVec, 
                                            Mdist, newEndTime1, tBlock2,
                                            iniPos1, endPos1, zeroVec2);

            cout << "Feasible insertion r1: " << feasible << endl;
            getchar();


            if (feasible){
                //testing piece of route 1 into route 2
                feasible = 0;
                cout << "Testing feasibility r2" << endl;
                feasible = sroute2.testBlockIns(inst, nodeVec, 
                                                Mdist, newEndTime2, tBlock1,
                                                iniPos2, endPos2, zeroVec1);

                cout << "Feasible insertion r2: " << feasible << endl;
                getchar();
            }
            else{
                continue;
            }

            if (feasible){
                double p1, p2;
                p1 = sroute1.blockProfit(inst, nodeVec, Mdist, iniPos1, endPos1);
                p2 = sroute2.blockProfit(inst, nodeVec, Mdist, iniPos2, endPos2);

                interpos1.first = sroute1.getReq(iniPos1 - 1);
                interpos1.second = sroute1.getReq(endPos1 + 1);

                interpos2.first = sroute2.getReq(iniPos2 - 1);
                interpos2.second = sroute2.getReq(endPos2 + 1);

                //delta is the whole new cost, not only the difference
                delta.first = sroute1.cost() - p1 + p2 - (inst->costkm*(Mdist[interpos1.first][zeroVec2[0]] + Mdist[zeroVec2.back()][interpos1.second]));

                delta.second = sroute2.cost() - p2 + p1 - (inst->costkm*(Mdist[interpos2.first][zeroVec1[0]] + Mdist[zeroVec1.back()][interpos2.second]));

                if ((delta.first + delta.second) > (bestDelta.first + bestDelta.second)){
                    bestDelta.first = delta.first;
                    bestDelta.second =  delta.second;

                    currPairPos.first = iniPos1;
                    currPairPos.second = iniPos2;

                    bestBlock1.clear();
                    bestBlock2.clear();

                    bestBlock1 = zeroVec1;
                    bestBlock2 = zeroVec2;

                    bestp1 = p1;
                    bestp2 = p2;

                    improve = 1;
                }
            }


        }
    }

    if (improve){
        sroute1.eraseBlock(inst, Mdist, currPairPos.first, loadSize1-1, bestp1);
        sroute2.eraseBlock(inst, Mdist, currPairPos.second, loadSize2-1, bestp2);

        sroute1.insertBlock(inst, Mdist, zeroVec2, currPairPos.first, bestp2);
        sroute2.insertBlock(inst, Mdist, zeroVec1, currPairPos.second, bestp1);
    }



    return bestDelta;
}

void sarpSolution::addunserved(instanceStat *inst, vector<nodeStat> &nodeVec, double **Mdist, probStat* problem){
    sarpRoute sroute(inst, 0);
    int candidate, candidate2;

    vector<int> inspositions, inspositions2;
    pair <int, double> cheapestpair;

    cheapestpair.first = -1;
    cheapestpair.second = -100000;

    vector< pair<int, double> > cheapestMove;
    pair <int, int> bestpositions;

    double best_cost;
    int best_route;
    bool inserted;
    bool feastime = 0;
    int vehicle;

    this->stats.setStart();

    cheapestMove.push_back(cheapestpair);
    cheapestMove.push_back(cheapestpair);

    for (int cand = 0; cand < unserved.size(); cand++){
        candidate = unserved[cand];
        candidate2 = candidate + inst->m;

        best_route = -1;
        best_cost = -100000;
        inserted = 0;
        feastime = 0;
        bestpositions.first = -1;
        bestpositions.second = -1;

        for (int rid = 0; rid < routes.size(); rid++){
            cheapestMove[0].first = -1;
            cheapestMove[0].second = -100000;
            cheapestMove[1].first = -1;
            cheapestMove[1].second = -100000;

            sroute = getRoute(rid);

            inspositions.clear();
            inspositions2.clear();
            
            sroute.availablePos(inst, nodeVec, candidate, problem, inspositions);

            sroute.cheapestInsertionParcel(inst, nodeVec, Mdist, candidate, candidate2, inspositions, inspositions2, cheapestMove, problem);
            
            if (cheapestMove[0].first == cheapestMove[1].first){
                if (cheapestMove[0].second > best_cost){
                    best_cost = cheapestMove[0].second;
                    best_route = rid;
                    bestpositions.first = cheapestMove[0].first;
                    bestpositions.second = cheapestMove[1].first;
                    inserted = 1;
                }                    
            }
            else{
                if (cheapestMove[0].second + cheapestMove[1].second > best_cost){
                    best_cost = cheapestMove[0].second + cheapestMove[1].second;
                    best_route = rid;
                    bestpositions.first = cheapestMove[0].first;
                    bestpositions.second = cheapestMove[1].first;
                    inserted = 1;
                }
            }
        }
        if (best_cost == -100000){
            if (problem->dParcel > 0){
                vehicle = getvehicle();
                sarpRoute newroute(inst, vehicle);
                feastime = newroute.fInsertionParcel(inst, nodeVec, Mdist, candidate, candidate2);
                newroute.calcCost(inst, nodeVec, Mdist);
                newroute.insert(inst, Mdist, candidate, 1, nodeVec[candidate].profit);
                newroute.insert(inst, Mdist, candidate2, 2, nodeVec[candidate2].profit);
                newroute.updateAll(inst, nodeVec, Mdist);
                addRoute(&newroute);

                for (int u = 0; u < unserved.size(); ++u) {
                    if (this->unserved[u] == candidate) {
                        this->unserved.erase(this->unserved.begin() + u);
                        break;
                    }
                }

                updateVehicles();
                inserted = 1;
                updateCost();
                break;
            }
            else{
                addtounserved(candidate);
                for (int u = 0; u < unserved.size(); ++u) {
                    if (this->unserved[u] == candidate) {
                        this->unserved.erase(this->unserved.begin() + u);
                        break;
                    }
                }
                break;
            }
        }

        else{
            sroute = getRoute(best_route);

            // cout << "::::******::::::******::::::" << endl;
            // cout << "Request1: " << candidate << " - Request2: " << candidate2 << endl;
            // cout << "Best position 1: " << bestpositions.first << "; Best position 2: " << bestpositions.second << endl;
            // cout << "::::******::::::******::::::" << endl;

            // getchar();

            sroute.insert(inst, Mdist, candidate, bestpositions.first, nodeVec[candidate].profit);
            sroute.insert(inst, Mdist, candidate2, bestpositions.second+1, nodeVec[candidate2].profit);
            sroute.updateAll(inst, nodeVec, Mdist);
            updateRoutes(&sroute, best_route);
            updateCost();

            inserted = 1;
            // cout << "*****Route with insertion: " << endl;
            // for (auto a: sroute){
            //     cout << a << " - ";
            // }
            // cout << "New route cost: " << sroute.cost() << endl;
            // cout << endl << endl;

            for (int u = 0; u < unserved.size(); ++u) {
                if (this->unserved[u] == candidate) {
                    this->unserved.erase(this->unserved.begin() + u);
                    break;
                }
            }
            break;
        }

    }

    cout << "After add unserved" << endl;
    printSol(inst);
    printCosts();

    this->stats.setEnd();
    cout << "\nAdd Unserved Time: " << std::setprecision(8) << this->stats.printTime() << endl;

}
