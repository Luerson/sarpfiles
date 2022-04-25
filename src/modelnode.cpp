#include "modelnode.h"
#include <cstdlib>
#include <stdio.h>

void initArcs (instanceStat *inst, nodeArcsStruct *nas){
    vector<bool> auxVec;
    vector< pair<int,int> > auxPairVec;

    vector<int> aux1d;
    vector< vector<int> > aux2d;

    nas->arcPlus.clear();
    nas->vArcPlus.clear();
    nas->vArcMinus.clear();
    nas->arcV.clear();
    nas->arcPlus.clear();
    nas->arcMinus.clear();
    auxVec.clear();


    for (int k = 0; k < inst->K; k++){
        nas->arcPlus.push_back(auxPairVec);
    }

    for(int i = 0; i < inst->V + inst->dummy; i++){
        aux2d.push_back(aux1d);

        nas->vArcPlus.push_back(nas->arcPlus);
        nas->vArcMinus.push_back(nas->arcPlus);
    }

    nas->arcPlus.clear();

    for(int i = 0; i < inst->V + inst->dummy; i++){
        for(int j = 0; j < inst->V + inst->dummy; j++){
            auxVec.push_back(false);
        }
        nas->arcs.push_back(auxVec);
        nas->arcPlus.push_back(auxPairVec);
        nas->arcMinus.push_back(auxPairVec);
        


        nas->arcV.push_back(aux2d);
        
        auxVec.clear();
    }

} 

void feasibleArcs (instanceStat *inst, nodeArcsStruct *nas, probStat* problem, vector<nodeStat> &nodeVec, double **mdist){
    int auxK;

    int fDepot = inst->n + 2*inst->m;
    int fDummy = inst->n + 2*inst->m + inst->K;

    //independently from sarp scenario, these arcs are always true

    for(int i = inst->n + 2*inst->m; i < inst->V; i++){//i is a starting point
        // for (int j = 0; j < inst->n + inst->m; j++){//j is a passenger or parcel pickup node
        //     nas->arcs[i][j] = true;
        //     nas->fArc.first = i;
        //     nas->fArc.second = j;
        //     nas->arcMinus[j].push_back(nas->fArc);
        //     nas->arcPlus[i].push_back(nas->fArc);
        //     nas->allArcs.push_back(nas->fArc);

        //     nas->arcnf.push_back(nas->fArc);

        //     auxK = i - fDepot;
        //     nas->arcV[i][j].push_back(auxK);
        // }

        int j = i + inst->K;

        nas->arcs[i][j] = true;
        nas->fArc.first = i;
        nas->fArc.second = j;
        nas->arcMinus[j].push_back(nas->fArc);
        nas->arcPlus[i].push_back(nas->fArc);
        nas->allArcs.push_back(nas->fArc);

        auxK = j - inst->V;
        nas->arcV[i][j].push_back(auxK);    
    }

    // for (int i = 0; i < inst->n; i++){//i is a passenger node
    //     for(int j = 0; j < inst->n; j++){// j is a passenger req
    //         if(i != j){
    //             double ttij = mdist[i][j]/inst->vmed;//travel time between requests i and j 
    //             //if lowest time for req i + travel time from i to j is lower or equal to
    //             //the latest point in time to serve request j. If latest time == T, it is always valid                        
            
    //             if (nodeVec[i].e + ttij <= nodeVec[j].l){
    //                 nas->arcs[i][j] = true;
    //                 nas->fArc.first = i;
    //                 nas->fArc.second = j;
    //                 nas->arcMinus[j].push_back(nas->fArc);
    //                 nas->arcPlus[i].push_back(nas->fArc);

    //                 nas->arcNN.push_back(nas->fArc);
    //                 nas->arcNplus.push_back(nas->fArc);

    //                 nas->allArcs.push_back(nas->fArc);
    //                 nas->arcnf.push_back(nas->fArc);
    //                 for (int k = 0; k < inst->K; k++){
    //                     nas->arcV[i][j].push_back(k);
    //                 }                        
    //             }
    //         }
    //     }

    //     for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
    //         nas->arcs[i][j] = true;
    //         nas->fArc.first = i;
    //         nas->fArc.second = j;
    //         nas->arcMinus[j].push_back(nas->fArc);
    //         nas->arcPlus[i].push_back(nas->fArc);

    //         nas->arcNplus.push_back(nas->fArc);

    //         nas->allArcs.push_back(nas->fArc);
    //         auxK = j - inst->V;
    //         nas->arcV[i][j].push_back(auxK);
    //     }
    // }

    // for (int i = inst->n + inst->m; i < inst->n + 2*inst->m; i++){//i is a parcel dl node           
    //     for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
    //         nas->arcs[i][j] = true;
    //         nas->fArc.first = i;
    //         nas->fArc.second = j;
    //         nas->arcMinus[j].push_back(nas->fArc);
    //         nas->arcPlus[i].push_back(nas->fArc);
    //         nas->allArcs.push_back(nas->fArc);

    //         auxK = j - inst->V;
    //         nas->arcV[i][j].push_back(auxK);
    //     }
    // }
 

    // if (problem->p1 < 0){ // arcs for baseline scenarios
    //     //passenger-passenger
    //     if (problem->scen == "BL2"){
    //         for (int i = 0; i < inst->n; i++){ //i is a passenger request
    //             for (int j = inst->n; j < inst->n + inst->m; j++){// j is a parcel pu req
    //                 double ttij = mdist[i][j]/inst->vmed;//travel time between requests i and j 
    //                 //if lowest time for req i + travel time from i to j is lower or equal to
    //                 //the latest point in time to serve request j. If latest time == T, it is always valid                        
    //                 if (nodeVec[i].e + ttij <= nodeVec[j].l){                
    //                     nas->arcs[i][j] = true;
    //                     nas->fArc.first = i;
    //                     nas->fArc.second = j;
    //                     nas->arcMinus[j].push_back(nas->fArc);
    //                     nas->arcPlus[i].push_back(nas->fArc);

    //                     nas->arcNplus.push_back(nas->fArc);

    //                     nas->allArcs.push_back(nas->fArc);
    //                     nas->arcnf.push_back(nas->fArc);
    //                     for (int k = 0; k < inst->K; k++){
    //                         nas->arcV[i][j].push_back(k);
    //                     }       
    //                 }                 
    //             }                
    //         }
    //     }
    //     //parcelpu-parcel
    //     for (int i = inst->n; i < inst->n + inst->m; i++){//i is a parcel pu node
    //         for (int j = inst->n; j < inst->n + 2*inst->m; j++){//j is a parcel pu or dl node
    //             if (i != j){
    //                 nas->arcs[i][j] = true;
    //                 nas->fArc.first = i;
    //                 nas->fArc.second = j;
    //                 nas->arcMinus[j].push_back(nas->fArc);
    //                 nas->arcPlus[i].push_back(nas->fArc);
    
    //                 nas->allArcs.push_back(nas->fArc);
    //                 nas->arcnf.push_back(nas->fArc);
    //                 for (int k = 0; k < inst->K; k++){
    //                     nas->arcV[i][j].push_back(k);
    //                 } 
    //             }
    //         }
    //     }
    //     //parceldl-parcel
    //     for (int i = inst->n + inst->m; i < inst->n + 2*inst->m; i++){//i is a parcel dl node
    //         for (int j = inst->n; j < inst->n + 2*inst->m; j++){//j is a parcel pu or dl node
    //             if (i != j && j != i - inst->m){//different nodes; no dl to its pickup
    //                 nas->arcs[i][j] = true;
    //                 nas->fArc.first = i;
    //                 nas->fArc.second = j;
    //                 nas->arcMinus[j].push_back(nas->fArc);
    //                 nas->arcPlus[i].push_back(nas->fArc);
    //                 nas->allArcs.push_back(nas->fArc);
    //                 nas->arcnf.push_back(nas->fArc);
    //                 for (int k = 0; k < inst->K; k++){
    //                     nas->arcV[i][j].push_back(k);
    //                 } 
    //             }
    //         }
    //         if (problem->scen == "BL2"){
    //             for (int j = 0; j < inst->n; j++){//j is a passenger node
    //                 nas->arcs[i][j] = true;
    //                 nas->fArc.first = i;
    //                 nas->fArc.second = j;
    //                 nas->arcMinus[j].push_back(nas->fArc);
    //                 nas->arcPlus[i].push_back(nas->fArc);

    //                 nas->arcPN.push_back(nas->fArc);

    //                 nas->allArcs.push_back(nas->fArc);
    //                 nas->arcnf.push_back(nas->fArc);
    //                 for (int k = 0; k < inst->K; k++){
    //                     nas->arcV[i][j].push_back(k);
    //                 }  
    //             } 
    //         } 
    //     }
    // }

    // else{
    //     for (int i = 0; i < inst->n; i++){//i is a passenger node
    //         for(int j = inst->n; j < inst->n + 2*inst->m; j++){// j is a parcel req (pu or del)
    //             double ttij = mdist[i][j]/inst->vmed;//travel time between requests i and j 
    //             //if lowest time for req i + travel time from i to j is lower or equal to
    //             //the latest point in time to serve request j. If latest time == T, it is always valid                        
    //             if (nodeVec[i].e + ttij <= nodeVec[j].l){
    //                 nas->arcs[i][j] = true;
    //                 nas->fArc.first = i;
    //                 nas->fArc.second = j;
    //                 nas->arcMinus[j].push_back(nas->fArc);
    //                 nas->arcPlus[i].push_back(nas->fArc);

    //                 nas->arcNplus.push_back(nas->fArc);

    //                 nas->allArcs.push_back(nas->fArc);
    //                 nas->arcnf.push_back(nas->fArc);
    //                 for (int k = 0; k < inst->K; k++){
    //                     nas->arcV[i][j].push_back(k);
    //                 }
    //             }
    //         }
    //     }

    //     for (int i = inst->n; i < inst->n + inst->m; i++){//i is a parcel pu node           
    //         for (int j = 0; j < inst->n; j++){ //j is a passenger node
    //             nas->arcs[i][j] = true;
    //             nas->fArc.first = i;
    //             nas->fArc.second = j;
    //             nas->arcMinus[j].push_back(nas->fArc);
    //             nas->arcPlus[i].push_back(nas->fArc);
                
    //             nas->arcPN.push_back(nas->fArc);

    //             nas->allArcs.push_back(nas->fArc);
    //             nas->arcnf.push_back(nas->fArc);
    //             for (int k = 0; k < inst->K; k++){
    //                 nas->arcV[i][j].push_back(k);
    //             }
    //         }
    //     }

    //     for (int i = inst->n + inst->m; i < inst->n + 2*inst->m; i++){//i is a parcel dl node           
    //         for (int j = 0; j < inst->n + inst->m; j++){ //j is a passenger node or a parcel pu node
    //             if (j + inst->m != i){
    //                 nas->arcs[i][j] = true;
    //                 nas->fArc.first = i;
    //                 nas->fArc.second = j;
    //                 nas->arcMinus[j].push_back(nas->fArc);
    //                 nas->arcPlus[i].push_back(nas->fArc);

    //                 if (j < inst->n){//j is a passenger node
    //                     nas->arcPN.push_back(nas->fArc);
    //                 }

    //                 nas->allArcs.push_back(nas->fArc);
    //                 nas->arcnf.push_back(nas->fArc);
    //                 for (int k = 0; k < inst->K; k++){
    //                     nas->arcV[i][j].push_back(k);
    //                 }
    //             }
    //         }
    //         for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
    //             nas->arcs[i][j] = true;
    //             nas->fArc.first = i;
    //             nas->fArc.second = j;
    //             nas->arcMinus[j].push_back(nas->fArc);
    //             nas->arcPlus[i].push_back(nas->fArc);

    //             nas->allArcs.push_back(nas->fArc);
    //             auxK = j - inst->V;
    //             nas->arcV[i][j].push_back(auxK);
    //         }
    //     }

    //     //%%%%%%%%%%%%%%%%%%%%%%%
    //     //specific cases
    //     //%%%%%%%%%%%%%%%%%%%%%%%

    //     if (problem->p2 > 0){ //multi parcel
    //         for (int i = inst->n; i < inst->n + 2*inst->m; i++){//i is a parcel pu or dl node                   
    //             for (int j = inst->n; j < inst->n + 2*inst->m; j++){//j is a parcel pu or dl node
    //                 if (j + inst->m != i && i != j && i + inst->m != j){ // no dl to its pu; no pu to its dl; not a node to itself
    //                     nas->arcs[i][j] = true;
    //                     nas->fArc.first = i;
    //                     nas->fArc.second = j;
    //                     nas->arcMinus[j].push_back(nas->fArc);
    //                     nas->arcPlus[i].push_back(nas->fArc);

    //                     nas->allArcs.push_back(nas->fArc);
    //                     nas->arcnf.push_back(nas->fArc);
    //                     for (int k = 0; k < inst->K; k++){
    //                         nas->arcV[i][j].push_back(k);
    //                     }
    //                 }                    
    //             }
    //         }
    //     }

    //     if (problem->dParcel > 0){//direct parcel delivery
    //         for (int i = inst->n; i < inst->n + inst->m; i++){//i is a parcel pu node
    //             int j = i + inst->m; //j is i's delivery location

    //             nas->arcs[i][j] = true;
    //             nas->fArc.first = i;
    //             nas->fArc.second = j;
    //             nas->arcMinus[j].push_back(nas->fArc);
    //             nas->arcPlus[i].push_back(nas->fArc);
    //             nas->allArcs.push_back(nas->fArc);
    //             nas->arcnf.push_back(nas->fArc);
    //             for (int k = 0; k < inst->K; k++){
    //                 nas->arcV[i][j].push_back(k);
    //             }                                                
    //         }
    //     }
    // }

    // for (int a = 0; a < nas->allArcs.size(); a++){
    //     int i = nas->allArcs[a].first;
    //     int j = nas->allArcs[a].second;

    //     for(int k1 = 0; k1 < nas->arcV[i][j].size(); k1++){
    //         int k = nas->arcV[i][j][k1];
    //         nas->vArcPlus[i][k].push_back(nas->allArcs[a]);
    //         nas->vArcMinus[j][k].push_back(nas->allArcs[a]);
    //     }

    // }

    // //building reqV
    // //for every request, for every vehicle, check if the vehicle start time + time of
    // // travel from the starting point to request is < starting time of request TW
    // //only parcel pickup and customer
    // vector <int> vehs;

    // for (int i = 0; i < inst->n + inst->m; i++){
    //     for (int k = 0; k < inst->K; k++){
    //         int startdepot = inst->n + 2*inst->m + k;
    //         double ttSI = mdist[startdepot][i]/inst->vmed;//travel time between starting depot of k and request i
    //         //if the online time of vehicle k + travel time to node i is lower or equal to
    //         //the latest point in time to serve request i. If latest time == T, it is always valid
    //         if (nodeVec[startdepot].e + ttSI <= nodeVec[i].l){
    //             vehs.push_back(k);
    //         }
    //     }
    //     nas->reqV.push_back(vehs);
    //     vehs.clear();
    // }

    // // cout << "allowed vehicles: " << endl;

    // // for (int i = 0; i < nas->arcV.size(); i++){
    // //     for (int j = 0; j < nas->arcV[i].size(); j++){
    // //         if (nas->arcs[i][j]){
    // //             printf("\narc %i, %i, Vehicles: ", i, j);
    // //             for (int k = 0; k < nas->arcV[i][j].size(); k++){
    // //                 cout << nas->arcV[i][j][k] << endl;
    // //             }
    // //             cout << endl;                
    // //         }
    // //     }
    // // }
    // // getchar();

    // int auxK;

    // int fDepot = inst->n + 2*inst->m;
    // int fDummy = inst->n + 2*inst->m + inst->K;

    if (problem->scen == "1A" || problem->scen == "2A"){ 
        for (int i = 0; i < inst->V; i++){
            if (i < inst->n){//i is a passenger req
                for(int j = 0; j < inst->n + 2*inst->m; j++){// j is a parcel req (pu or del)
                    if(i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        if (j < inst->n){
                            nas->arcNN.push_back(nas->fArc);
                        }
                        nas->arcNplus.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }
                    }
                }

                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->arcNplus.push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->n + inst->m){//i is a parcel pickup node
                for (int j = 0; j < inst->n; j++){ //j is a passenger node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->arcPN.push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    nas->arcnf.push_back(nas->fArc);
                    for (int k = 0; k < inst->K; k++){
                        nas->arcV[i][j].push_back(k);
                    }
                }
            }

            else if (i < inst->n + 2*inst->m){// i is a parcel delivery node
                for (int j = 0; j < inst->n + inst->m; j++){//j is a passenger node or parcel pickup node
                    if (j + inst->m != i){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }
                    }                    
                }

                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->V + inst->dummy){ // i is a starting node
                for (int j = 0; j < inst->n + inst->m; j++){//j is a passenger or parcel pickup node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    nas->arcnf.push_back(nas->fArc);
                    auxK = i - fDepot;
                    nas->arcV[i][j].push_back(auxK);
                }
            }          
        }
    }
    else if (problem->scen == "1B" || problem->scen == "2B"){
        for (int i = 0; i < inst->V; i++){
            if (i < inst->n){//i is a passenger req
                for(int j = 0; j < inst->n + 2*inst->m; j++){// j is a parcel req (pu or del)
                    if(i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        if (j < inst->n){
                            nas->arcNN.push_back(nas->fArc);
                        }
                        nas->arcNplus.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }                        
                    }
                }
                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->arcNplus.push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->n + inst->m){//i is a parcel pickup node
                for (int j = 0; j < inst->n + inst->m; j++){ //j is a passenger or parcel pickup node
                    if (i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        nas->arcPN.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }                          
                    }
                }
            }
            else if (i < inst->n + 2*inst->m){// i is a parcel delivery node
                for (int j = 0; j < inst->n + 2*inst->m; j++){//j is a passenger node or parcel node (pu or del)
                    if (i != j){
                        if (j + inst->m != i){
                            nas->arcs[i][j] = true;
                            nas->fArc.first = i;
                            nas->fArc.second = j;
                            nas->arcMinus[j].push_back(nas->fArc);
                            nas->arcPlus[i].push_back(nas->fArc);
                            nas->allArcs.push_back(nas->fArc);
                            nas->arcnf.push_back(nas->fArc);
                            for (int k = 0; k < inst->K; k++){
                                nas->arcV[i][j].push_back(k);
                            }
                        }                         
                    }
                }

                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->V + inst->dummy){ // i is a starting node
                for (int j = 0; j < inst->n + inst->m; j++){//j is a passenger or parcel pickup node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    nas->arcnf.push_back(nas->fArc);
                    auxK = i - fDepot;
                    nas->arcV[i][j].push_back(auxK);
                }
            }          
        }
    }
    else if (problem->scen == "P"){//serving only parcels
        for (int i = inst->n; i < inst->V; i++){//skip passenger nodes
            if (i < inst->n + inst->m){//i is a parcel pickup node
                for (int j = inst->n; j < inst->n + 2*inst->m; j++){ //j is a passenger or parcel pickup node
                    if (i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        nas->arcPN.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }                           
                    }
                }
            }
            else if (i < inst->n + 2*inst->m){// i is a parcel delivery node
                for (int j = inst->n; j < inst->n + 2*inst->m; j++){//j is a passenger node or parcel node (pu or del)
                    if (i != j){
                        if (j + inst->m != i){
                            nas->arcs[i][j] = true;
                            nas->fArc.first = i;
                            nas->fArc.second = j;
                            nas->arcMinus[j].push_back(nas->fArc);
                            nas->arcPlus[i].push_back(nas->fArc);
                            nas->allArcs.push_back(nas->fArc);
                            nas->arcnf.push_back(nas->fArc);
                            for (int k = 0; k < inst->K; k++){
                                nas->arcV[i][j].push_back(k);
                            }                             
                        }
                    }
                }

                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->V + inst->dummy){ // i is a starting node
                for (int j = inst->n; j < inst->n + inst->m; j++){//j is a parcel node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    nas->arcnf.push_back(nas->fArc);
                    auxK = i - fDepot;
                    nas->arcV[i][j].push_back(auxK);
                }
                int j = i + inst->K;
                nas->arcs[i][j] = true;
                nas->fArc.first = i;
                nas->fArc.second = j;
                nas->arcMinus[j].push_back(nas->fArc);
                nas->arcPlus[i].push_back(nas->fArc);
                nas->allArcs.push_back(nas->fArc);
                auxK = i - fDepot;
                nas->arcV[i][j].push_back(auxK);
            }          
        }        
    }

    else if (problem->scen == "C"){//serving only customers
         for (int i = 0; i < inst->V; i++){
            if (i < inst->n){//i is a passenger req
                for(int j = 0; j < inst->n; j++){// j is a pass req
                    if(i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        if (j < inst->n){
                            nas->arcNN.push_back(nas->fArc);
                        }
                        nas->arcNplus.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }                        
                    }
                }
                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->arcNplus.push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i >= inst->n + 2*inst->m){
                if (i < inst->V + inst->dummy){ // i is a starting node
                    for (int j = 0; j < inst->n; j++){//j is a passenger or parcel pickup node
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        auxK = i - fDepot;
                        nas->arcV[i][j].push_back(auxK);
                    }
                }    
            }
        }
    }

    else if (problem->scen == "PC"){//serving both parcels and passengers but one type on each vehicle
        for (int i = 0; i < inst->V; i++){
            if (i < inst->n){//i is a passenger req
                for(int j = 0; j < inst->n; j++){// j is a passenger req
                    if(i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        nas->arcNN.push_back(nas->fArc);

                        nas->arcNplus.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }                        
                    }
                }
                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->arcNplus.push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->n + inst->m){//i is a parcel pickup node
                for (int j = inst->n; j < inst->n + 2*inst->m; j++){ //j is a parcel node (pu or del)
                    if (i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        nas->arcPN.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        } 
                    }                 
                }
            }
            else if (i < inst->n + 2*inst->m){// i is a parcel delivery node
                for (int j = inst->n; j < inst->n + 2*inst->m; j++){//j is a parcel node (pu or del)
                    if (i != j){
                        if (j + inst->m != i){
                            nas->arcs[i][j] = true;
                            nas->fArc.first = i;
                            nas->fArc.second = j;
                            nas->arcMinus[j].push_back(nas->fArc);
                            nas->arcPlus[i].push_back(nas->fArc);
                            nas->allArcs.push_back(nas->fArc);
                            nas->arcnf.push_back(nas->fArc);
                            for (int k = 0; k < inst->K; k++){
                                nas->arcV[i][j].push_back(k);
                            }
                        }                         
                    }
                       
                }

                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->V + inst->dummy){ // i is a starting node
                for (int j = 0; j < inst->n + inst->m; j++){//j is a passenger or parcel pickup node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    nas->arcnf.push_back(nas->fArc);
                    auxK = i - fDepot;
                    nas->arcV[i][j].push_back(auxK);
                }
                int j = i + inst->K;
                nas->arcs[i][j] = true;
                nas->fArc.first = i;
                nas->fArc.second = j;
                nas->arcMinus[j].push_back(nas->fArc);
                nas->arcPlus[i].push_back(nas->fArc);
                nas->allArcs.push_back(nas->fArc);
                auxK = i - fDepot;
                nas->arcV[i][j].push_back(auxK);
            }          
        }
    }

    else if (problem->scen == "BL2"){//serving both parcels and passengers but the latter only when vehicle is empty
        for (int i = 0; i < inst->V; i++){
            if (i < inst->n){//i is a passenger req
                for(int j = 0; j < inst->n + inst->m; j++){// j is a passenger or a parcel pu
                    if(i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        if (j < inst->n){
                            nas->arcNN.push_back(nas->fArc);
                        }
                        nas->arcNplus.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }                        
                    }
                }
                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->arcNplus.push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->n + inst->m){//i is a parcel pickup node
                for (int j = inst->n; j < inst->n + 2*inst->m; j++){ //j is a passenger or parcel node
                    if (i != j){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        nas->arcPN.push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }                          
                    }                
                }
            }
            else if (i < inst->n + 2*inst->m){// i is a parcel delivery node
                for (int j = 0; j < inst->n + 2*inst->m; j++){//j is a passenger node or parcel node (pu or del)
                    if (j + inst->m != i){
                        nas->arcs[i][j] = true;
                        nas->fArc.first = i;
                        nas->fArc.second = j;
                        nas->arcMinus[j].push_back(nas->fArc);
                        nas->arcPlus[i].push_back(nas->fArc);
                        nas->allArcs.push_back(nas->fArc);
                        nas->arcnf.push_back(nas->fArc);
                        for (int k = 0; k < inst->K; k++){
                            nas->arcV[i][j].push_back(k);
                        }
                    }                    
                }

                for (int j = inst->V; j < inst->V + inst->dummy; j++){//j is the dummy node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    auxK = j - inst->V;
                    nas->arcV[i][j].push_back(auxK);
                }
            }

            else if (i < inst->V + inst->dummy){ // i is a starting node
                for (int j = 0; j < inst->n + inst->m; j++){//j is a passenger or parcel pickup node
                    nas->arcs[i][j] = true;
                    nas->fArc.first = i;
                    nas->fArc.second = j;
                    nas->arcMinus[j].push_back(nas->fArc);
                    nas->arcPlus[i].push_back(nas->fArc);
                    nas->allArcs.push_back(nas->fArc);
                    nas->arcnf.push_back(nas->fArc);
                    auxK = i - fDepot;
                    nas->arcV[i][j].push_back(auxK);
                }
            }          
        }
    }

    for (int a = 0; a < nas->allArcs.size(); a++){
        int i = nas->allArcs[a].first;
        int j = nas->allArcs[a].second;

        for(int k1 = 0; k1 < nas->arcV[i][j].size(); k1++){
            int k = nas->arcV[i][j][k1];
            nas->vArcPlus[i][k].push_back(nas->allArcs[a]);
            nas->vArcMinus[j][k].push_back(nas->allArcs[a]);
        }

    }

    // cout << "allowed vehicles: " << endl;

    // for (int i = 0; i < nas->arcV.size(); i++){
    //     for (int j = 0; j < nas->arcV[i].size(); j++){
    //         if (nas->arcs[i][j]){
    //             printf("\narc %i, %i, Vehicles: ", i, j);
    //             for (int k = 0; k < nas->arcV[i][j].size(); k++){
    //                 cout << nas->arcV[i][j][k] << endl;
    //             }
    //             cout << endl;                
    //         }
    //     }
    // }
    // getchar();
}

void viewSol (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, solStats *sStat){
    bool inserted;

    vector< pair <int, int> > auxVec;
    pair<int, int> auxPair;
    vector<int> auxSolOrder;
    // int setN = bStat->bundleVec.size() - inst->K - 1;
    int currSP;
    vector<int> orderVec;

	// solStatIni(sStat);

    for (int k = 0; k < inst->K; k++){
        sStat->solOrder.push_back(auxSolOrder);
    }

    for (int k = 0; k < inst->K; k++){
        currSP = inst->V - inst->K + k;

        for (int i = 0; i < sStat->solvec[k].size(); i++){
            auxPair.first = sStat->solvec[k][i].first;
            auxPair.second = sStat->solvec[k][i].second;            
            auxVec.push_back(auxPair);
        }
        // cout<< "here1";
        // getchar();
        // cout << "auxVec: " << endl;
        // for (int i = 0; i < auxVec.size(); i++){
        //     cout << auxVec[i].first << " " << auxVec[i].second << endl;
        // }

        while(!auxVec.empty()){
            if (sStat->solOrder[k].empty()){

                for (int i = 0; i < auxVec.size(); i++){
                    if (auxVec[i].first == currSP){
                        sStat->solOrder[k].push_back(auxVec[i].first);
                        sStat->solOrder[k].push_back(auxVec[i].second);

                        auxVec.erase(auxVec.begin()+i);

                    }
                }
            }
            else{

                for (int j = 0; j < auxVec.size(); j++){
                    if(auxVec[j].first == sStat->solOrder[k].back()){
                        sStat->solOrder[k].push_back(auxVec[j].second);

                        auxVec.erase(auxVec.begin()+j);
                    }
                }
            }       
        // cout<< "auxvec size: " << auxVec.size();
        // getchar();
        }
        // cout<< "here3";
        // getchar();
    }

    cout << "\nNumber of Vehicles: " << inst->K << endl;

    cout << "\nSolution: " << endl;
    for (int k = 0; k < inst->K; k++){
        cout << "Vehicle " << k << ": ";
        for (int i = 0; i < sStat->solOrder[k].size(); i++){
            if (i < sStat->solOrder[k].size() - 1){
                cout << sStat->solOrder[k][i] << " - ";
            }
            else{
                cout << sStat->solOrder[k][i];
            }
        }
        cout << endl;
    }
    cout << endl;

    cout << "\nSolution structure: " << endl;
    for (int k = 0; k < inst->K; k++){
        cout << "Vehicle " << k << ": ";
        for (int i = 0; i < sStat->solOrder[k].size(); i++){
            if (i < sStat->solOrder[k].size() - 1){
                if (sStat->solOrder[k][i] < inst->n){
                    cout << "d" << " - ";
                }
                else if (sStat->solOrder[k][i] < inst->n + inst->m){
                    cout << "P" << " - ";
                    sStat->servedParcels++;
                }
                else if (sStat->solOrder[k][i] < inst->n + 2*inst->m){
                    cout << "D" << " - ";
                }
                else if (sStat->solOrder[k][i] < inst->n + 2*inst->m + inst->K){
                    cout << "S" << " - ";
                }                                      
            }
            else{

                cout << "f";
            }
        }
        cout << endl;
    }
    cout << endl;   
    // getchar();    
}

void output(instanceStat *inst, vector<nodeStat> &nodeVec,  solStats *sStat, probStat* problem){

    //output
    string btoutputname;

    btoutputname = "bt-" + inst->InstName + "-" + problem->scen + ".txt";
    // cout << "output bt: " << btoutputname << endl;

    ofstream ofile;

    ofile.open(btoutputname);
    
    // ofile << K << "\t" << 5 << "\t" << dimVec[i].first << "\t" << dimVec[i].second << endl;

    for (int i = 0; i < inst->n; i++){
        ofile << i << "\t" << setw(9) << fixed << setprecision(4) << sStat->solBegin[i] << endl;
    }
    // for (int i = 0; i < inst->n; i++){
    //     ofile << i << "\t" << setw(9) << fixed << setprecision(4) << sStat->solBegin[i] << endl;

    // }
}

void nodeMethod (nodeStat *node, instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, probStat* problem, solStats *sStat){
	
	nodeArcsStruct nas;
	
	// for (int i = 0; i < inst->V; i++){
	// 	cout << "load " << i << ": " << nodeVec[i].load << endl;
	// }
	
	// for (int i = 0; i < inst->n; i++){
	// 	cout << "delta " << i << ": " << nodeVec[i].delta << endl;
	// }

    // cout << "\nDistance Matrix: " << endl;

    // for (int i = 0; i < inst->V + inst->dummy; i++){
    // 	for (int j = 0; j < inst->V + inst->dummy; j++){
    // 		cout << setw(8) << setprecision(5) << mdist[i][j] << " ";
    // 	}
    // 	cout << endl;
    // }
    // getchar();

    // cout << "\nDelta vector: " << endl;

    // for (int i = 0; i < nodeVec.size(); i++){
    //     cout << i << ": " << nodeVec[i].delta << endl;
    // }
    // cout << endl;
    // getchar();

    // cout << "\nProfit vector: " << endl;

    // for (int i = 0; i < nodeVec.size(); i++){
    //     cout << i << ": " << nodeVec[i].profit << endl;
    // }
    // cout << endl;
    // getchar();

	initArcs(inst, &nas);
	feasibleArcs (inst, &nas, problem, nodeVec, mdist);

	cout<< "\nFeasible arcs between nodes:" << endl;
    for (int i = 0; i < nas.arcs.size(); i++){
        if (i == 0){
            cout << setw(3) << " ";
        }
        cout << setw(3) << std::right << i << " ";
    }
    cout << endl;
    for(int i = 0; i < nas.arcs.size(); i++){
        cout << setw(3) << std::right << i;
        for(int j = 0; j < nas.arcs[i].size(); j++){
            cout << setw(3) <<  nas.arcs[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;
    getchar();

 //    cout << "arcs NN: " << endl;
	// for (int i = 0; i < nas.arcNN.size(); i++){
	// 	cout << nas.arcNN[i].first << " - " << nas.arcNN[i].second << " | | ";
	// }
	// getchar();
	// cout << "Arcs that leave a pickup: " << endl;

	// for (int i = 0; i < nas.arcPN.size(); i++){
	// 	cout<< nas.arcPN[i].first << "-" << nas.arcPN[i].second << "  |  ";
	// }
	// cout << endl;

    mipnode(inst, nodeVec, mdist, problem, &nas, sStat);
    
    // mtznode(inst, nodeVec, mdist, problem, &nas, sStat);

	if(sStat->feasible){
		viewSol (inst, mdist, nodeVec, sStat);

		mipSolStats (inst, mdist, nodeVec, sStat);

		printStats(inst, sStat);

        if (inst->preInst == 1) {
            output(inst, nodeVec,  sStat, problem);
        }
	}
    
	for ( int i = 0; i < inst->V + inst->dummy; i++) {
		delete[] mdist[i];
	}
	delete[] mdist;
}

