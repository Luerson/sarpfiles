#include "readdata.h"
#include "functions.h"
#include "modelnode.h"
#include "modeltwostage.h"
#include <cstdlib>
#include <stdio.h>
//maybe later it is necessary to add the maxiparkmm driving time.

void readData (int argc, char** argv, nodeStat *node, instanceStat *inst, vector<nodeStat> &nodeVec, double ***Mdist, probStat* problem, int trialK, double trialMulti){
    
    if (argc < 5) {
        cout << "\nMissing parameters\n";
        // cout << " ./exeSARP [Instance] [Optimization strategy] [Scenario]"<< endl;
        cout << " ./exeSARP [Instance] [Scenario] [Parcel percentage] [model type]"<< endl;
        exit(1);
    }
    
    if (argc > 5) {
        cout << "\nToo many parameters\n";
        cout << " ./exeSARP [Instance] [Scenario] [Parcel percentage] [model type]" << endl;
        exit(1);
    }  


    inst->preInst = 0;
    inst->InstName = getInstName(argv);

    // cout << "Instance Name: " << inst->InstName;
    // getchar();
    // if (argv[2] == "sim"){
    //     problem->sim = true;
    // }
    // else if (argv[2] == "seq"){
    //     problem->seq = true;
    // }

    problem->scen = argv[2];
    problem->model = argv[4];

    string file, ewf;
    int n;
    int m;
    int K;
    double service;
    double T;
    int V;
    int originalV;
    int dummy;
    string parcelArg;
    double parcelP;

    // string instType;

    char *instance; 
    instance = argv[1];

    ifstream in(instance, ios::in);

    
    if( !in ) {
        cout << "the file could not be opened\n";
        exit (1);
    }

    inst->instType = getInstanceType(argv);
    parcelArg = argv[3];
    parcelP = stod(parcelArg)/100; //string to double
    
    if (inst->instType == "sf_data"){

        in >> K;
        in >> service;
        in >> n;
        in >> m;

        // K = 11;

        if (trialK >= K){
            K = trialK;
        }
        else{
            trialK = K;
        }
        
        cout << "K: " << K << endl;
        // getchar();

        service = service/60;
        V = n + 2*m + K;

        originalV = 2*n + 2*m + 2; 
        // inst->vmed = 19.3;
        inst->dummy = 1;

        double *delta = new double[V + inst->dummy];
        double *slatitude = new double [V + inst->dummy];
        double *slongitude = new double [V + inst->dummy];
        double *flatitude = new double [V + inst->dummy];
        double *flongitude = new double [V + inst->dummy];
        double *profit = new double[V+inst->dummy];
        double *trip = new double[V+inst->dummy];

        double **dist = new double*[V + inst->dummy];
        for (int i= 0; i < V + inst->dummy; i++){
            dist[i] = new double [V + inst->dummy];
        }

        vector<double> vxs;
        vector<double> vys;
        vector<double> vload;
        vector<double> ve;
        vector<double> vxf;
        vector<double> vyf;
        vector<double> vl;

        int tempNode;

        for (int i = 0; i < originalV; i++){
            vxs.push_back(0);
            vys.push_back(0);
            vload.push_back(0);
            ve.push_back(0);
            vl.push_back(0);
        }

        for (int i = 0; i < originalV; i++){
            in >> tempNode >> vxs[i] >> vys[i] >> tempNode >> vload[i] >> ve[i] >> vl[i];
        }

        // ve[ve.size()-1] = 0;

        for (int i = 0; i < vxs.size(); i++){
            vxf.push_back(vxs[i]);
            vyf.push_back(vys[i]);

            if (vload[i] < -2.0){
                vxf[i - m - n] = vxs[i];
                vyf[i - m - n] = vys[i];
            }
        }

        vxs.erase(vxs.begin());
        vys.erase(vys.begin());
        vload.erase(vload.begin());
        vxf.erase(vxf.begin());
        vyf.erase(vyf.begin());
        ve.erase(ve.begin());
        vl.erase(vl.begin());

        for (int i = 0; i < n; i++){
            vxs.erase(vxs.begin() + n + m);
            vys.erase(vys.begin() + n + m);
            vload.erase(vload.begin() + n + m);
            ve.erase(ve.begin() + n + m);
            vl.erase(vl.begin() + n + m);
            vxf.erase(vxf.begin() + n + m);
            vyf.erase(vyf.begin() + n + m);
        }

        for (int i = 0; i < n; i++){
            vl[i] = ve[i];
        }

        for (int i = n; i < n + 2*m; i++){
            ve[i] = 0;
   
        }

        for (int i = 1; i < K; i++){
            vxs.push_back(vxs[vxs.size()-1]);
            vys.push_back(vys[vys.size()-1]);
            vload.push_back(vload[vload.size()-1]);
            ve.push_back(ve[ve.size()-1]);
            vl.push_back(vl[vl.size()-1]);
            vxf.push_back(vxf[vxf.size()-1]);
            vyf.push_back(vyf[vyf.size()-1]);
        }

        for (int i = 0; i < n; i++){
            vload[i] = 0;
        }


        cout << "\nve: " << endl;
        for (int i = 0; i < ve.size(); i++){
            cout << ve[i] << endl;
        }
        // getchar();

        // Calculate distance matrix (Geolocation)

        CalcLatLong ( vxs, vys, vxf, vyf, V, slatitude, slongitude, flatitude, flongitude );

        double singleProfit;
        for (int i = 0; i < V + inst->dummy; i++){
            if (i < n){ 
                delta[i] = (2 * (service)) + (CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i))/inst->vmed;
                profit[i] = inst->minpas + inst->paskm*CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i) - inst->costkm*CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i);
            }
            else if (i < V - K){ 
                delta[i] = service;
                if (i < n + m){
                    profit[i] = inst->minpar + inst->parkm*CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i+m);
                }
                else{
                    profit[i] = 0;
                }
            }
            else if (i >= V - K){
                delta[i] = 0;
                profit[i] = 0;
            }
            for (int j = 0; j < V + inst->dummy; j++){
                if(i == j){
                   dist[i][j] = 0;
                }
                else{
                    if (i < V){
                        if (j < V){
                            dist[i][j] = CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, j);
                        }
                        else if (j >= V){
                            dist[i][j] = 0;
                        }
                    }
                    else{
                        dist[i][j] = 0;
                    }
                }
            }
        }

        for (int i = 0; i < V; i++){
            node->xs = vxs[i];
            node->ys = vys[i];
            node->load = vload[i];
            node->e = ve[i]/60;
            node->l = vl[i]/60;
            node->xf = vxf[i];
            node->yf = vyf[i];
            node->delta = delta[i];
            node->profit = profit[i];
            node->index = i;
            nodeVec.push_back(*node);
        }

        // Adding dummy nodes
        for (int i = 0; i < inst->dummy; i++){
            node->xs = 0;
            node->ys = 0;
            node->load = 0;
            node->e = 0;
            // node->l = 14*60;
            node->l = 24;
            node->xf = 0;
            node->yf = 0;
            node->delta = 0;
            node->profit = 0;
            node->index = V + i;
            nodeVec.push_back(*node);
        }
        

        // cout << "Earlier: \n{";

        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << nodeVec[i].e << " } {";
        // }
        // cout << endl;
        // getchar();

        *Mdist = dist;
        inst->K = K;
        inst->n = n;
        inst->m = m;
        inst->V = V;
        inst->service = service;
        // inst->T = nodeVec[V + inst->dummy - 1].l;
        // inst->T = 1440/60;

        inst->totalCustomProfit = 0;

        // if (problem->model == "node"){
        //     for (int i = 0; i < n; i++){
        //         inst->totalCustomProfit += nodeVec[i].profit;
        //     }
        // }


        for (int i = 0; i < n; i++){
            inst->totalCustomProfit += nodeVec[i].profit;
        }

        // cout << "\ntotal profit: " << inst->totalCustomProfit;

        // getchar();

        delete[] profit;
        delete[] delta;
        delete[] slatitude;
        delete[] slongitude;
        delete[] flatitude;
        delete[] flongitude;   
        delete[] trip;
    }

    if (inst->instType == "debug"){

        in >> K;
        in >> service;
        in >> n;
        in >> m;

        // K = 11;

        if (trialK >= K){
            K = trialK;
        }
        else{
            trialK = K;
        }
        
        // cout << "K: " << K << endl;
        // getchar();

        service = service/60;
        V = n + 2*m + K;

        originalV = 2*n + 2*m + 2; 
        // inst->vmed = 19.3;
        inst->dummy = 1;

        double *delta = new double[V + inst->dummy];
        double *slatitude = new double [V + inst->dummy];
        double *slongitude = new double [V + inst->dummy];
        double *flatitude = new double [V + inst->dummy];
        double *flongitude = new double [V + inst->dummy];
        double *profit = new double[V+inst->dummy];
        double *trip = new double[V+inst->dummy];

        double **dist = new double*[V + inst->dummy];
        for (int i= 0; i < V + inst->dummy; i++){
            dist[i] = new double [V + inst->dummy];
        }

        vector<double> vxs;
        vector<double> vys;
        vector<double> vload;
        vector<double> ve;
        vector<double> vxf;
        vector<double> vyf;
        vector<double> vl;

        int tempNode;

        for (int i = 0; i < originalV; i++){
            vxs.push_back(0);
            vys.push_back(0);
            vload.push_back(0);
            ve.push_back(0);
            vl.push_back(0);
        }

        for (int i = 0; i < originalV; i++){
            in >> tempNode >> vxs[i] >> vys[i] >> tempNode >> vload[i] >> ve[i] >> vl[i];
        }

        // ve[ve.size()-1] = 0;

        for (int i = 0; i < vxs.size(); i++){
            vxf.push_back(vxs[i]);
            vyf.push_back(vys[i]);

            if (vload[i] < -2.0){
                vxf[i - m - n] = vxs[i];
                vyf[i - m - n] = vys[i];
            }
        }

        vxs.erase(vxs.begin());
        vys.erase(vys.begin());
        vload.erase(vload.begin());
        vxf.erase(vxf.begin());
        vyf.erase(vyf.begin());
        ve.erase(ve.begin());
        vl.erase(vl.begin());

        for (int i = 0; i < n; i++){
            vxs.erase(vxs.begin() + n + m);
            vys.erase(vys.begin() + n + m);
            vload.erase(vload.begin() + n + m);
            ve.erase(ve.begin() + n + m);
            vl.erase(vl.begin() + n + m);
            vxf.erase(vxf.begin() + n + m);
            vyf.erase(vyf.begin() + n + m);
        }

        for (int i = 0; i < n; i++){
            vl[i] = ve[i];
        }

        for (int i = n; i < n + 2*m; i++){
            ve[i] = 0;
   
        }

        for (int i = 1; i < K; i++){
            vxs.push_back(vxs[vxs.size()-1]);
            vys.push_back(vys[vys.size()-1]);
            vload.push_back(vload[vload.size()-1]);
            ve.push_back(ve[ve.size()-1]);
            vl.push_back(vl[vl.size()-1]);
            vxf.push_back(vxf[vxf.size()-1]);
            vyf.push_back(vyf[vyf.size()-1]);
        }

        for (int i = 0; i < n; i++){
            vload[i] = 0;
        }

        // cout << "Vector of vxf:" << endl;

        // for (int i = 0; i < vxf.size(); i++){
        //     cout << i << ": " << vxf[i] << endl;
        // }
        // cout << endl;

        // getchar();

        // cout << "\nve: " << endl;
        // for (int i = 0; i < ve.size(); i++){
        //     cout << ve[i] << endl;
        // }
        // getchar();

        // Calculate distance matrix (Geolocation)

        // CalcLatLong ( vxs, vys, vxf, vyf, V, slatitude, slongitude, flatitude, flongitude );

        double singleProfit;
        for (int i = 0; i < V + inst->dummy; i++){
            if (i < n){ 
                delta[i] = (2 * (service)) + (calcEucDist(vxs, vys, vxf, vyf, i, i))/inst->vmed;

                trip[i] = calcEucDist(vxs, vys, vxf, vyf, i, i);
                // cout << "Trip:" << endl;
                // cout << i << ": " << trip[i] << endl;
                profit[i] = inst->minpas + inst->paskm*calcEucDist(vxs, vys, vxf, vyf, i, i) - inst->costkm*calcEucDist(vxs, vys, vxf, vyf, i, i);
                // cout << "Profit: " << endl;
                // cout << i << ": " << profit[i] << endl;
            }
            else if (i < V - K){ 
                delta[i] = service;
                if (i < n + m){
                    profit[i] = inst->minpar + inst->parkm*calcEucDist(vxs, vys, vxf, vyf, i, i+m);
                }
                else{
                    profit[i] = 0;
                }
            }
            else if (i >= V - K){
                delta[i] = 0;
                profit[i] = 0;
            }
            for (int j = 0; j < V + inst->dummy; j++){
                if(i == j){
                   dist[i][j] = 0;
                }
                else{
                    if (i < V){
                        if (j < V){
                            dist[i][j] = calcEucDist(vxs, vys, vxf, vyf, i, j);
                        }
                        else if (j >= V){
                            dist[i][j] = 0;
                        }
                    }
                    else{
                        dist[i][j] = 0;
                    }
                }
            }
        }
        // getchar();
        for (int i = 0; i < V; i++){
            node->xs = vxs[i];
            node->ys = vys[i];
            node->load = vload[i];
            node->e = ve[i]/60;
            node->l = vl[i]/60;
            node->xf = vxf[i];
            node->yf = vyf[i];
            node->delta = delta[i];
            node->profit = profit[i];
            node->index = i;
            nodeVec.push_back(*node);
        }

        // Adding dummy nodes
        for (int i = 0; i < inst->dummy; i++){
            node->xs = 0;
            node->ys = 0;
            node->load = 0;
            node->e = 9;
            // node->l = 14*60;
            node->l = 24;
            node->xf = 0;
            node->yf = 0;
            node->delta = 0;
            node->profit = 0;
            node->index = V + i;
            nodeVec.push_back(*node);
        }
        
        // cout << "Earlier: \n{";

        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << nodeVec[i].e << " } {";
        // }
        // cout << endl;
        // getchar();

        *Mdist = dist;
        inst->K = K;
        inst->n = n;
        inst->m = m;
        inst->V = V;
        inst->service = service;
        // inst->T = nodeVec[V + inst->dummy - 1].l;
        // inst->T = 1440/60;

        inst->totalCustomProfit = 0;

        // if (problem->model == "node"){
        //     for (int i = 0; i < n; i++){
        //         inst->totalCustomProfit += nodeVec[i].profit;
        //     }
        // }
        for (int i = 0; i < n; i++){
            inst->totalCustomProfit += nodeVec[i].profit;
        }
        // cout << "\ntotal profit: " << inst->totalCustomProfit;

        // getchar();

        delete[] profit;
        delete[] delta;
        delete[] slatitude;
        delete[] slongitude;
        delete[] flatitude;
        delete[] flongitude;   
        delete[] trip;
    }


    else if (inst->instType == "csarp"){
        in >> K;
        in >> service;
        in >> n;
        in >> m;

        inst->preInst = 1;
        // K = 11;

        if (trialK < K){
            K = trialK;
        }
        else{
            trialK = K;
        }

        cout << "K: " << K << endl;
        // getchar();
        service = service/60;
        V = n + 2*m + K;

        originalV = 2*n + 2*m + K; 
        // inst->vmed = 19.3;
        inst->dummy = 1;

        double *delta = new double[V + inst->dummy];
        double *profit = new double[V+inst->dummy];
        double *trip = new double[V+inst->dummy];

        double **dist = new double*[V + inst->dummy];
        for (int i= 0; i < V + inst->dummy; i++){
            dist[i] = new double [V + inst->dummy];
        }

        vector<double> vxs;
        vector<double> vys;
        vector<double> vload;
        vector<double> ve;
        vector<double> vxf;
        vector<double> vyf;
        vector<double> vl;

        int tempNode;

        for (int i = 0; i < originalV; i++){
            vxs.push_back(0);
            vys.push_back(0);
            vload.push_back(0);
            ve.push_back(0);
            vl.push_back(0);
        }

        for (int i = 0; i < originalV; i++){
            in >> tempNode >> vxs[i] >> vys[i] >> vload[i] >> ve[i] >> vl[i];
        }

        // ve[ve.size()-1] = 0;

        for (int i = 0; i < vxs.size(); i++){
            vxf.push_back(vxs[i]);
            vyf.push_back(vys[i]);

            if (vload[i] < -2.0){
                vxf[i - n] = vxs[i];
                vyf[i - n] = vys[i];
            }
        }

        for (int i = 0; i < n; i++){
            vxs.erase(vxs.begin() + n);
            vys.erase(vys.begin() + n);
            vload.erase(vload.begin() + n);
            ve.erase(ve.begin() + n);
            vl.erase(vl.begin() + n);
            vxf.erase(vxf.begin() + n);
            vyf.erase(vyf.begin() + n);
        }


        // for (int i = 0; i < n; i++){
        //     vl[i] = ve[i];
        // }

        for (int i = n; i < n + 2*m; i++){
            ve[i] = 0;
   
        }

        for (int i = n + m; i < n + 2*m; i++){
            vl[i] = 1440;
   
        }

        // for (int i = 1; i < K; i++){
        //     vxs.push_back(vxs[vxs.size()-1]);
        //     vys.push_back(vys[vys.size()-1]);
        //     vload.push_back(vload[vload.size()-1]);
        //     ve.push_back(ve[ve.size()-1]);
        //     vl.push_back(vl[vl.size()-1]);
        //     vxf.push_back(vxf[vxf.size()-1]);
        //     vyf.push_back(vyf[vyf.size()-1]);
        // }

        // cout << "Vector of vxf:" << endl;

        // for (int i = 0; i < vxf.size(); i++){
        //     cout << i << ": " << vxf[i] << endl;
        // }
        // cout << endl;       

        // getchar();

        for (int i = 0; i < n; i++){
            vload[i] = 0;
        }


        // cout << "\nve: " << endl;
        // for (int i = 0; i < ve.size(); i++){
        //     cout << ve[i] << endl;
        // }
        // getchar();

        // Calculate distance matrix (Euclidian)

        double singleProfit;
        for (int i = 0; i < V + inst->dummy; i++){
            if (i < n){ 
                delta[i] = (2 * (service)) + (calcEucDist(vxs, vys, vxf, vyf, i, i))/inst->vmed;
                profit[i] = inst->minpas + inst->paskm*calcEucDist(vxs, vys, vxf, vyf, i, i) - inst->costkm*calcEucDist(vxs, vys, vxf, vyf, i, i);
            }
            else if (i < V - K){ 
                delta[i] = service;
                if (i < n + m){
                    profit[i] = inst->minpar + inst->parkm*calcEucDist(vxs, vys, vxf, vyf, i, i+m);
                }
                else{
                    profit[i] = 0;
                }
            }
            else if (i >= V - K){
                delta[i] = 0;
                profit[i] = 0;
            }
            for (int j = 0; j < V + inst->dummy; j++){
                if(i == j){
                   dist[i][j] = 0;
                }
                else{
                    if (i < V){
                        if (j < V){
                            dist[i][j] = calcEucDist(vxs, vys, vxf, vyf, i, j);
                        }
                        else if (j >= V){
                            dist[i][j] = 0;
                        }
                    }
                    else{
                        dist[i][j] = 0;
                    }
                }
            }
        }

        for (int i = 0; i < V; i++){
            node->xs = vxs[i];
            node->ys = vys[i];
            node->load = vload[i];
            node->e = ve[i]/60;
            node->l = vl[i]/60;
            node->xf = vxf[i];
            node->yf = vyf[i];
            node->delta = delta[i];
            node->profit = profit[i];
            node->index = i;
            nodeVec.push_back(*node);
        }

        // Adding dummy nodes
        for (int i = 0; i < inst->dummy; i++){
            node->xs = 0;
            node->ys = 0;
            node->load = 0;
            node->e = 9;
            // node->l = 14*60;
            node->l = 24;
            node->xf = 0;
            node->yf = 0;
            node->delta = 0;
            node->profit = 0;
            node->index = V + i;
            nodeVec.push_back(*node);
        }
        
        // cout << "Earlier // Later: " << endl;

        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << i << ": {" << nodeVec[i].e << "}-{" << nodeVec[i].l << "}" << endl;
        // }
        // cout << endl;
        // getchar();

        // cout << "Service times: " << endl;
        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << i << ": " << nodeVec[i].delta << endl;
        // }

        // getchar();

        // cout << "Profits: " << endl;
        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << i << ": " << nodeVec[i].profit << endl;
        // }

        // getchar();

        *Mdist = dist;
        inst->K = K;
        inst->n = n;
        inst->m = m;
        inst->V = V;
        inst->service = service;
        // inst->T = nodeVec[V + inst->dummy - 1].l;
        // inst->T = 1020/60;

        inst->totalCustomProfit = 0;

        // if (problem->model == "node"){
        //     for (int i = 0; i < n; i++){
        //         inst->totalCustomProfit += nodeVec[i].profit;
        //     }
        // }

        for (int i = 0; i < n; i++){
            inst->totalCustomProfit += nodeVec[i].profit;
        }
        // cout << "\ntotal profit: " << inst->totalCustomProfit;

        // getchar();

        delete[] profit;
        delete[] delta;
        delete[] trip;

    }

    else if (inst->instType == "grubhub"){

        K = 1;

        int scale;
        int scCounter = 0;

        // bool increaseK = false;
        // K = 2;
        service = 5; //for some reason, service = 5/60 did not work
        service = service/60;
        int refpoint = K + 1;
        int instV;
        dummy = 1;
        inst->dummy = dummy;
        // inst->vmed = 19.3; //(km/h)

        vector<int> instParam;

        getInstParam (inst, instParam);
        long power = pow(2, instParam[1]);
        // cout << "power: " << power << endl;
        long int seed;
        seed = (instParam[0]*power)*trialMulti;
        
        // if (instParam[0] < 14){
        //     if (instParam[0] > 12 && instParam[1] == 7){
        //         seed = (instParam[0]*power)*trialMulti;//old seed
        //     }
        //     else{
        //         seed = (instParam[0]*instParam[0]*power)*trialMulti;//new seed
        //     }                       
        // }

        // else if (instParam[0] > 14) {
        //     seed = (instParam[0]*power)*trialMulti;//old seed
        // }

        // else{
        //     if (instParam[1] > 1){
        //         seed = (instParam[0]*power)*trialMulti;//old seed         
        //     }
        //     else{
        //         seed = (instParam[0]*instParam[0]*power)*trialMulti;//new seed
        //     }
        // }

        for (int i = 0; i < instParam.size(); i++){
            inst->instParam.push_back(instParam[i]);
        }

        // int seed = time(NULL);
        srand(seed);
        // cout << "Instance seed: " << seed << endl;
        // getchar();

        vector <vector <double> > tempData;
        vector<double> auxtempdata;

        vector <vector <double> > realData;

        while ( file.compare("DIMENSION:") != 0 && file.compare("DIMENSION") != 0 ){
            in >> file;
        }
        
        in >> instV;

        m = floor(instV * parcelP);
        
        if (m % 2 != 0){
            m--;
        } 

        n = (instV - refpoint - m)/2;

        m /= 2; 

        V = n + 2*m + K;
        // cout << "V: " << V << endl;
        cout << "\nn: " << n << " m: " << m << endl; 


        while ( file.compare("EDGE_WEIGHT_FORMAT") != 0){
            in >> file;
        }

        in >> file;
        in >> ewf;

        while (file.compare("EDGE_WEIGHT_SECTION") != 0){
            in >> file;
        }
        
        for (int i = 0; i < instV + refpoint; i++){
            for(int j = 0; j < instV + refpoint; j++){
                auxtempdata.push_back(0);
            }
            tempData.push_back(auxtempdata);
            auxtempdata.clear();
        }

        for (int i = 0; i < V + inst->dummy; i++){
            for(int j = 0; j < V + inst->dummy; j++){
                auxtempdata.push_back(0);
            }
            realData.push_back(auxtempdata);
            auxtempdata.clear();
        }

        if (ewf == "LOWER_DIAG_ROW"){
           for (int i = 0; i < instV; i++) {
                for (int j = 0; j < i + 1; j++) {
                    in >> tempData[i][j];
                    if (i > 0){
                        tempData[j][i] = tempData[i][j];                        
                    }
                }
            }
        }
        
        //putting starting point in the end
        //adjusting rows 
        for (int i = 0; i < instV; i++){
            for (int j = 0; j < refpoint; j++){
                tempData[instV + j][i] = tempData[j][i];
            }
        }

        //adjusting columns
        for (int i = 0; i < instV; i++){
            for (int j = 0; j < refpoint; j++){
                tempData[i][instV + j] = tempData[i][j];
            }
        }


        // //maybe we needed to adjust the corner (relying on the -0 being f)

        //erase unused
        for (int j = 0; j < refpoint; j++){
            tempData.erase(tempData.begin());
        }

        for (int i = 0; i < instV; i++){
            for (int j = 0; j < refpoint; j++){
                tempData[i].erase(tempData[i].begin());
            }   
        }
        // cout << "\nDistance Matrix BEFORE (1): " << endl;
        // for (int i = 0; i < tempData.size(); i++){
        //     for (int j = 0; j < tempData[i].size(); j++){
        //         cout << setw(5) << tempData[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // getchar();

        tempData.erase(tempData.begin() + tempData.size() - 1);
        
        for (int i = 0; i < tempData.size(); i++){
            tempData[i].erase(tempData[i].begin() + tempData[i].size() - 1);
        }

        for(int i = 0; i < instV; i++){
            auxtempdata.push_back(0);
        }
        tempData.push_back(auxtempdata);
        auxtempdata.clear();

        for (int i = 0; i < tempData.size() - 1; i++){
            tempData[i].push_back(0);
        }

        // cout << "\nDistance Matrix BEFORE (2): " << endl;
        // for (int i = 0; i < tempData.size(); i++){
        //     for (int j = 0; j < tempData[i].size(); j++){
        //         cout << setw(5) << tempData[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // getchar();


        vector<double> dRow;
        vector<double> pRow;

        vector< vector<double> > deliveries;
        vector< vector<double> > pickups;

        for (int i = 0; i < 2*m; i++){
            for (int j = 0; j < tempData.size(); j++){
                if (i % 2 != 0){
                    dRow.push_back(tempData[2*n + i][j]);
                }
                else{
                    pRow.push_back(tempData[2*n + i][j]);
                }
            }
            if (i % 2 != 0){
                deliveries.push_back(dRow);
                dRow.clear();
            }
            else{
                pickups.push_back(pRow);
                pRow.clear();
            }
        }

        //Organize columns in dettached rows
        int count;
        double chosen;

        for (int j = 0; j < pickups.size(); j++){
            count = 1;
            for (int i = 2; i < 2*m; i++){
                if (i % 2 == 0){
                    // cout << "First i: " << i << "- j: " << j << endl;
                    // getchar();
                    chosen = pickups[j][2*n + i];
                    // getchar();
                    // cout << "chosen " << chosen << " i: " << i << "- j: " << j << endl;
                    // getchar();
                    pickups[j].erase(pickups[j].begin() + 2*n + i);
                    pickups[j].insert(pickups[j].begin() + 2*n + count, chosen);
                    // cout << "\nPickups verify: ";
                    // for (int k = 0; k < pickups[j].size(); k++){
                    //     cout << pickups[j][k] << " "; 
                    // }
                    // cout << endl;
                    // getchar();
                    count++;
                }
            }
        }


        for (int j = 0; j < deliveries.size(); j++){
            count = 1;
            for (int i = 2; i < 2*m; i++){
                if (i % 2 == 0){
                    chosen = deliveries[j][2*n + i];
                    deliveries[j].erase(deliveries[j].begin() + 2*n + i);
                    deliveries[j].insert(deliveries[j].begin() + 2*n + count, chosen);
                    count++;
                }
            }
        }

        for (int j = 0; j < tempData.size(); j++){
            for (int i = 0; i < 2*m; i++){
                tempData[j].erase(tempData[j].begin()+2*n); 
            }
        }

        tempData.erase(tempData.begin()+2*n, tempData.begin()+2*n+2*m);  

        for (int i = deliveries.size() - 1; i >= 0; i--){
            tempData.insert(tempData.begin()+2*n, deliveries[i]);
        }

        for (int i = pickups.size() - 1; i >= 0; i--){
            tempData.insert(tempData.begin()+2*n, pickups[i]);
        }

        for (int k = deliveries.size() - 1; k >= 0; k--){
            for (int j = 0; j < deliveries[k].size(); j++){
                if (j < 2*n){
                    tempData[j].insert(tempData[j].begin()+2*n, deliveries[k][j]);
                }
                else if (j < 2*n+2*m){
                    continue;
                }
                else{
                    tempData[j].insert(tempData[j].begin()+2*n, deliveries[k][j]);
                }
            }
        }

        for (int k = pickups.size() - 1; k >= 0; k--){
            for (int j = 0; j < pickups[k].size(); j++){
                if (j < 2*n){
                    tempData[j].insert(tempData[j].begin()+2*n, pickups[k][j]);
                }
                else if (j < 2*n+2*m){
                    continue;
                }
                else{
                    tempData[j].insert(tempData[j].begin()+2*n, pickups[k][j]);
                }
            }
        }
        
        // cout << "\nDistance Matrix (Pre-adapting): " << endl;

        // for (int i = 0; i < tempData.size(); i++){
        //     for (int j = 0; j < tempData[i].size(); j++){
        //         cout << setw(5) << setprecision(5) << tempData[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // getchar();

        double curAvg = 0;
        // double curStddv = 99999999999;
        scCounter = 0;

        // while (curAvg < inst->realAvg){
            // scale = 100 - 10*scCounter;
        scale = 50;
        // scale = 1000;
        // distScale(inst, &instV, tempData, &curAvg, &scale);

        //     scCounter++;
        // }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for (int i = 0; i < V + inst->dummy; i++){
            for (int j = 0; j < V + inst->dummy; j++){
                if (i == j){
                    realData[i][j] = 0;
                }
                else{
                    if (i < n){
                        if (j < n){
                            realData[i][j] = (tempData[2*i+1][2*j])/scale;
                        }
                        else{
                            realData[i][j] = (tempData[2*i+1][n+j])/scale;
                        }
                    }
                    else{
                        if (j < n){
                            realData[i][j] = (tempData[n+i][2*j])/scale;
                        }
                        else{
                            realData[i][j] = (tempData[n+i][n+j])/scale;
                        }
                    }
                }
            }
        }


        // cout << "Testing matrix(2): " << endl;
        // for (int i = 0; i < realData.size(); i++){
        //     for (int j = 0; j < realData[i].size(); j++){
        //         cout << setw(5) << realData[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // getchar();

        // if (trialK > 1){
        //     if (trialK >= K){
        //         K = trialK;
        //     }
        //     else{
        //         trialK = K;
        //     } 
        //     for (int l = 0; l < K -1; l++){
        //         vector<double> distRow;
        //         vector<double> dummyRow;

        //         double valueDist;
                
        //         for (int i = 0; i < V + inst->dummy; i++){
        //             valueDist = realData[i][realData[i].size() - 2];
        //             realData[i].insert(realData[i].begin() + realData[i].size() - 1, valueDist);
        //         }

        //         for (int i = 0; i < V + inst->dummy; i++){
        //             distRow.push_back(realData[V - 1][i]);
        //             dummyRow.push_back(realData[V][i]);
        //         }
        //         distRow.push_back(0);
        //         dummyRow.push_back(0); 

        //         realData.pop_back();

        //         realData.push_back(distRow);
        //         realData.push_back(dummyRow);
        //         V++;
        //     }
        // }

        vector<double> distRow;
        vector<double> dummyRow;
        vector<double> distColumn;

        // cout << "\nV: " << V << endl;
        // getchar();
        if (trialK > 1){
            if (trialK >= K){
                K = trialK;
            }
            else{
                trialK = K;
            }

            for (int i = 0; i < V - 1; i++){
                distColumn.push_back(realData[i][realData[i].size() - 2]);
                dummyRow.push_back(realData[V][i]);
                // if (i < V - 1){
                    distRow.push_back(realData[V - 1][i]);
                // }
                realData[i].pop_back();
                realData[i].pop_back();
            }

            realData.pop_back();
            realData.pop_back();
            V--;
            inst->dummy--;

            for (int l = 0; l < K; l++){
                distRow.push_back(0);
                dummyRow.push_back(0);
            }

            for (int l = 0; l < K; l++){
                for (int i = 0; i < distColumn.size(); i++){
                    realData[i].push_back(distColumn[i]);
                }

                realData.push_back(distRow);

                V++;
            }

            for (int l = 0; l < K; l++){
                for (int i = 0; i < dummyRow.size(); i++){
                    realData[i].push_back(dummyRow[i]);
                }
                // if (l > 0){
                    dummyRow.push_back(0);
                // }
                realData.push_back(dummyRow);
                inst->dummy++;
            }

        }

        // cout << "Testing matrix(3): " << endl;
        // for (int i = 0; i < realData.size(); i++){
        //     for (int j = 0; j < realData[i].size(); j++){
        //         cout << setw(5) << realData[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // getchar();

        double *delta = new double[V + inst->dummy];
        double *profit = new double[V + inst->dummy];
        double *e = new double[V + inst->dummy];
        double *l = new double[V + inst->dummy];
        int *w = new int[V + inst->dummy];      

        int reference = n;

        //calculate deltas
        for(int i = 0; i < V + inst->dummy; i++){
            if (i < n){
                // cout << i << ": " << (tempData[2*i][2*i+1]);
                // delta[i] = 2 * service + (((tempData[2*i][2*i+1])/1000)/inst->vmed);
                delta[i] = 2 * service + (((tempData[2*i][2*i+1])/scale)/inst->vmed);

                // cout << "i: " << i << " - " << ((tempData[2*i][2*i+1])/1000)/inst->vmed << endl;

                // profit[i] = inst->minpas + inst->paskm*(tempData[2*i][2*i+1]/1000) - inst->costkm*(tempData[2*i][2*i+1]/1000);    
                profit[i] = inst->minpas + inst->paskm*(tempData[2*i][2*i+1]/scale) - inst->costkm*(tempData[2*i][2*i+1]/scale);    
                w[i] = 0;
            }
            else if (i < V - K){
                delta[i] = service;
                if (i < n + m){
                    // profit[i] = inst->minpar + inst->parkm*(tempData[i + n][i + n + m]/1000);
                    profit[i] = inst->minpar + inst->parkm*(tempData[i + n][i + n + m]/scale);
                    w[i] = 1;
                }
                else if (i < n + 2*m){
                    profit[i] = 0;
                    w[i] = -1;
                }
                else{
                   profit[i] = 0;
                   w[i] = 0;
                }
                
            }
            else if (i >= V - K){
                delta[i] = 0;
                profit[i] = 0;
                w[i] = 0;
            }
        }

        for (int i = 0; i < V + inst->dummy; i++){

            if(i < n){
                e[i] = rand() % 1440;

                if (e[i] < 60){
                    e[i] += (rand() % 60)*2;
                }

                else if (e[i] + (delta[i]*60) > 1440){
                    e[i] -= delta[i]+(rand() % 60);
                }

                // e[i] += (rand()%240)/2;
                // if (e[i] + delta[i] > 1020){
                //     e[i] = 1020;
                // }
                // cout << "e[" << i << "]: " << e[i]/60 << endl;
                // getchar();
                l[i] = e[i];
                // e[i] = 540;
                // l[i] = 1020;
            }
            else if (i >= n && i < n + 2*m){
                e[i] = 0;
                l[i] = 1440;
            }
            else if (i >= n + 2*m && i < V + inst->dummy - 1){
                e[i] = 0;
                l[i] = 1440;
            }
            else if (i >= V + inst->dummy - 1){
                e[i] = 0;
                l[i] = 1440;
            } 
        }

        // for (int i = 0; i < V + inst->dummy; i++){
        //     if(i < n){
        //         if(i == 0){
        //             e[i] = 620;
        //             l[i] = e[i];                    
        //         }
        //         else if(i == 1){
        //             e[i] = 700;
        //             l[i] = e[i];                    
        //         }
        //     }
        //     else if (i < n + 2*m){
        //         e[i] = 0;
        //         l[i] = 1020;                
        //     }
        //     else if (i < V + inst->dummy - 1){
        //         e[i] = 540;
        //         l[i] = 1020;                    
        //     }
        //     else{
        //         e[i] = 540;
        //         l[i] = 1020;
        //     } 
        // }

        for (int i = 0; i < V + inst->dummy; i++){
            node->e = e[i]/60;
            node->l = l[i]/60;
            node->delta = delta[i];
            node->profit = profit[i];
            node->load = w[i];
            node->index = i;
            nodeVec.push_back(*node);
        }

        // cout << "\nLoad vector: " << endl;

        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << i << ": " << nodeVec[i].load << " ||| " << endl;
        // }
        // getchar();

        // cout << "\nDelta vector: " << endl;

        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << i << ": " << nodeVec[i].delta << " || ";
        // }
        // cout << endl;
        // getchar();

        double **dist = new double*[V + inst->dummy];
        for (int i= 0; i < V + inst->dummy; i++){
            dist[i] = new double [V + inst->dummy];
        }

        for(int i = 0; i < V + inst->dummy; i++){
            for (int j = 0; j < V + inst->dummy; j++){
                dist[i][j] = realData[i][j];
            }
        }

        *Mdist = dist;
        inst->K = K;
        inst->n = n;
        inst->m = m;
        inst->V = V;
        inst->service = service;
        // inst->T = nodeVec[V + inst->dummy - 1].l;
        // inst->T = 1440/60;

        inst->totalCustomProfit = 0;

        // if (problem->model == "node"){
        //     for (int i = 0; i < n; i++){
        //         inst->totalCustomProfit += nodeVec[i].profit;
        //     }
        // }

        for (int i = 0; i < n; i++){
            inst->totalCustomProfit += nodeVec[i].profit;
        }
        // cout << "\ntotal profit: " << inst->totalCustomProfit;

        // getchar();

        delete[] profit;
        delete[] delta;
        
        // cout << "\nStarting times: " << endl;

        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << i << ": " << nodeVec[i].e << " || ";
        // }
        // cout << endl << endl;
        // getchar();

    }

    //Print starting and end times: (debugging)
    // cout << "\nEarlier times: " << endl;

    // for (int i = 0; i < nodeVec.size(); i++){
    //     cout << i << ": " << nodeVec[i].e << endl;
    // }

    // cout << "\nLater times: " << endl;

    // for (int i = 0; i < nodeVec.size(); i++){
    //     cout << i << ": " << nodeVec[i].l << endl;
    // }
    cout << "\nDist Multiplier: " << trialMulti << endl;
    // getchar();
}

double calcEucDist (vector<double> &Xs, vector<double> &Ys, vector<double> &Xf, vector<double> &Yf, int I, int J){
    return sqrt(pow(Xf[I] - Xs[J], 2) + pow(Yf[I] - Ys[J], 2));
}

double CalcMan (vector<double> &Xs, vector<double> &Ys, vector<double> &Xf, vector<double> &Yf, int I, int J){
    return abs(Xf[I] - Xs[J]) + abs(Yf[I] - Ys[J]);
}

double CalcLatLong (vector<double> &Xs, vector<double> &Ys, vector<double> &Xf, vector<double> &Yf, int n, double *slatit, double* slongit, double *flatit, double* flongit){
    double PI = 3.141592, min;
    int deg;
    
    for (int i = 0; i < n; i++) {
        deg = (int) Xs[i];
        min = Xs[i] - deg;
        // slatit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
        slatit[i] = Xs[i] * PI / 180.0;
      
        deg = (int) Xf[i];
        min = Xf[i] - deg;
        // flatit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
        flatit[i] = Xf[i] * PI/ 180.0;
    }
    
    for (int i = 0; i < n; i++) {
        deg = (int) Ys[i];
        min = Ys[i] - deg;
        // slongit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
        slongit[i] = Ys[i] * PI / 180.0;

        deg = (int) Yf[i];
        min = Yf[i] - deg;
        // flongit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
        flongit[i] = Yf[i] * PI / 180;
    }
    return 0;
}


double CalcDistGeo (double *slatit, double* slongit, double *flatit, double* flongit, int I, int J){
    double q1, q2, q3, RRR = 6378.388;
    
    q1 = cos(flongit[I] - slongit[J]);
    q2 = cos(flatit[I] - slatit[J]);
    q3 = cos(flatit[I] + slatit[J]);

    return (RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3)));
    // (int) (RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3)));
}

string getInstanceType (char **argv){

    string filename(argv[1]);

    string::size_type loc = filename.find_first_of("/");
    string::size_type loc2 = filename.find_last_of("/", filename.size());
    string InstanceType;

    InstanceType.append(filename, loc+1, loc2-loc-1 );

    return InstanceType;
}

string getInstName (char **argv){

    string filename(argv[1]);

    string::size_type loc = filename.find_last_of("/");
    string::size_type loc2 = filename.find_first_of(".");
    string InstanceName;

    InstanceName.append(filename, loc+1, loc2-loc-1 );



    return InstanceName;
}

void getInstParam (instanceStat *inst, vector<int> &instParam){

    string::size_type loc = inst->InstName.find_first_of("-");
    string::size_type loc2 = inst->InstName.find_last_of("-");
    string param1;
    string param2;

    param1.append(inst->InstName, loc+1, loc2-loc-1);
    // param2.append(inst->InstName.c_str(), loc2+1);
    param2.append(inst->InstName, loc2+1, inst->InstName.back());

    // cout << "\nParam2: " << param2 << endl;
    // getchar();
 
    instParam.push_back(stoi(param1));
    instParam.push_back(stoi(param2));
}