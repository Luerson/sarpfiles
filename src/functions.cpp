#include "functions.h"
#include <cstdlib>
#include <stdio.h>

void readData (int argc, char** argv, nodeStat *node, instanceStat *inst, vector<nodeStat> &nodeVec, double ***Mdist, probStat* problem, vector<double> &passProfit){
    
    if (argc < 4) {
        cout << "\nMissing parameters\n";
        cout << " ./exeSARP [Instance] [Optimization strategy] [Scenario]"<< endl;
        exit(1);
    }
    
    if (argc > 4) {
        cout << "\nToo many parameters\n";
        cout << " ./exeSARP [Instance] [Optimization strategy] [Scenario]" << endl;
        exit(1);
    }  

    if (argv[2] == "sim"){
        problem->sim = true;
    }
    else if (argv[2] == "seq"){
        problem->seq = true;
    }

    problem->scen = argv[3];

    string file;
    int n;
    int m;
    int K;
    double service;
    double T;
    int V;
    int originalV;
    int dummy;

    string instType;

    char *instance; 
    instance = argv[1];

    ifstream in(instance, ios::in);

    if( !in ) {
        cout << "the file could not be opened\n";
        exit (1);
    }

    instType = getInstanceType(argv);

    if (instType == "myinstances"){

        in >> K;
        in >> service;
        in >> n;
        in >> m;
        in >> T;

        V = n + 2*m + K;
        // inst->dummy = K;
        inst->dummy = 1;

        double *xs = new double[V + inst->dummy];
        double *ys = new double[V + inst->dummy];
        char *label = new char[V + inst->dummy];
        double *load = new double[V + inst->dummy];
        double *e = new double[V + inst->dummy];
        double *l = new double[V + inst->dummy];
        double *xf = new double[V + inst->dummy];
        double *yf = new double[V + inst->dummy];
        double *delta = new double[V + inst->dummy];
        double *profit = new double[V+inst->dummy];

        double **dist = new double*[V + inst->dummy];
        for (int i= 0; i < V + inst->dummy; i++){
            dist[i] = new double [V + inst->dummy];
        }

        int tempNode;
        // for (int i = 0; i < V; i++){
        //     nodeVec.push_back(*node);
        // }

        for (int i = 0; i < V; i++){
            in >> tempNode >> xs[i] >> ys[i] >> label[i] >> load[i] >> e[i] >> l[i] >> xf[i] >> yf[i];
        }

        // Calculate distance matrix (Euclidian)
        double singleProfit;
        for (int i = 0; i < V + inst->dummy; i++){
            if (i < n){ 
               delta[i] = (2 * (service/60)) + (floor(calcEucDist(xs, ys, xf, yf, i, i) + 0.5))/inst->vmed;
               profit[i] = inst->gamma2 + inst->mu2*floor(calcEucDist(xs, ys, xf, yf, i, i) + 0.5) - floor(calcEucDist(xs, ys, xf, yf, i, i) + 0.5);
               // cout << "delta " << i << ": " << delta[i] << endl;
            }
            else if (i < V - K){ 
               delta[i] = service/60;
               profit[i] = 0;
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
                            dist[i][j] = floor(calcEucDist(xs, ys, xf, yf, i, j) + 0.5);
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
            node->xs = xs[i];
            node->ys = ys[i];
            node->label = label[i];
            node->load = load[i];
            node->e = e[i]/60;
            node->l = l[i]/60;
            node->xf = xf[i];
            node->yf = yf[i];
            node->delta = delta[i];
            node->profit = profit[i];
            nodeVec.push_back(*node);

        }

        //Adding dummy nodes
        for (int i = 0; i < inst->dummy; i++){
            node->xs = 0;
            node->ys = 0;
            node->label = 'f';
            node->load = 0;
            node->e = 0;
            node->l = 240/60;
            node->xf = 0;
            node->yf = 0;
            node->delta = 0;
            node->profit = 0;
            nodeVec.push_back(*node);
        }

        *Mdist = dist;
        inst->K = K;
        inst->n = n;
        inst->m = m;
        inst->T = T/60;
        inst->V = V;
        // cout << "\nNode vec: ";
        // for (int i = 0; i < nodeVec.size(); i++){
        //     cout << nodeVec[i].label << " ";
        // }
        // cout << endl;

        if( problem->scen == "1A" ){
            inst->nCluster = inst->n + inst->K + inst->dummy;
        }

        delete[] profit;
        delete[] xs;
        delete[] ys;
        delete[] label;
        delete[] load;
        delete[] e;
        delete[] l;
        delete[] xf;
        delete[] yf;
    }

    else if (instType == "sarpdata"){

        string filename(argv[1]);
        string InstanceName;

        string::size_type loc = filename.find_last_of("r", filename.size());
        string::size_type loc2 = filename.find_last_of(".", filename.size());
       
        InstanceName.append(filename, loc+1, loc2-loc-1);
        std::replace(InstanceName.begin(), InstanceName.end(), '_', ' ');
        vector<int> instanceData;
        stringstream ss(InstanceName);
        int temp;
        while (ss >> temp){
            instanceData.push_back(temp);
        }

        K = instanceData[0];
        n = instanceData[1]/2;
        m = instanceData[2]/2;
        service = 5;
        V = n + 2*m + K;
        originalV = 2*n + 2*m + 2*K;

        inst->dummy = 1;

        double *delta = new double[V + inst->dummy];
        double *slatitude = new double [V + inst->dummy];
        double *slongitude = new double [V + inst->dummy];
        double *flatitude = new double [V + inst->dummy];
        double *flongitude = new double [V + inst->dummy];
        double *profit = new double[V+inst->dummy];


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

        for (int i = 0; i < originalV; i++){
            vxs.push_back(0);
            vys.push_back(0);
            vload.push_back(0);
            ve.push_back(0);
        }

        for (int i = 0; i < originalV; i++){
            in >> vxs[i] >> vys[i] >> vload[i] >> ve[i];
        }

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

        for (int i = 0; i < n; i++){
            vxs.erase(vxs.begin() + n + m);
            vys.erase(vys.begin() + n + m);
            vload.erase(vload.begin() + n + m);
            ve.erase(ve.begin() + n + m);
            vxf.erase(vxf.begin() + n + m);
            vyf.erase(vyf.begin() + n + m);
        }

        for (int i = 0; i < ve.size(); i++){
            vl.push_back(ve[i]);
        }

        for (int i = n; i < vl.size(); i++){
            vl[i] = 14*60;
        }

        // cout << "vl: " << endl;
        // for(int i = 0; i < vl.size(); i++){
        //     cout << vl[i] << " ";
        // }

        CalcLatLong ( vxs, vys, vxf, vyf, V, slatitude, slongitude, flatitude, flongitude );
        
        for (int i = 0; i < V + inst->dummy; i++){
            if (i < n){ 
               delta[i] = (2 * (service/60)) + (CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i))/inst->vmed;
               profit[i] = inst->gamma2 + inst->mu2*CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i) - CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i);

               // cout << "delta " << i << ": " << delta[i] << endl;
                // delta[i] = (2 * (service/60)) + (CalcMan(vxs, vys, vxf, vys, i, i))/inst->vmed;

            }
            else if (i < V - K){ 
               delta[i] = service/60;
               profit[i] = 0;
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
                            dist[i][j] = CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i);
                            // dist[i][j] = CalcMan(vxs, vys, vxf, vys, i, i);
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


        *Mdist = dist;
        inst->K = K;
        inst->n = n;
        inst->m = m;
        inst->V = V;
        
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
            nodeVec.push_back(*node);
        }

        for (int i = 0; i < inst->dummy; i++){
            node->xs = 0;
            node->ys = 0;
            node->load = 0;
            node->e = 0;
            node->l = 14*60;
            node->xf = 0;
            node->yf = 0;
            node->delta = 0;
            node->profit = 0;
            nodeVec.push_back(*node);
        }

        if( problem->scen == "1A" ){
            inst->nCluster = inst->n + inst->K + inst->dummy;
        }
        
        delete[] profit;
        delete[] delta;
        delete[] slatitude;
        delete[] slongitude;
        delete[] flatitude;
        delete[] flongitude;
    }
 
    else if (instType == "sf_data"){

        in >> K;
        in >> service;
        in >> n;
        in >> m;

        V = n + 2*m + K;

        originalV = 2*n + 2*m + 2; 

        inst->dummy = 1;

        double *delta = new double[V + inst->dummy];
        double *slatitude = new double [V + inst->dummy];
        double *slongitude = new double [V + inst->dummy];
        double *flatitude = new double [V + inst->dummy];
        double *flongitude = new double [V + inst->dummy];
        double *profit = new double[V+inst->dummy];

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

        for (int i = 0; i < n; i++){
            vxs.erase(vxs.begin() + n + m);
            vys.erase(vys.begin() + n + m);
            vload.erase(vload.begin() + n + m);
            ve.erase(ve.begin() + n + m);
            vxf.erase(vxf.begin() + n + m);
            vyf.erase(vyf.begin() + n + m);
        }

        for (int i = 0; i < n; i++){
            vl[i] = ve[i];
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


        // cout << "Testing: " << endl;
        // for (int i = 0; i < vxs.size(); i++){
        //     cout << "i: " << i << " vxs: " << vxs[i] << " vys: " << vys[i] << " vload: " << vload[i] << " ve: " << ve[i] << " vl: " << vl[i] << " vxf: " << vxf[i] << " vyf: " << vyf[i] << endl; 
        // }
        // getchar();


        // Calculate distance matrix (Geolocation)

        CalcLatLong ( vxs, vys, vxf, vyf, V, slatitude, slongitude, flatitude, flongitude );
        
        double singleProfit;
        for (int i = 0; i < V + inst->dummy; i++){
            if (i < n){ 
                delta[i] = (2 * (service/60)) + (CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i))/inst->vmed;
                profit[i] = inst->gamma2 + inst->mu2*CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i) - CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i);
               // cout << "delta " << i << ": " << delta[i] << endl;
                // delta[i] = (2 * (service/60)) + (CalcMan(vxs, vys, vxf, vys, i, i))/inst->vmed;

            }
            else if (i < V - K){ 
               delta[i] = service/60;
               profit[i] = 0;
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
                            dist[i][j] = CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i);
                            // dist[i][j] = CalcMan(vxs, vys, vxf, vys, i, i);
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


        *Mdist = dist;
        inst->K = K;
        inst->n = n;
        inst->m = m;
        inst->V = V;

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
            nodeVec.push_back(*node);
        }

        //Adding dummy nodes
        for (int i = 0; i < inst->dummy; i++){
            node->xs = 0;
            node->ys = 0;
            node->load = 0;
            node->e = 0;
            node->l = 14*60;
            node->xf = 0;
            node->yf = 0;
            node->delta = 0;
            node->profit = 0;
            nodeVec.push_back(*node);
        }

        if( problem->scen == "1A" ){
            inst->nCluster = inst->n + inst->K + inst->dummy;
        }

        delete[] profit;
        delete[] delta;
        delete[] slatitude;
        delete[] slongitude;
        delete[] flatitude;
        delete[] flongitude;
   
    }

}

double calcEucDist (double *Xs, double *Ys, double *Xf, double *Yf, int I, int J){
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

    return
    (int) (RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3)) + 1.0);
    // (int) (RRR * acos( 0.5*((1.0+q1)*q2 - (1.0-q1)*q3)));
}

// void feasibleArcs (instanceStat *inst, vector<nodeStat> &nodeVec, vector< vector<bool> > &arcs, pair<int, int> &fArc, vector< vector< pair<int,int> > > &arcPlus, vector< vector< pair<int,int> > > &arcMinus){

//     for (int i = 0; i < inst->V; i++){
//         if (i < inst->n){
//             for(int j = 0; j < inst->n + 2*inst->m; j++){
//                 if(i != j){
//                     arcs[i][j] = true;
//                     fArc.first = i;
//                     fArc.second = j;
//                     arcMinus[j].push_back(fArc);
//                     arcPlus[i].push_back(fArc);
//                 }
//             }
//             arcs[i][inst->V] = true;
//             fArc.first = i;
//             fArc.second = inst->V;
//             arcMinus[inst->V].push_back(fArc);
//             arcPlus[i].push_back(fArc);
//         }

//         else if (i > inst->n - 1){
//             if (i < inst->n + inst->m){
//                 for (int j = 0; j < inst->n + 2*inst->m; j++){
//                     if (i != j){
//                         if (j != i + inst->m){
//                             arcs[i][j] = true;
//                             fArc.first = i;
//                             fArc.second = j;
//                             arcMinus[j].push_back(fArc);
//                             arcPlus[i].push_back(fArc);
//                         }
//                     }
//                 }
//             }
//             else if (i > inst->n + inst->m - 1){
//                 if (i < inst->n + 2*inst->m){
//                     for (int j = 0; j < inst->n + 2*inst->m; j++){
//                         if (i != j){
//                             if (j + inst->m != i){
//                                 arcs[i][j] = true;
//                                 fArc.first = i;
//                                 fArc.second = j;
//                                 arcMinus[j].push_back(fArc);
//                                 arcPlus[i].push_back(fArc);
//                             }
//                         }
//                     }
//                     arcs[i][inst->V] = true;
//                     fArc.first = i;
//                     fArc.second = inst->V;
//                     arcMinus[inst->V].push_back(fArc);
//                     arcPlus[i].push_back(fArc);                   
//                 }

//                 else if (i > inst->n + 2* inst->m - 1){
//                     for (int j = 0; j < inst->n + inst->m; j++){
//                         arcs[i][j] = true;
//                         fArc.first = i;
//                         fArc.second = j;
//                         arcMinus[j].push_back(fArc);
//                         arcPlus[i].push_back(fArc);
//                     } 
//                     arcs[i][inst->V] = true;
//                     fArc.first = i;
//                     fArc.second = inst->V;
//                     arcMinus[inst->V].push_back(fArc);
//                     arcPlus[i].push_back(fArc);             
//                 }
//             }
//         }
//     }
// }

void feasibleArcs (instanceStat *inst, vector<nodeStat> &nodeVec, vector< vector<bool> > &arcs, pair<int, int> &fArc, vector< vector< pair<int,int> > > &arcPlus, vector< vector< pair<int,int> > > &arcMinus, probStat* problem, vector< pair<int,int> > &arcNN){

    if (problem->scen == "1A" || problem->scen == "2A"){
           for (int i = 0; i < inst->V; i++){
            if (i < inst->n){
                for(int j = 0; j < inst->n + 2*inst->m; j++){
                    if(i != j){
                        arcs[i][j] = true;
                        fArc.first = i;
                        fArc.second = j;
                        arcMinus[j].push_back(fArc);
                        arcPlus[i].push_back(fArc);
                        if (j < inst->n){
                            arcNN.push_back(fArc);
                        }
                    }
                }
                for (int j = inst->V; j < inst->V + inst->dummy; j++){
                    arcs[i][j] = true;
                    fArc.first = i;
                    fArc.second = j;
                    arcMinus[j].push_back(fArc);
                    arcPlus[i].push_back(fArc);
                }
            }

            else if (i > inst->n - 1){
                if (i < inst->n + inst->m){
                    for (int j = 0; j < inst->n + 2*inst->m; j++){
                        if (i != j){
                            if (j != i + inst->m){
                                arcs[i][j] = true;
                                fArc.first = i;
                                fArc.second = j;
                                arcMinus[j].push_back(fArc);
                                arcPlus[i].push_back(fArc);
                            }
                        }
                    }
                }
                else if (i > inst->n + inst->m - 1){
                    if (i < inst->n + 2*inst->m){
                        for (int j = 0; j < inst->n + 2*inst->m; j++){
                            if (i != j){
                                if (j + inst->m != i){
                                    arcs[i][j] = true;
                                    fArc.first = i;
                                    fArc.second = j;
                                    arcMinus[j].push_back(fArc);
                                    arcPlus[i].push_back(fArc);
                                }
                            }
                        }
                        for (int j = inst->V; j < inst->V + inst->dummy; j++){
                            arcs[i][j] = true;
                            fArc.first = i;
                            fArc.second = j;
                            arcMinus[j].push_back(fArc);
                            arcPlus[i].push_back(fArc);
                        }
                    }

                    else if (i > inst->n + 2* inst->m - 1){
                        for (int j = 0; j < inst->n + inst->m; j++){
                            arcs[i][j] = true;
                            fArc.first = i;
                            fArc.second = j;
                            arcMinus[j].push_back(fArc);
                            arcPlus[i].push_back(fArc);
                        }
                        // for (int j = inst->V; j < inst->V + inst->dummy; j++){
                        //     arcs[i][j] = true;
                        //     fArc.first = i;
                        //     fArc.second = j;
                        //     arcMinus[j].push_back(fArc);
                        //     arcPlus[i].push_back(fArc);
                        // }
                    }
                }
            }
        }
    }
    if (problem->scen == "1B" || problem->scen == "2B"){
        for (int i = 0; i < inst->V; i++){
            if (i < inst->n){
                for(int j = 0; j < inst->n + 2*inst->m; j++){
                    if(i != j){
                        arcs[i][j] = true;
                        fArc.first = i;
                        fArc.second = j;
                        arcMinus[j].push_back(fArc);
                        arcPlus[i].push_back(fArc);
                    }
                }
                for (int j = inst->V; j < inst->V + inst->dummy; j++){
                    arcs[i][j] = true;
                    fArc.first = i;
                    fArc.second = j;
                    arcMinus[j].push_back(fArc);
                    arcPlus[i].push_back(fArc);
                }
            }

            else if (i > inst->n - 1){
                if (i < inst->n + inst->m){
                    for (int j = 0; j < inst->n + inst->m; j++){
                        if (i != j){
                            if (j != i + inst->m){
                                arcs[i][j] = true;
                                fArc.first = i;
                                fArc.second = j;
                                arcMinus[j].push_back(fArc);
                                arcPlus[i].push_back(fArc);
                            }
                        }
                    }
                }
                else if (i > inst->n + inst->m - 1){
                    if (i < inst->n + 2*inst->m){
                        for (int j = 0; j < inst->n + 2*inst->m; j++){
                            if (i != j){
                                if (j + inst->m != i){
                                    arcs[i][j] = true;
                                    fArc.first = i;
                                    fArc.second = j;
                                    arcMinus[j].push_back(fArc);
                                    arcPlus[i].push_back(fArc);
                                }
                            }
                        }
                        for (int j = inst->V; j < inst->V + inst->dummy; j++){
                            arcs[i][j] = true;
                            fArc.first = i;
                            fArc.second = j;
                            arcMinus[j].push_back(fArc);
                            arcPlus[i].push_back(fArc);
                        }
                    }

                    else if (i > inst->n + 2* inst->m - 1){
                        for (int j = 0; j < inst->n + inst->m; j++){
                            arcs[i][j] = true;
                            fArc.first = i;
                            fArc.second = j;
                            arcMinus[j].push_back(fArc);
                            arcPlus[i].push_back(fArc);
                        }
                        // for (int j = inst->V; j < inst->V + inst->dummy; j++){
                        //     arcs[i][j] = true;
                        //     fArc.first = i;
                        //     fArc.second = j;
                        //     arcMinus[j].push_back(fArc);
                        //     arcPlus[i].push_back(fArc);
                        // }
                    }
                }
            }
        }
    }
}

string getInstanceType (char **argv){

    string filename(argv[1]);

    string::size_type loc = filename.find_first_of("/");
    string::size_type loc2 = filename.find_last_of("/", filename.size());
    string InstanceType;

    InstanceType.append(filename, loc+1, loc2-loc-1 );

    return InstanceType;
}

void makeBundles (instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat, vector<int> &clusters, vector< vector<int> > &clusterVec){

//For 1A:
    int counter = 0;
    for (int i = 0; i < inst->n; i++){
        bStat->bundle.push_back(i);
        bStat->bundleVec.push_back(bStat->bundle);
        // clusters.push_back(bStat->bundle);
        clusters.push_back(bStat->bundleVec.size()-1);
        bStat->bundle.clear();
        for (int j = inst->n; j < inst->m + inst->n; j++){
            bStat->bundle.push_back(j);
            bStat->bundle.push_back(i);
            bStat->bundle.push_back(j+inst->m);
            bStat->bundleTimes.push_back(nodeVec[j].delta + nodeVec[i].delta + nodeVec[j+inst->m].delta);

            // cout << "\nBundle: [";
            // for (int k = 0; k < bStat->bundle.size(); k++){
            //     cout << bStat->bundle[k];
            //     if (k < bStat->bundle.size() - 1){
            //         cout << ", ";
            //     }
            // }
            // cout << "]" << endl;

            bStat->bundleVec.push_back(bStat->bundle);
            // clusters.push_back(bStat->bundle);
            clusters.push_back(bStat->bundleVec.size()-1);
            bStat->bundle.clear();
        }
        clusterVec.push_back(clusters);
        clusters.clear();
    }

    for (int i = 2*inst->m + inst->n; i < nodeVec.size(); i++){
        bStat->bundle.push_back(i);
        bStat->bundleVec.push_back(bStat->bundle);
        // clusters.push_back(bStat->bundle);
        clusters.push_back(bStat->bundleVec.size()-1);
        bStat->bundle.clear();
        
        clusterVec.push_back(clusters);
        clusters.clear();
    }


}

void bundleProfit(instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat, vector<double> &passProfit){
    double cost;

    // for (int i = 0; i < bStat->bundleVec.size(); i++){
    //     if (bStat->bundleVec[i].size() <= 1){
    //         cost = 0;
    //         bStat->bundleProfVec.push_back(cost);
    //     }
    //     else{
    //         for (int j = 0; j < bStat->bundleVec[i].size() - 1; j++){
    //             cost = 0;
    //             if (bStat->bundleVec[i][j] > inst->n - 1 && bStat->bundleVec[i][j] < inst->n + inst->m){
    //                 cost += inst->gamma + inst->mu*mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j] + inst->m];
    //                 for (int k = j + 1; k < bStat->bundleVec[i].size(); k++){
    //                     if (bStat->bundleVec[i][k] == bStat->bundleVec[i][j] + inst->m){
    //                         break;
    //                     }
    //                     else if (bStat->bundleVec[i][k] < inst->n + inst->m){
    //                         cost += - mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][k]] - mdist[bStat->bundleVec[i][k]][bStat->bundleVec[i][j] + inst->m];
                            
    //                     }
    //                 }
    //                 bStat->bundleProfVec.push_back(cost);
    //             }
    //         }            
    //     }
    // }

    for (int i = 0; i < bStat->bundleVec.size(); i++){
        if (bStat->bundleVec[i].size() <= 1){
            cost = nodeVec[bStat->bundleVec[i][0]].profit;
            bStat->bundleProfVec.push_back(cost);
        }
        else{
            for (int j = 0; j < bStat->bundleVec[i].size() - 1; j++){
                cost = 0;
                if (bStat->bundleVec[i][j] > inst->n - 1 && bStat->bundleVec[i][j] < inst->n + inst->m){
                    cost += inst->gamma + inst->mu*mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j] + inst->m];
                    cost += - mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]];
                    bStat->bundleProfVec.push_back(cost);
                    continue;
                }
                else if (bStat->bundleVec[i][j] < inst->n){
                    cost += nodeVec[bStat->bundleVec[i][j]].profit;
                    cost += - mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]];
                    bStat->bundleProfVec.push_back(cost);
                    continue;
                } 
                else if (bStat->bundleVec[i][j] > inst->n + inst->m - 1){
                    cost += nodeVec[bStat->bundleVec[i][j]].profit;
                    cost += - mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]];
                    bStat->bundleProfVec.push_back(cost);
                    continue;
                }               
            }
        }
    } 
}

void feasibleBundleArcs (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat){
    //1-A
    int setN = bStat->bundleVec.size() - inst->K - 1;
    int currentCluster = 0;

    for(int i = 0; i < bStat->bundleVec.size() - 1; i++){
        if (i > currentCluster*(inst->m + 1) + inst->m){
            currentCluster++;
        }
        if(i < setN){
            for (int j = 0; j < setN; j++){
                if (i != j){
                    if (j > currentCluster*(inst->m + 1) + inst->m || j < currentCluster*(inst->m + 1)){
                        if (bStat->bundleStart[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                            bStat->bArcs[i][j] = true;
                            // bStat->bFArc.first = i;
                            // bStat->bFArc.second = j;
                            // bStat->bArcMinus[j].push_back(bStat->bFArc);
                            // bStat->bArcPlus[i].push_back(bStat->bFArc);                            
                        }
                    }
                } 
            }
            bStat->bArcs[i][bStat->bundleVec.size()-1] = true;
            // bStat->bFArc.first = i;
            // bStat->bFArc.second = bStat->bundleVec.size()-1;
            // bStat->bArcMinus[bStat->bundleVec.size()-1].push_back(bStat->bFArc);
            // bStat->bArcPlus[i].push_back(bStat->bFArc);
        }
        else if (i >= setN){
            for (int j = 0; j < setN; j++){
                if (bStat->bundleStart[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                    bStat->bArcs[i][j] = true;
                    // bStat->bFArc.first = i;
                    // bStat->bFArc.second = j;
                    // bStat->bArcMinus[j].push_back(bStat->bFArc);
                    // bStat->bArcPlus[i].push_back(bStat->bFArc);
                }
            }
        }
    }

    //remove arcs to bundles with same parcel requests
    for (int i = 0; i < bStat->parcelBundleVec.size(); i++){
        for (int j = 0; j < bStat->parcelBundleVec[i].size(); j++){
            for (int k = 0; k < bStat->parcelBundleVec[i].size(); k++){
                if (bStat->parcelBundleVec[i][j] != bStat->parcelBundleVec[i][k]){
                    bStat->bArcs[bStat->parcelBundleVec[i][j]][bStat->parcelBundleVec[i][k]] = false;
                }
            }
        }
    }

    //remove arcs from/to bundles with negative start times
    for (int i = 0; i < bStat->bundleStart.size(); i++){
        if (bStat->bundleStart[i] < 0){
            for (int j = 0; j < bStat->bundleVec.size(); j++){
                bStat->bArcs[j][i] = false;
                bStat->bArcs[i][j] = false;
            }
        }
    }



    for (int i = 0; i < bStat->bundleVec.size(); i++){
        for (int j = 0; j < bStat->bundleVec.size(); j++){
            if (bStat->bArcs[i][j]){
                bStat->bFArc.first = i;
                bStat->bFArc.second = j;
                bStat->bArcMinus[j].push_back(bStat->bFArc);
                bStat->bArcPlus[i].push_back(bStat->bFArc); 
            }
        }
    }

}

void feasibleClusterArcs (instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat, vector< vector<int> > &clusterVec, pair<int, int> &cFArc, vector< vector<bool> > &cArcs, vector< vector< pair<int,int> > > &cArcPlus, vector< vector< pair<int,int> > > &cArcMinus){
    
    int reqClusters = clusterVec.size() - inst->K - 1;
    
    // for(int i = 0; i < clusterVec.size() - 1; i++){
    //     if(i < reqClusters){
    //         for (int j = 0; j < reqClusters; j++){
    //             if (i != j){
    //                 cArcs[i][j] = true;
    //                 cFArc.first = i;
    //                 cFArc.second = j;
    //                 cArcMinus[j].push_back(cFArc);
    //                 cArcPlus[i].push_back(cFArc);
    //             } 
    //         }
    //         cFArc.first = i;
    //         cFArc.second = clusterVec.size()-1;
    //         cArcs[i][clusterVec.size()-1] = true;
    //         cArcMinus[clusterVec.size()-1].push_back(cFArc);
    //         cArcPlus[i].push_back(cFArc);
    //     }
    //     else if (i >= reqClusters){
    //         for (int j = 0; j < reqClusters; j++){
    //             cArcs[i][j] = true;
    //             cFArc.first = i;
    //             cFArc.second = j;
    //             cArcMinus[j].push_back(cFArc);
    //             cArcPlus[i].push_back(cFArc);
    //         }
    //     }
    // }
    
    //1-A
    
    // for (int i = 0; i < reqClusters; i++){
    //     for (int j = 0; j < reqClusters; j++){
    //         if (i != j){
    //             if ()
    //         }
    //     }
    // }
    int clusterA = 0;
    int clusterB;
    int setN = bStat->bundleVec.size() - inst->K - 1;

    for(int i = 0; i < bStat->bundleVec.size(); i++){
        if (i < setN){
            if (i > clusterA*(inst->m + 1) + inst->m){
                clusterA++;
            }           
        }
        else{
            clusterA++;
        }
        clusterB = 0;
        for (int j = 0; j < bStat->bundleVec.size(); j++){
            if (j < setN){
                if (j > clusterB*(inst->m + 1) + inst->m){
                    clusterB++;
                }                
            }
            else{
                clusterB++;
            }
            if (clusterA == clusterB){
                continue;
            }
            else{
                if (bStat->bArcs[i][j]){
                    cArcs[clusterA][clusterB] = true;
                    cFArc.first = i;
                    cFArc.second = j;
                    cArcMinus[clusterB].push_back(cFArc);
                    cArcPlus[clusterA].push_back(cFArc);
                }
            }
        }
    }
}

void makeParcelBundles(instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat){
    int parcelReq;
    for (int i = 0; i < bStat->bundleVec.size(); i++){
        for (int j = 0; j < bStat->bundleVec[i].size(); j++){
            if (bStat->bundleVec[i][j] < inst->n){
                break;
            }
            else if (bStat->bundleVec[i][j] > inst->n + inst->m - 1){
                break;
            }
            else{
                parcelReq = bStat->bundleVec[i][j];
                bStat->parcelBundleVec[parcelReq - inst->n].push_back(i);
            }

        }
    }

    // cout << "\nParcel bundles: ";
    
    // for (int i = 0; i < parcelBundleVec.size(); i++){
    //     cout << "\np req: " << i + inst->n << " - bundles: ";
    //     for(int j = 0; j < parcelBundleVec[i].size(); j++){         
    //         if (j == parcelBundleVec[i].size()-1){
    //             cout << parcelBundleVec[i][j];
    //         }
    //         else{
    //             cout << parcelBundleVec[i][j] << ", ";
    //         }
    //     }
    // }
    // getchar();

}

void makeStartTimes (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat){

    double parcelTime;
    double bundleTime;
    bool firstPassenger;

    for (int i = 0; i < bStat->bundleVec.size(); i++){
        parcelTime = 0;
        bundleTime = 0;
        firstPassenger = false;
        for (int j = 0; j < bStat->bundleVec[i].size(); j++){
            if (bStat->bundleVec[i].size() > 1){
                if (bStat->bundleVec[i][j] >= inst->n){
                    parcelTime += ((mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]])/inst->vmed) + nodeVec[bStat->bundleVec[i][j]].delta;

                    if (bStat->bundleVec[i][j + 1] < inst->n){
                        bundleTime = nodeVec[bStat->bundleVec[i][j + 1]].e - parcelTime;
                        firstPassenger = true;
                        break;
                    }
                }
            }
            else{
                firstPassenger = true;
                bundleTime = nodeVec[bStat->bundleVec[i][j]].e;
                break;
            }
        }
        if (firstPassenger){
            bStat->bundleStart.push_back(bundleTime);
        }
    }
}

void makeBundleReference(instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat){

    for (int i = 0; i < bStat->bundleVec.size(); i++){

        bStat->lastElement.push_back(bStat->bundleVec[i][bStat->bundleVec[i].size()-1]);

    }
    // cout << "\nLast Element Vector: " << endl;

    // for (int i = 0; i < bStat->lastElement.size(); i++){
    //     cout << bStat->lastElement[i] << ", ";
    // }

    for (int i = 0; i < bStat->bundleVec.size(); i++){

        bStat->firstElement.push_back(bStat->bundleVec[i][0]);

    }
    // cout << "\nFirst Element Vector: " << endl;

    // for (int i = 0; i < bStat->firstElement.size(); i++){
    //     cout << bStat->firstElement[i] << ", ";
    // }
}
