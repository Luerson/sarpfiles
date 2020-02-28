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
               // cout << "profit " << i << ": " << profit[i] << endl;
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

        ve[ve.size()-1] = 0;

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
                            dist[i][j] = CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, j);
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

        ve[ve.size()-1] = 0;

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
        // Calculate distance matrix (Geolocation)

        CalcLatLong ( vxs, vys, vxf, vyf, V, slatitude, slongitude, flatitude, flongitude );
        
        double singleProfit;
        for (int i = 0; i < V + inst->dummy; i++){
            if (i < n){ 
                delta[i] = (2 * (service/60)) + (CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i))/inst->vmed;
                profit[i] = inst->gamma2 + inst->mu2*CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i) - CalcDistGeo(slatitude, slongitude, flatitude, flongitude, i, i);
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

        delete[] profit;
        delete[] delta;
        delete[] slatitude;
        delete[] slongitude;
        delete[] flatitude;
        delete[] flongitude;
   
    }

    if(problem->scen == "1A" || "1B"){
        inst->nCluster = inst->n + inst->K + inst->dummy;
    }
    else if(problem->scen == "2A" || "2B"){
        inst->nCluster = inst->n + inst->K + inst->dummy;
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

void makeBundles (instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat, vector<int> &clusters, vector< vector<int> > &clusterVec, vector< vector<int> > &clsParcel, probStat* problem){
//*******
//For 1A:
//*******    
    if (problem->scen == "1A"){
        int counter = 0;
        for (int i = 0; i < inst->n; i++){
            bStat->bundle.push_back(i);
            bStat->bundleVec.push_back(bStat->bundle);
            // clusters.push_back(bStat->bundle);
            clusters.push_back(bStat->bundleVec.size()-1);
            bStat->bundle.clear();
            for (int j = 0; j < clsParcel[i].size(); j++){
                bStat->bundle.push_back(clsParcel[i][j]);
                bStat->bundle.push_back(i);
                bStat->bundle.push_back(clsParcel[i][j]+inst->m);
                bStat->bundleTimes.push_back(nodeVec[clsParcel[i][j]].delta + nodeVec[i].delta + nodeVec[clsParcel[i][j]+inst->m].delta);
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


//*******
//For 2A:
//******* 

    if (problem->scen == "2A"){
        int counter = 0;
        for (int i = 0; i < inst->n; i++){
            bStat->bundle.push_back(i);
            bStat->label.push_back(0);
            bStat->mainNode.push_back(i);
            bStat->bundleVec.push_back(bStat->bundle);
            // clusters.push_back(bStat->bundle);
            clusters.push_back(bStat->bundleVec.size()-1);
            bStat->bundle.clear();
            for (int j = 0; j < clsParcel[i].size(); j++){
                if (clsParcel[i][j] > inst->n + 2*inst->m){
                    bStat->bundle.push_back(clsParcel[i][j]/inst->m);
                    bStat->bundle.push_back(i);
                    bStat->label.push_back(2);
                    bStat->mainNode.push_back(clsParcel[i][j]/inst->m);
                    bStat->bundleTimes.push_back(nodeVec[clsParcel[i][j]].delta + nodeVec[i].delta);
                    bStat->bundleVec.push_back(bStat->bundle);
                    clusters.push_back(bStat->bundleVec.size()-1);
                    bStat->activeDL.push_back(bStat->bundleVec.size()-1);
                    bStat->bundle.clear();
                } 
                
                else if (clsParcel[i][j] < inst->n + inst->m){
                    bStat->bundle.push_back(clsParcel[i][j]);
                    bStat->bundle.push_back(i);
                    bStat->label.push_back(1);
                    bStat->mainNode.push_back(clsParcel[i][j]);
                    bStat->bundleTimes.push_back(nodeVec[clsParcel[i][j]].delta + nodeVec[i].delta);
                    bStat->bundleVec.push_back(bStat->bundle);
                    clusters.push_back(bStat->bundleVec.size()-1);
                    bStat->activePU.push_back(bStat->bundleVec.size()-1);
                    bStat->bundle.clear();
                }

                else{
                    bStat->bundle.push_back(i);
                    bStat->bundle.push_back(clsParcel[i][j]);
                    bStat->mainNode.push_back(clsParcel[i][j]);
                    bStat->label.push_back(2);
                    bStat->bundleTimes.push_back(nodeVec[clsParcel[i][j]].delta + nodeVec[i].delta);
                    bStat->bundleVec.push_back(bStat->bundle);
                    clusters.push_back(bStat->bundleVec.size()-1);
                    bStat->activeDL.push_back(bStat->bundleVec.size()-1);
                    bStat->bundle.clear(); 
                }               
            }
            clusterVec.push_back(clusters);
            clusters.clear();
        }

        for (int i = 2*inst->m + inst->n; i < nodeVec.size(); i++){
            bStat->bundle.push_back(i);
            bStat->bundleVec.push_back(bStat->bundle);
            bStat->label.push_back(3);
            bStat->mainNode.push_back(i);
            // clusters.push_back(bStat->bundle);
            clusters.push_back(bStat->bundleVec.size()-1);
            bStat->bundle.clear();
            
            clusterVec.push_back(clusters);
            clusters.clear();
        }
    }

//*******
//For 1B:
//*******
    // else if (problem->scen == "1B"){
    //     int counter = 0;
    //     for (int i = 0; i < inst->n; i++){
    //         bStat->bundle.push_back(i);
    //         bStat->bundleVec.push_back(bStat->bundle);
    //         // clusters.push_back(bStat->bundle);
    //         clusters.push_back(bStat->bundleVec.size()-1);
    //         bStat->bundle.clear();
    //         for (int j = 0; j < clsParcel[i].size(); j++){
    //             bStat->bundle.push_back(clsParcel[i][j]);
    //             bStat->bundle.push_back(i);
    //             bStat->bundle.push_back(clsParcel[i][j]+inst->m);
    //             bStat->bundleTimes.push_back(nodeVec[clsParcel[i][j]].delta + nodeVec[i].delta + nodeVec[clsParcel[i][j]+inst->m].delta);
    //             bStat->bundleVec.push_back(bStat->bundle);
    //             // clusters.push_back(bStat->bundle);
    //             clusters.push_back(bStat->bundleVec.size()-1);
    //             bStat->bundle.clear();
    //         }
    //         clusterVec.push_back(clusters);
    //         clusters.clear();
    //     }

    //     for (int i = 2*inst->m + inst->n; i < nodeVec.size(); i++){
    //         bStat->bundle.push_back(i);
    //         bStat->bundleVec.push_back(bStat->bundle);
    //         // clusters.push_back(bStat->bundle);
    //         clusters.push_back(bStat->bundleVec.size()-1);
    //         bStat->bundle.clear();
            
    //         clusterVec.push_back(clusters);
    //         clusters.clear();
    //     }
    // }


}

void bundleProfit(instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat, vector<double> &passProfit){
    double cost;
    double service;

    for (int i = 0; i < bStat->bundleVec.size(); i++){
        cost = 0;
        service = 0;
        if (bStat->bundleVec[i].size() <= 1){
            cost = nodeVec[bStat->bundleVec[i][0]].profit;
            bStat->bundleProfVec.push_back(cost);
            service = nodeVec[bStat->bundleVec[i][0]].delta;
            bStat->bundleServVec.push_back(service);
        }
        else{
            for (int j = 0; j < bStat->bundleVec[i].size() - 1; j++){
                if (bStat->bundleVec[i][j] > inst->n - 1 && bStat->bundleVec[i][j] < inst->n + inst->m){ //parcel pickup
                    cost += inst->gamma + inst->mu*mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j] + inst->m];
                    cost += - mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]];
                    service += nodeVec[bStat->bundleVec[i][j]].delta + (mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]]/inst->vmed);
                    continue;
                }
                else if (bStat->bundleVec[i][j] < inst->n){//passenger request
                    cost += nodeVec[bStat->bundleVec[i][j]].profit;
                    cost += - mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]];
                    service += nodeVec[bStat->bundleVec[i][j]].delta + (mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]]/inst->vmed);
                    continue;
                } 
                else if (bStat->bundleVec[i][j] > inst->n + inst->m - 1){//parcel delivery
                    cost += nodeVec[bStat->bundleVec[i][j]].profit;
                    cost += - mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]];
                    service += nodeVec[bStat->bundleVec[i][j]].delta + (mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]]/inst->vmed);
                    continue;
                }               
            }
            cost += nodeVec[bStat->bundleVec[i][bStat->bundleVec[i].size() - 1]].profit;
            service += nodeVec[bStat->bundleVec[i][bStat->bundleVec[i].size() - 1]].delta;
            // cout << bStat->bundleVec[i][bStat->bundleVec[i].size() - 1] << " - delta: " << nodeVec[bStat->bundleVec[i][bStat->bundleVec[i].size() - 1]].delta;
            // getchar();
            // cout << "bundle: " << i << " - service(final): " << service << endl;
            // getchar(); 
            bStat->bundleProfVec.push_back(cost);
            bStat->bundleServVec.push_back(service);
        }
    } 
}

void feasibleBundleArcs (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat, int p, probStat* problem){
    int setN = bStat->bundleVec.size() - inst->K - 1;
    int currentCluster = 0;
    int ref;

    if (problem->scen == "1A"){
        if (p < 0){
            ref = inst->m;
        }
        else{
            ref = p;
        }

        for(int i = 0; i < bStat->bundleVec.size() - 1; i++){
            if (i > currentCluster*(ref + 1) + ref){
                currentCluster++;
            }
            if(i < setN){
                for (int j = 0; j < setN; j++){
                    if (i != j){
                        if (j > currentCluster*(ref + 1) + ref || j < currentCluster*(ref + 1)){
                            if (bStat->bundleStart[i] + bStat->bundleServVec[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                                bStat->bArcs[i][j] = true;
                            }
                        }
                    } 
                }
                bStat->bArcs[i][bStat->bundleVec.size()-1] = true;
            }
            else if (i >= setN){
                for (int j = 0; j < setN; j++){
                    if (bStat->bundleStart[i] + bStat->bundleServVec[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                        bStat->bArcs[i][j] = true;
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
        for (int i = 0; i < bStat->bundleStart.size() - 1; i++){
            if (bStat->bundleStart[i] < 0){
                for (int j = 0; j < bStat->bundleVec.size() - 1; j++){
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

    else if (problem->scen == "2A"){
        if (p < 0){
            ref = 3*inst->m;
        }
        else{
            ref = p;
        }

        for(int i = 0; i < bStat->bundleVec.size() - 1; i++){
            if (i > currentCluster*(ref + 1) + ref){
                currentCluster++;
            }
            if(i < setN){
                for (int j = 0; j < setN; j++){
                    if (i != j){
                        if (j > currentCluster*(ref + 1) + ref || j < currentCluster*(ref + 1)){
                            if (bStat->label[i] == 1){
                                if (bStat->label[j] == 0){
                                    if (bStat->bundleStart[i] + bStat->bundleServVec[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                                        bStat->bArcs[i][j] = true;
                                    }                                     
                                }
                                // else if (bStat->label[j] == 1){
                                //     if(bStat->mainNode[i] != bStat->mainNode[j]){  
                                //         bStat->bArcs[i][j] = true;
                                //     } 
                                // }
                                else if (bStat->label[j] == 2){
                                    if (bStat->mainNode[j] == bStat->mainNode[i] + inst->m){
                                        if (bStat->bundleStart[i] + bStat->bundleServVec[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                                            bStat->bArcs[i][j] = true;
                                        }
                                    }  
                                }
                            }
                            else if (bStat->label[i] == 2){
                                if (bStat->label[j] == 0 || bStat->label[j] == 1){
                                    if (bStat->bundleStart[i] + bStat->bundleServVec[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                                        bStat->bArcs[i][j] = true;
                                    }  
                                }
                            }
                            else{
                                if (bStat->bundleStart[i] + bStat->bundleServVec[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                                    bStat->bArcs[i][j] = true;
                                }                                   
                            }
                        }
                    } 
                }
                if (bStat->label[i] == 0 || bStat->label[i] == 2){
                    bStat->bArcs[i][bStat->bundleVec.size()-1] = true;
                }
            }
            else if (i >= setN){
                for (int j = 0; j < setN; j++){
                    if (bStat->label[j] == 0 || bStat->label[j] == 1){
                        if (bStat->bundleStart[i] + bStat->bundleServVec[i] + (mdist[bStat->lastElement[i]][bStat->firstElement[j]]/inst->vmed) < bStat->bundleStart[j]){
                            bStat->bArcs[i][j] = true;
                        }                        
                    }
                }
            }
        }

        //remove arcs from/to bundles with negative start times
        for (int i = 0; i < bStat->bundleStart.size() - 1; i++){
            if (bStat->bundleStart[i] < 0){
                for (int j = 0; j < bStat->bundleVec.size() - 1; j++){
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

    
}

void feasibleClusterArcs (instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat, vector< vector<int> > &clusterVec, pair<int, int> &cFArc, vector< vector<bool> > &cArcs, vector< vector< pair<int,int> > > &cArcPlus, vector< vector< pair<int,int> > > &cArcMinus, int p, probStat* problem){
    
    int reqClusters = clusterVec.size() - inst->K - 1;
    int clusterA = 0;
    int clusterB;
    int setN = bStat->bundleVec.size() - inst->K - 1;
    int ref;

    if (p < 0){
        if(problem->scen == "1A"){
            ref = inst->m;
        }
        else if (problem->scen == "2A"){
            ref = 3*inst->m;
        }
    }

    else{
        ref = p;
    }

    for (int i = 0; i < bStat->bundleVec.size(); i++){
        if (i < setN){
            if (i > clusterA*(ref + 1) + ref){
                clusterA++;
            }           
        }
        else{
            clusterA++;
        }
        clusterB = 0;
        for (int j = 0; j < bStat->bundleVec.size(); j++){
            if (j < setN){
                if (j > clusterB*(ref + 1) + ref){
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

void makeParcelBundles(instanceStat *inst, vector<nodeStat> &nodeVec, bundleStat *bStat, probStat* problem){
    int parcelReq;
    if (problem->scen == "1A"){
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
    }

    else if (problem->scen == "2A"){
        int setN = bStat->bundleVec.size() - inst->K - 1;
        for (int i = 0; i < setN; i++){
            if (bStat->bundleVec[i].size() <= 1){
                continue;
            }
            else{
                for (int j = 0; j < bStat->bundleVec[i].size(); j++){
                    if (bStat->bundleVec[i][j] < inst->n){
                        continue;
                    }
                    else{
                        parcelReq = bStat->bundleVec[i][j];
                        bStat->parcelBundleVec[parcelReq - inst->n].push_back(i);
                    }
                }                
            }   
        }         
    }
}

void makeStartTimes (instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat){

    double parcelTime;
    double bundleTime;
    bool firstPassenger;

    for (int i = 0; i < bStat->bundleVec.size(); i++){
        parcelTime = 0;
        bundleTime = 0;
        firstPassenger = false;
        if (bStat->bundleVec[i].size() > 1){
            if (bStat->bundleVec[i][0] < inst->n){
                firstPassenger = true;
                bundleTime = nodeVec[bStat->bundleVec[i][0]].e;
                bStat->bundleStart.push_back(bundleTime);
                continue;
            }
            for (int j = 0; j < bStat->bundleVec[i].size() - 1; j++){
                if (bStat->bundleVec[i][j] >= inst->n){
                    parcelTime += ((mdist[bStat->bundleVec[i][j]][bStat->bundleVec[i][j + 1]])/inst->vmed) + nodeVec[bStat->bundleVec[i][j]].delta;

                    if (bStat->bundleVec[i][j + 1] < inst->n){
                        bundleTime = nodeVec[bStat->bundleVec[i][j + 1]].e - parcelTime;
                        firstPassenger = true;
                        break;
                    }
                }
            }
            if (firstPassenger){
                bStat->bundleStart.push_back(bundleTime);
            }            
        }

        else{
            firstPassenger = true;
            bundleTime = nodeVec[bStat->bundleVec[i][0]].e;
            bStat->bundleStart.push_back(bundleTime);
            continue;
        }
    }

}

void makeBundleReference(instanceStat *inst, double **mdist, vector<nodeStat> &nodeVec, bundleStat *bStat){

    for (int i = 0; i < bStat->bundleVec.size(); i++){

        bStat->lastElement.push_back(bStat->bundleVec[i][bStat->bundleVec[i].size()-1]);

    }
    for (int i = 0; i < bStat->bundleVec.size(); i++){

        bStat->firstElement.push_back(bStat->bundleVec[i][0]);

    }
}

bool compareCosts(const bParcelStruct &a, const bParcelStruct &b){
    return a.cost < b.cost;
}

void makeSmallerProblem(instanceStat *inst, vector<nodeStat> &nodeVec, double **mdist, int p, vector< vector<int> > &clsParcel, probStat* problem, int Q){

    if (problem->scen == "1A"){
        vector<bParcelStruct> vecOfDist;
        double singleCost;
        int counter;
        bParcelStruct bps;

        if (p > -1 && p != 0){
            for (int i = 0; i < inst->n; i++){

                for (int j = inst->n; j < inst->n + inst->m; j++){
                    bps.cost = mdist[j][i];
                    bps.parcelreq = j;
                    vecOfDist.push_back(bps);
                }
                counter = 0;
                for (int j = inst->n + inst->m; j < inst->n + 2*inst->m; j++){

                    vecOfDist[counter].cost += mdist[i][j];
                    counter++;
                }

                sort(vecOfDist.begin(), vecOfDist.end(), compareCosts);

                for (int j = 0; j < p; j++){
                    clsParcel[i].push_back(vecOfDist[j].parcelreq);
                }
                vecOfDist.clear();
            }
        }

        else if (p == 0){
            
        }

        else{
            for (int i = 0; i < inst->n; i++){
                for (int j = inst->n; j < inst->n + inst->m; j++){
                    bps.cost = mdist[j][i];
                    bps.parcelreq = j;
                    vecOfDist.push_back(bps);
                }
                counter = 0;
                for (int j = inst->n + inst->m; j < inst->n + 2*inst->m; j++){
                    vecOfDist[counter].cost += mdist[i][j];
                    counter++;
                }
                for (int j = 0; j < inst->m; j++){
                    clsParcel[i].push_back(vecOfDist[j].parcelreq);
                }
                vecOfDist.clear();
            }
        }
    }

    else if (problem->scen == "2A"){
        vector<bParcelStruct> vecOfDist;
        double singleCost;
        bParcelStruct bps;

        if (p > -1 && p != 0){
            //*********************************************
            //find p best of pickup and p best of delivery
            //*********************************************
            // for (int i = 0; i < inst->n; i++){
            //     for (int j = inst->n; j < inst->n + inst->m; j++){
            //         bps.cost = mdist[j][i];
            //         bps.parcelreq = j;
            //         vecOfDist.push_back(bps);
            //     }
            //     sort(vecOfDist.begin(), vecOfDist.end(), compareCosts);

            //     for (int j = 0; j < p; j++){
            //         clsParcel[i].push_back(vecOfDist[j].parcelreq);
            //     }
            //     vecOfDist.clear(); 

            //     for (int j = inst->n + inst->m; j < inst->n + 2*inst->m; j++){
            //         bps.cost = mdist[i][j];
            //         bps.parcelreq = j;
            //         vecOfDist.push_back(bps);
            //         bps.cost = mdist[j][i];
            //         bps.parcelreq = j*inst->m;
            //         vecOfDist.push_back(bps);                   
            //     }
            //     sort(vecOfDist.begin(), vecOfDist.end(), compareCosts);

            //     for (int j = 0; j < p; j++){
            //         clsParcel[i].push_back(vecOfDist[j].parcelreq);
            //     }
            //     vecOfDist.clear(); 
            // }
            //*********************************************
            //find p best of pickup and delivery together
            //*********************************************
            for (int i = 0; i < inst->n; i++){
                for (int j = inst->n; j < inst->n + inst->m; j++){
                    bps.cost = mdist[j][i];
                    bps.parcelreq = j;
                    vecOfDist.push_back(bps);
                }
                for (int j = inst->n + inst->m; j < inst->n + 2*inst->m; j++){
                    bps.cost = mdist[i][j];
                    bps.parcelreq = j;
                    vecOfDist.push_back(bps);
                    bps.cost = mdist[j][i];
                    bps.parcelreq = j*inst->m;
                    vecOfDist.push_back(bps);                   
                }
                sort(vecOfDist.begin(), vecOfDist.end(), compareCosts);

                for (int j = 0; j < p; j++){
                    clsParcel[i].push_back(vecOfDist[j].parcelreq);
                }
                vecOfDist.clear();
            }
        }

        else if (p == 0){
            // break;
        }

        else{
            for (int i = 0; i < inst->n; i++){
                for (int j = inst->n; j < inst->n + inst->m; j++){
                    clsParcel[i].push_back(j);
                }
                vecOfDist.clear(); 

                for (int j = inst->n + inst->m; j < inst->n + 2*inst->m; j++){
                    
                    clsParcel[i].push_back(j);
                    clsParcel[i].push_back(j*inst->m);

                }
                vecOfDist.clear();
            }
        }  
    }    


    // else if (problem->scen == "1B"){

    //     vector< vector<int> > ppuSets;
    //     vector< vector<int> > pdoSets;
    //     vector<int> prep;
    //     vector<bParcelStruct> vecOfDist;
    //     double singleCost;
    //     int counter;
    //     bParcelStruct bps;
        
    //     for (int j = inst->n; j < inst->n + inst->m; j++){
    //         ppuSets.push_back(prep);
    //         pdoSets.push_back(prep);
    //     }
    //     counter = 0;
    //     for (int j = inst->n; j < inst->n + inst->m; j++){
    //         ppu[counter].push_back(j);
    //         for (int k = inst->n; k < inst->n + inst->m; k++){
    //             if(j != k){

                    
    //                 ppu[counter].push_back(k);
    //                 counter++;
    //             }
    //         }
    //     }
    // }



}