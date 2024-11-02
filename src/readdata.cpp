#include "readdata.h"
#include "functions.h"
#include "modelnode.h"
#include "modeltwostage.h"
#include <cstdlib>
#include <stdio.h>

void readData (int argc, char** argv, nodeStat *node, instanceStat *inst, vector<nodeStat> &nodeVec, double ***Mdist, probStat* problem){
    
    if (argc < 4) {
        exit(1);
    }
    
    if (argc > 4) {
        exit(1);
    }  

    inst->preInst = 0;
    inst->instFolder = getOrigin(argv);
    inst->InstName   = getInstName(argv);
    inst->instType   = getInstanceType(argv); // instâcnia sf_data, csarp ou ghsarp
    inst->instModel  = getInstModel(argv);
    problem->model   = argv[3];


    // Convertendo as definiçõpes padrão para minutos
    /*---------------------------------------------------*/
    double kmPerMin = inst->vmed/double(60);
    int beginMin = inst->B*60;
    int endMin = inst->dayEnd*60;
    /*---------------------------------------------------*/


    string file, ewf;

    int n;      // Número de clientes
    int m;      // Número de pacotes
    int K;      // Número de veículos
    int S = -1; // Número de shifts
    
    int R;      // quantidade de requests
    int V;      // quantidade original de nós (desconsiderando dummy)
    int dummy;  // índice inicial dos dummy nodes
    int full;   // total de nós da instância

    int tempNode;   // Variável lixo
    double service; // Tempo de serviço para cada nó
    double T;       // Final da jornada de trabalho
    double B;       // Inicio da jornada de trabalho

    char *instance; 
    instance = argv[1];

    ifstream in(instance, ios::in);
    
    if( !in ) {
        exit (1);
    }

    // Lendo a primeira linha das instâncias
    in >> K >> service >> n >> m;

    R = 2*n + 2*m;
    V = R + K;
    dummy = K;
    full = V + dummy;

    // Criando vectors importantes
    /*---------------------------------------------------*/
    vector<double> vxs;
    vector<double> vys;
    vector<double> vloadCustomer;
    vector<double> vloadParcel;
    vector<double> ve;
    vector<double> vxf;
    vector<double> vyf;
    vector<double> vl;

    inst->service = service;
    inst->n = n;
    inst->m = m;
    inst->K = K;
    /*---------------------------------------------------*/


    // Lendo todos os dados
    /*---------------------------------------------------*/
    resizeStructures(vxs, vys, vloadCustomer, vloadParcel, ve, vl, R);

    if (inst->instFolder == "InstancesZTest") {
        for (int i = 0; i < R; i++){
            in >> tempNode >> vxs[i] >> vys[i] >> vloadCustomer[i] >> vloadParcel[i] >> ve[i] >> vl[i];
        }
    } else if (inst->instType != "sf_data") {
        for (int i = 0; i < R; i++){
            in >> tempNode >> vxs[i] >> vys[i] >> vloadCustomer[i] >> ve[i] >> vl[i];
        }
    } else {
        for (int i = 0; i < R; i++){
            in >> tempNode >> vxs[i] >> vys[i] >> tempNode >> vloadCustomer[i] >> ve[i] >> vl[i];
        }
    }

    if (inst->instFolder == "InstancesZTest") {
        S = max(S, readNewZTestsCsarp(inst, in, tempNode, vxs, vys, vloadCustomer, vloadParcel, ve, vl, R, V));
    } else {
        S = max(S, readDepotCsarp(inst, in, tempNode, vxs, vys, vloadCustomer, vloadParcel, ve, vl, R, V));
        S = max(S, readDepotGhsarp(inst, in, tempNode, vxs, vys, vloadCustomer, vloadParcel, ve, vl, R, V));
        S = max(S, readDepotSf_data(inst, in, tempNode, vxs, vys, vloadCustomer, vloadParcel, ve, vl, R, V));
        vloadParcel = vloadCustomer;
    }

    fillDummy(vxs, vys, vloadCustomer, vloadParcel, ve, vl, S, inst->B, inst->dayEnd);
    /*---------------------------------------------------*/

    in.close();

    int Sdummy = S;
    inst->dummy = Sdummy;

    int sV = R + Sdummy;
    full = sV + Sdummy;

    // Calcula a matriz de distâncias
    /*---------------------------------------------------*/
    double **dist = new double*[full];
    for (int i= 0; i < full; i++){
        dist[i] = new double [full];
    }

    calcDistCsarp(dist, full, V, vxs, vys, vxf, vyf, inst->instType);
    calcDistGhsarp(dist, full, V, vxs, vys, vxf, vyf, inst->instType);
    calcDistSfsarp(dist, full, V, vxs, vys, vxf, vyf, inst->instType);
    /*---------------------------------------------------*/
    
    // Calculando os profits de cada nó
    /*---------------------------------------------------*/
    double *delta = new double[full]; // tempo de serviço de cada nó
    double *profit = new double[full];  // profit de cada nó

    double singleProfit;
    double mandist;

    if (inst->instFolder == "InstancesZTest") {}


    for (int i = 0; i < full; i++){
        delta[i] = service/double(60);
        profit[i] = 0;

        if (inst->instFolder == "InstancesZTest") {
            double discount = 1;

            if (vloadCustomer[i] == 1) {
                mandist = dist[i][i+n]; 
                profit[i] += inst->minpas + inst->paskm*mandist;
            }

            if (vloadParcel[i] == 1) {
                mandist = dist[i][i+m];
                double discount = inst->minpar + inst->parkm*mandist;
                if (profit[i] - 0.00001 > 0) {
                    discount -= inst->parkm*mandist;
                    discount /= 2;
                }

                profit[i] += discount;
            }
        } else {
            if (i < n){ 
                mandist = dist[i][i+n]; 
                profit[i] = inst->minpas + inst->paskm*mandist;
            }
            else if (i < R) {
                if (i < 2*n || i >= 2*n + m){
                    profit[i] = 0;
                }
                else {
                    mandist = dist[i][i+m];
                    profit[i] =  inst->minpar + inst->parkm*mandist;
                }
            }
        }

        if (i >= sV - K){
            delta[i] = 0;
            profit[i] = 0;
        }
    }
    /*---------------------------------------------------*/


    // Folgar a janela de tempo do delivery
    /*---------------------------------------------------*/
    for (int i = n; i < 2*n; i++) {
        vl[i] = inst->dayEnd*60;
    }
    /*---------------------------------------------------*/


    // Ajsutando janela de tempo dos customers
    /*---------------------------------------------------*/
    tightWindowDETOUR1(dist, n, m, ve, vl, kmPerMin, inst->instModel);
    /*---------------------------------------------------*/

    // Preenchendo os dados de nodeVec
    /*---------------------------------------------------*/
     for (int i = 0; i < sV; i++){
        node->xs = vxs[i];
        node->ys = vys[i];
        node->load = vloadCustomer[i];

        if (inst->instFolder != "InstancesZTest") {
            if (i < n){
                node->customerLoad = 1;
            }
            else if (i < 2*n){
                node->customerLoad = -1;
            }
            else{
                node->customerLoad = 0;
            }    

            if (i >= 2*n && i < 2*n + m){
                node->parcelLoad = 1;
            }
            else if (i >= 2*n + m && i < 2*n + 2*m){
                node->parcelLoad = -1;
            }
            else{
                node->parcelLoad = 0;
            }   
        } else {
            node->parcelLoad = vloadParcel[i];
            node->customerLoad = vloadCustomer[i];
        }
                    
        node->e = ve[i]/60;

        node->l = vl[i]/60;

        // cout << i << " " << node->e << " " << node->l << endl;

        node->xf = vxs[i];
        node->yf = vys[i];
        node->delta = delta[i];
        node->profit = profit[i];
        node->index = i;
        nodeVec.push_back(*node);
    }
    // getchar();

    // Adding dummy nodes
    for (int i = 0; i < inst->dummy; i++){
        node->xs = 0;
        node->ys = 0;
        node->load = 0;    
        node->customerLoad = 0;
        node->parcelLoad = 0;    
        node->e = ve[sV + i];
        node->l = vl[sV + i];      
        node->xf = 0;
        node->yf = 0;
        node->delta = 0;
        node->profit = 0;
        node->index = sV + i;
        nodeVec.push_back(*node);
    }
    /*---------------------------------------------------*/

    // Preenchendo o objeto inst e Mdist
    /*---------------------------------------------------*/
    *Mdist = dist;
    inst->K = K;
    inst->n = n;
    inst->m = m;
    inst->V = sV;
    inst->Ks = inst->S.size();
    inst->service = service;
    inst->totalCustomProfit = 0;
    for (int i = 0; i < n; i++){
        inst->totalCustomProfit += nodeVec[i].profit;
    }
    /*---------------------------------------------------*/

    delete[] profit;
    delete[] delta;
}

void fillDummy(vector<double> &vxs, vector<double> &vys, vector<double> &vloadCustomer, vector<double> &vloadParcel, vector<double> &ve, vector<double> &vl, int S, int B, int T) {
    for (int i = 0; i < S; i++){
        vxs.push_back(0);
        vys.push_back(0);
        vloadCustomer.push_back(0);
        vloadParcel.push_back(0);
        ve.push_back(ve[ve.size() - S]/60.0);
        vl.push_back(vl[vl.size() - S]/60.0);
    }
}

void resizeStructures(vector<double> &vxs, vector<double> &vys, vector<double> &vloadCustomer, vector<double> &vloadParcel, vector<double> &ve, vector<double> &vl, int _size) {
    vxs.resize(vxs.size() + _size, 0);
    vys.resize(vys.size() + _size, 0);
    vloadCustomer.resize(vloadCustomer.size() + _size, 0);
    vloadParcel.resize(vloadParcel.size() + _size, 0);
    ve.resize(ve.size() + _size, 0);
    vl.resize(vl.size() + _size, 0);
}

int readDepotCsarp(instanceStat *inst, ifstream &in, int tempNode, vector<double> &vxs, vector<double> &vys, vector<double> &vloadCustomer, vector<double> &vloadParcel, vector<double> &ve, vector<double> &vl, int startDepot, int startDummy) {
    if (inst->instType != "csarp") {
        return -1;
    }

    resizeStructures(vxs, vys, vloadCustomer, vloadParcel, ve, vl, startDummy - startDepot);

    int shiftIndex = 0;
    int shiftDepotIndex = 2*inst->n + 2*inst->m;

    for (int i = startDepot; i < startDummy; i++){

        // Each vehicle in old csarp has only one shift
        inst->S.push_back(shiftIndex++);
        inst->vehicleShifts.push_back(vector<int>(1, shiftDepotIndex++));

        in >> tempNode >> vxs[i] >> vys[i] >> vloadCustomer[i] >> ve[i] >> vl[i];
    }

    return shiftIndex;
}

int readDepotGhsarp(instanceStat *inst, ifstream &in, int tempNode, vector<double> &vxs, vector<double> &vys, vector<double> &vloadCustomer, vector<double> &vloadParcel, vector<double> &ve, vector<double> &vl, int startDepot, int startDummy) {
    if (inst->instType != "ghsarp") {
        return -1;
    }

    resizeStructures(vxs, vys, vloadCustomer, vloadParcel, ve, vl, startDummy - startDepot);

    in >> tempNode >> vxs[startDepot] >> vys[startDepot] >> vloadCustomer[startDepot] >> ve[startDepot] >> vl[startDepot];

    // Each vehicle in old ghsarp has only one shift
    int shiftIndex = 0;
    int shiftDepotIndex = 2*inst->n + 2*inst->m;

    inst->S.push_back(shiftIndex++);
    inst->vehicleShifts.push_back(vector<int>(1, shiftDepotIndex++));

    for (int i = startDepot + 1; i < startDummy; i++) {
        vxs[i] = vxs[startDepot];
        vys[i] = vys[startDepot];
        vloadCustomer[i] = vloadCustomer[startDepot];
        ve[i] = ve[startDepot];
        vl[i] = vl[startDepot];

        inst->S.push_back(shiftIndex++);
        inst->vehicleShifts.push_back(vector<int>(1, shiftDepotIndex++));
    }
    
    return shiftIndex;
}

int readDepotSf_data(instanceStat *inst, ifstream &in, int tempNode, vector<double> &vxs, vector<double> &vys, vector<double> &vloadCustomer, vector<double> &vloadParcel, vector<double> &ve, vector<double> &vl, int startDepot, int startDummy) {
    if (inst->instType != "sf_data") {
        return -1;
    }

    resizeStructures(vxs, vys, vloadCustomer, vloadParcel, ve, vl, 1);
    
    in >> tempNode >> vxs[startDepot] >> vys[startDepot] >> tempNode >> vloadCustomer[startDepot] >> ve[startDepot] >> vl[startDepot];

    rotate(vxs.begin(), vxs.begin() + 1, vxs.end());
    rotate(vys.begin(), vys.begin() + 1, vys.end());
    rotate(vloadCustomer.begin(), vloadCustomer.begin() + 1, vloadCustomer.end());
    rotate(vloadParcel.begin(), vloadParcel.begin() + 1, vloadParcel.end());
    rotate(ve.begin(), ve.begin() + 1, ve.end());
    rotate(vl.begin(), vl.begin() + 1, vl.end());

    rotate(vxs.begin() + inst->n, vxs.begin() + inst->n + inst->m, vxs.begin() + 2*inst->n + inst->m);
    rotate(vys.begin() + inst->n, vys.begin() + inst->n + inst->m, vys.begin() + 2*inst->n + inst->m);
    rotate(vloadCustomer.begin() + inst->n, vloadCustomer.begin() + inst->n + inst->m, vloadCustomer.begin() + 2*inst->n + inst->m);
    rotate(vloadParcel.begin() + inst->n, vloadParcel.begin() + inst->n + inst->m, vloadParcel.begin() + 2*inst->n + inst->m);
    rotate(ve.begin() + inst->n, ve.begin() + inst->n + inst->m, ve.begin() + 2*inst->n + inst->m);
    rotate(vl.begin() + inst->n, vl.begin() + inst->n + inst->m, vl.begin() + 2*inst->n + inst->m); 

    resizeStructures(vxs, vys, vloadCustomer, vloadParcel, ve, vl, startDummy - startDepot - 1);

    // Each vehicle in old ghsarp has only one shift
    int shiftIndex = 0;
    int shiftDepotIndex = 2*inst->n + 2*inst->m;

    inst->S.push_back(shiftIndex++);
    inst->vehicleShifts.push_back(vector<int>(1, shiftDepotIndex++));

    for (int i = startDepot + 1; i < startDummy; i++){
        vxs[i] = vxs[startDepot];
        vys[i] = vys[startDepot];
        vloadCustomer[i] = vloadCustomer[startDepot];
        vloadParcel[i] = vloadParcel[startDepot];
        ve[i] = ve[startDepot];
        vl[i] = vl[startDepot];

        inst->S.push_back(shiftIndex++);
        inst->vehicleShifts.push_back(vector<int>(1, shiftDepotIndex++));
    }

    return shiftIndex;
}

int readNewZTestsCsarp(instanceStat *inst, ifstream &in, int tempNode, vector<double> &vxs, vector<double> &vys, vector<double> &vloadCustomer, vector<double> &vloadParcel, vector<double> &ve, vector<double> &vl, int startDepot, int startDummy) {
    int n;
    
    // Each vehicle can have multiple shifts
    int shiftIndex = 0;
    int shiftDepotIndex = 2*inst->n + 2*inst->m;
 
    while (in >> n) {
        inst->vehicleShifts.push_back(vector<int>());

        for (int i = 0; i < n; i++) {
            resizeStructures(vxs, vys, vloadCustomer, vloadParcel, ve, vl, 1);
            vloadParcel.resize(vloadParcel.size() + 1);

            in >> tempNode >> vxs[shiftDepotIndex] >> vys[shiftDepotIndex] >> vloadCustomer[shiftDepotIndex] >> vloadParcel[i] >> ve[shiftDepotIndex] >> vl[shiftDepotIndex];

            inst->S.push_back(shiftIndex++);
            inst->vehicleShifts.back().push_back(shiftDepotIndex++);
        }
    }

    return shiftIndex;
}

// Read the shifts of the vehicles in the instances
// K is the total number of vehicles
// N is the total number of requests
vector<vector<pair<double, double>>> vehicleShifts(const int K, ifstream &in) {
    vector<vector<pair<double, double>>> shifts;

    shifts.resize(K);

    for (int i = 0; i < K; i++) {
        int n;
        in >> n;

        for (int j = 0; j < n; j++) {
            shifts[i].push_back(make_pair(-1, -1));
            in >> shifts[i][j].first >> shifts[i][j].second;
            shifts[i][j].first /= 60.0;
            shifts[i][j].second /= 60.0;
        }
    }

    return shifts;
}

// Ajusta as janelas de tempo para o modo DETOUR1
void tightWindowDETOUR1(double **dist, int n, int m, vector<double> &ve, vector<double> &vl, double kmPerMin, string instModel) {
    if (instModel != "DETOUR1") {
        return;
    }

    for (int i = n; i < 2*n; i++){
        ve[i] = ve[i-n] + dist[i-n][i]/kmPerMin + double(5);

        double maxArrive = ve[i];

        // Na pior hipótese, um nó de pacote será inserido no serviço do customer
        for (int j = 2*n; j < 2*n + 2*m; j++) {
            double newTime = ve[i-n] + double(5) + dist[i-n][j]/kmPerMin;
            newTime += double(5) + dist[j][i]/kmPerMin;

            maxArrive = max(maxArrive, newTime);
        }

        vl[i] = maxArrive;
        vl[i] += 0.0001;
    }
}

void calcDistCsarp(double **dist, int full, int V, const vector<double> &vxs, const vector<double> &vys, const vector<double> &vxf, const vector<double> &vyf, string instType) {
    if (instType != "csarp") {
        return;
    }

    double manDist;

    for (int i = 0; i < full; i++) {
        for (int j = 0; j < full; j++){
            if (i < V && j < V && i != j) {
                manDist = calcEucDist2(vxs, vys, vxs, vys, i, j);              

                dist[i][j] = manDist;
            } else {
                dist[i][j] = 0;
            }
        }
    }
}

void calcDistGhsarp(double **dist, int full, int V, const vector<double> &vxs, const vector<double> &vys, const vector<double> &vxf, const vector<double> &vyf, string instType) {
    if (instType != "ghsarp") {
        return;
    }

    double manDist;
    double conversor = 50;

    for (int i = 0; i < full; i++) {
        for (int j = 0; j < full; j++){
            if (i < V && j < V && i != j) {
                manDist = calcEucDist2(vxs, vys, vxs, vys, i, j);              

                dist[i][j] = manDist/conversor;
            } else {
                dist[i][j] = 0;
            }
        }
    }
}

void calcDistSfsarp(double **dist, int full, int V, const vector<double> &vxs, const vector<double> &vys, const vector<double> &vxf, const vector<double> &vyf, string instType) {
    if (instType != "sf_data") {
        return;
    }

    double manDist;

    for (int i = 0; i < full; i++) {
        for (int j = 0; j < full; j++){
            if (i < V && j < V && i != j) {
                manDist = CalcMan(vxs, vys, vxs, vys, i, j);              

                dist[i][j] = manDist;
            } else {
                dist[i][j] = 0;
            }
        }
    }
}