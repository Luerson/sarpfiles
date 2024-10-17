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
    inst->InstName = getInstName(argv);
    inst->instModel = getInstModel(argv);
    problem->model = argv[3];


    // Convertendo as definiçõpes padrão para minutos
    /*---------------------------------------------------*/
    double kmPerMin = inst->vmed/double(60);
    int beginMin = inst->B*60;
    int endMin = inst->T*60;
    /*---------------------------------------------------*/


    string file, ewf;

    int n; // Número de clientes
    int m; // Número de pacotes
    int K; // Número de veículos
    
    int R;     // quantidade de requests
    int V;     // quantidade original de nós (desconsiderando dummy)
    int dummy; // índice inicial dos dummy nodes
    int full;  // total de nós da instância

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

    inst->instType = getInstanceType(argv); // instâcnia sf_data, csarp ou ghsarp

    // Lendo os dados da instância
    in >> K >> service >> n >> m;

    R = 2*n + 2*m;
    V = R + K;
    dummy = K;
    full = V + dummy;

    // Atualizando os dados da instância na struct
    inst->service = service;
    inst->dummy = dummy;
    inst->n = n;
    inst->m = m;
    inst->K = K;


    // Lendo todos os dados de requests
    /*---------------------------------------------------*/
    vector<double> vxs;
    vector<double> vys;
    vector<double> vload;
    vector<double> ve;
    vector<double> vxf;
    vector<double> vyf;
    vector<double> vl;

    for (int i = 0; i < R; i++){
        vxs.push_back(0);
        vys.push_back(0);
        vload.push_back(0);
        ve.push_back(0);
        vl.push_back(0);
    }
    /*---------------------------------------------------*/

    // Lendo todos os dados em csarp
    /*---------------------------------------------------*/
    for (int i = 0; i < R; i++){
        in >> tempNode >> vxs[i] >> vys[i] >> vload[i] >> ve[i] >> vl[i];
    }

    if (inst->instType == "csarp") {
        for (int i = R; i < V; i++){
            vxs.push_back(0);
            vys.push_back(0);
            vload.push_back(0);
            ve.push_back(0);
            vl.push_back(0);
        }

        for (int i = R; i < V; i++){
            in >> tempNode >> vxs[i] >> vys[i] >> vload[i] >> ve[i] >> vl[i];
        }
    }
    /*---------------------------------------------------*/


    // Lendo todos os dados em ghsarp ou sf_Data
    /*---------------------------------------------------*/

    if (inst->instType == "ghsarp" || inst->instType == "sf_data") {
        for (int i = 0; i < R; i++){
            in >> tempNode >> vxs[i] >> vys[i] >> tempNode >> vload[i] >> ve[i] >> vl[i];
        }

        for (int i = R; i < V; i++){
            vxs.push_back(0);
            vys.push_back(0);
            vload.push_back(0);
            ve.push_back(0);
            vl.push_back(0);
        }

        if (inst->instType == "ghsarp") {
            in >> tempNode >> vxs[R] >> vys[R] >> vload[R] >> ve[R] >> vl[R];

            for (int i = R + 1; i < V; i++){
                in >> tempNode >> vxs[i] >> vys[i] >> vload[i] >> ve[i] >> vl[i];
            }
        } else {
            in >> tempNode >> vxs[R] >> vys[R] >> tempNode >> vload[R] >> ve[R] >> vl[R];
            for (int i = R + 1; i < V; i++){
                in >> tempNode >> vxs[i] >> vys[i] >> tempNode >> vload[i] >> ve[i] >> vl[i];
            }
        }
    }
    /*---------------------------------------------------*/


    // Criando os dados dos dummy nodes
    /*---------------------------------------------------*/
    for (int i = 0; i < dummy; i++){
        vxs.push_back(0);
        vys.push_back(0);
        vload.push_back(0);
        ve.push_back(inst->B);
        vl.push_back(inst->T);
    }
    /*---------------------------------------------------*/

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

    for (int i = 0; i < full; i++){
        delta[i] = service/double(60);
        if (i < n){ 
            mandist = dist[i][i+n]; 
            profit[i] = inst->minpas + inst->paskm*mandist;
        }
        else if (i < R){ 
            if (i < 2*n || i >= 2*n + m){
                profit[i] = 0;
            }
            else {
                mandist = dist[i][i+m];
                profit[i] =  inst->minpar + inst->parkm*mandist;
            }
        }
        else if (i >= V - K){
            delta[i] = 0;
            profit[i] = 0;
        }
    }
    /*---------------------------------------------------*/


    // Folgar a janela de tempo do delivery
    /*---------------------------------------------------*/
    for (int i = n; i < 2*n; i++) {
        vl[i] = inst->T*60;
    }
    /*---------------------------------------------------*/


    // Ajsutando janela de tempo dos customers
    /*---------------------------------------------------*/
    tightWindowDETOUR1(dist, n, m, ve, vl, kmPerMin, inst->instModel);
    /*---------------------------------------------------*/

    // Preenchendo os dados de nodeVec
    /*---------------------------------------------------*/
     for (int i = 0; i < V; i++){
            node->xs = vxs[i];
            node->ys = vys[i];
            node->load = vload[i];

            if (i < n){
                node->load2 = -1;
            }
            else if (i < n + m){
                node->load2 = 1;
            }
            else{
                node->load2 = 0;
            }    
                        
            node->e = ve[i]/60;

            node->l = vl[i]/60;

            node->xf = vxs[i];
            node->yf = vys[i];
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
            node->e = inst->B;
            node->l = inst->T;         
            node->xf = 0;
            node->yf = 0;
            node->delta = 0;
            node->profit = 0;
            node->index = V + i;
            nodeVec.push_back(*node);
        }
    /*---------------------------------------------------*/


    // Preenchendo o objeto inst e Mdist
    /*---------------------------------------------------*/
    *Mdist = dist;
    inst->K = K;
    inst->n = n;
    inst->m = m;
    inst->V = V;
    inst->service = service;
    inst->totalCustomProfit = 0;
    for (int i = 0; i < n; i++){
        inst->totalCustomProfit += nodeVec[i].profit;
    }
    /*---------------------------------------------------*/

    delete[] profit;
    delete[] delta;
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