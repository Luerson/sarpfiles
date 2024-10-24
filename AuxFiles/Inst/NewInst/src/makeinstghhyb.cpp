#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <list>
#include <vector>
#include <algorithm>	
#include <iterator>
#include <math.h>
#include <cmath>
#include <limits>
#include <float.h>
#include <iomanip>
#include <ctime>
#include <iostream>
#include <locale.h>
#include <sys/time.h>
#include <unistd.h>
#include <cstdlib>
#include <stdio.h>

using namespace std;


double CalcDistEuc (double X1, double Y1, double X2, double Y2);

void readData (int argc, char** argv, double ***Mdist);

string getInstanceType (char **argv, int i);

string getInstanceType2 (char **argv, int i);

double CalcDistEuc (double X1, double Y1, double X2, double Y2){
    return sqrt ( pow ( X1 - X2, 2 ) + pow ( Y1 - Y2, 2 ) );
}

void readData (int argc, char** argv, double ***Mdist){

    if (argc < 3) {
        cout << "\nMissing parameters\n";
        cout << " ./exeHybIns [Instance1] [Instance2]" << endl;
        exit(1);
    }
    
    if (argc > 3) {
        cout << "\nToo many parameters\n";
        cout << " ./exeHybIns [Instance1] [Instance2]" << endl;
        exit(1);
    }  

 //   	vector<coordSt> coordVec1;
	// vector<coordSt> coordVec2;

	string file1;
	string file2;

    int s1;
    int s2;

    pair<double, double> coordinate;

    vector< pair <double, double> > coordVec;

    int refpoint = 2;

    int n, m;

    int seed = 1234;
    srand(seed);

    vector <vector <double> > tempData;
    vector<double> auxtempdata;

    vector< pair <double, double> > start;

    vector <vector <double> > realData;

    vector<double> service;

    vector< vector<double> > dist;
    vector< vector<double> > realdist;

    vector<double> rowvec;
    char *instance1; 
    char *instance2; 

    string it1;
    string it2;

    instance1 = argv[1];
    instance2 = argv[2];

    it1 = getInstanceType(argv, 1);
    it2 = getInstanceType2(argv, 2);

    // cout << "\nit1: " << it1 << " - it2: " << it2 << endl;
    // getchar();



    if (it1 != it2){
        ifstream in1(instance1, ios::in);

        
        if( !in1 ) {
            cout << "the file could not be opened\n";
            exit (1);
        }

        while ( file1.compare("DIMENSION:") != 0 && file1.compare("DIMENSION") != 0 ){
            in1 >> file1;
        }

        in1 >> s1;

        n = (s1 - 2)/2;

        while (file1.compare("NODE_COORD_SECTION") != 0){
            in1 >> file1;
        }

        for (int i = 0; i < s1 - 2; i++){
            coordVec.push_back(coordinate);
        }


        for (int i = 0; i < s1; i++){

            if (i == 0){
                in1 >> file1;
                in1 >> coordinate.first;
                in1 >> coordinate.second;
                start.push_back(coordinate);
            }

            else if (i == 1){
                in1 >> file1;
                in1 >> file1;
                in1 >> file1;
            }
            
            else{
                in1 >> file1;
                in1 >> coordVec[i - 2].first;
                in1 >> coordVec[i - 2].second;

            }

        }

        ifstream in2(instance2, ios::in);

        
        if( !in2 ) {
            cout << "the file could not be opened\n";
            exit (1);
        }

        while ( file2.compare("DIMENSION:") != 0 && file2.compare("DIMENSION") != 0 ){
            in2 >> file2;
        }

        in2 >> s2;

        // m = (s2 - 2)/2;
        m = 3;
        s2 = m*2 + 2;

        while (file2.compare("NODE_COORD_SECTION") != 0){
            in2 >> file2;
        }

        for (int i = 0; i < s2 - 2; i++){
            coordVec.push_back(coordinate);
        }

        for (int i = s1 - 2; i < s1+s2 - 2; i++){

            if (i == s1 - 2){
                in2 >> file2;
                in2 >> coordinate.first;
                in2 >> coordinate.second;
                start.push_back(coordinate);
            }
            else if (i == s1 - 1){
                in2 >> file2;
                in2 >> file2;
                in2 >> file2;
            }

            else{
                in2 >> file2;
                in2 >> coordVec[i - 2].first;
                in2 >> coordVec[i - 2].second;            
            }
        }

        int count;
        pair<double, double> chosen;

        //rearranging pickup and delivery
        count = 1;
        for (int i = 2; i < 2*m; i++){
            if (i % 2 == 0){
                chosen = coordVec[2*n + i];
                coordVec.erase(coordVec.begin() + 2*n + i);
                coordVec.insert(coordVec.begin() + 2*n + count, chosen);
                count++;
            }
        }

        for(int i = 0; i < start.size(); i++){
            coordVec.push_back(start[i]);
        }

        for (int i = 0; i < coordVec.size(); i++){
            if (i < n){
                service.push_back(0);
            }
            for (int j = 0; j < coordVec.size(); j++){
                rowvec.push_back(0);
            }
            dist.push_back(rowvec);
            rowvec.clear();
        }

        for (int i = 0; i < n + 2*m + 2; i++){
            for (int j = 0; j < n + 2*m + 2; j++){
                rowvec.push_back(0);
            }
            realdist.push_back(rowvec);
            rowvec.clear();
        }


        for (int i = 0; i < 2*n; i++){
            if(i % 2 == 0){
                service[i/2] = floor( CalcDistEuc ( coordVec[i].first, coordVec[i].second, coordVec[i+1].first, coordVec[i+1].second ) + 0.5 );
            }
        }

        for (int i = 0; i < coordVec.size(); i++){
            for (int j = 0; j < coordVec.size(); j++){
                dist[i][j] = floor ( CalcDistEuc (coordVec[i].first, coordVec[i].second, coordVec[j].first, coordVec[j].second) + 0.5 );                   
            }
        }
        cout << endl;

        //Shrink matrix

        for (int i = 0; i < n + 2*m + 2; i++){
            for (int j = 0; j < n + 2*m + 2; j++){
                if (i == j){
                    realdist[i][j] = 0;
                }
                else{
                    if (i < n){
                        if (j < n){
                            realdist[i][j] = (dist[2*i+1][2*j]);
                        }
                        else{
                            realdist[i][j] = (dist[2*i+1][n+j]);
                        }
                    }
                    else{
                        if (j < n){
                            realdist[i][j] = (dist[n+i][2*j]);
                        }
                        else{
                            realdist[i][j] = (dist[n+i][n+j]);
                        }
                    }
                }
            }
        }


        //adding dummy node

        for (int i = 0; i < realdist.size(); i++){
            realdist[i].push_back(0);
            rowvec.push_back(0);
        }

        rowvec.push_back(0);
        realdist.push_back(rowvec);

        // //to screen

        // cout << "DIMENSION: " << n + 2*m + 2;
        // cout << "\nN: " << n;
        // cout << "\nM: " << m;
        // cout << "\nK: " << 2;

        // cout << endl;

        // cout << "\nSERVICE: " << endl;
        // for (int i = 0; i < service.size(); i++){
        //     cout << service[i] << " ";
        // }
        // cout << endl;

        // cout << "\nDIST MATRIX: " << endl;

        // for (int i = 0; i < realdist.size(); i++){
        //     for(int j = 0; j < realdist[i].size(); j++){
        //         cout << setw(5) << realdist[i][j] << " "; 
        //     }
        //     cout << endl;
        // }
        // getchar();

        //output
        string outputname;

        outputname = "grubhub-" + it1 + "-" + it2 + ".tsp";
        // cout << "output: " << outputname << endl;

        ofstream ofile;

        ofile.open(outputname);
        
        ofile << "DIMENSION: " << n + 2*m + 2;
        ofile << "\nN: " << n;
        ofile << "\nM: " << m;
        ofile << "\nK: " << 2;

        ofile << endl;

        ofile << "\nSERVICE: " << endl;
        for (int i = 0; i < service.size(); i++){
            ofile << service[i] << " ";
        }
        ofile << endl;

        ofile << "\nDIST MATRIX: " << endl;

        for (int i = 0; i < realdist.size(); i++){
            for(int j = 0; j < realdist[i].size(); j++){
                ofile << setw(5) << realdist[i][j] << " "; 
            }
            ofile << endl;
        }

    }


}

string getInstanceType (char **argv, int i){

    string filename(argv[i]);

    string::size_type loc = filename.find_first_of("-");

    string::size_type loc2 = filename.find_last_of(".", filename.size());
    string InstanceType;

    InstanceType.append(filename, loc+1, loc2-loc-1 );

    return InstanceType;
}

string getInstanceType2 (char **argv, int i){

    string filename(argv[i]);

    string::size_type loc = filename.find_last_of("-");

    string::size_type loc2 = filename.find_last_of(".", filename.size());
    string InstanceType;
    
    InstanceType.append("03-");
    
    InstanceType.append(filename, loc+1, loc2-loc-1 );


    return InstanceType;
}

int main (int argc, char *argv[]) {

	double **distMatrix;
	double **auxdist;
	double cost;

	
	readData(argc, argv, &distMatrix);

return 0;

}

