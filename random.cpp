#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>

using namespace std;

double freq_A = 7519429.0 / (7519429 * 2 + 4637676 * 2);
double freq_C = 4637676.0 / (7519429 * 2 + 4637676 * 2);
double freq_G = 4637676.0 / (7519429 * 2 + 4637676 * 2);
double freq_T = 7519429.0 / (7519429 * 2 + 4637676 * 2);

vector <double> q = {freq_A, freq_C, freq_G, freq_T};

int main(){ 
    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<int> dist({q.at(0),q.at(1),q.at(2),q.at(3)});

    string bases = "ACGT";
    int length = 500;
    vector <string> results(100);

    ofstream fout ("random_seq.txt");
    for(int i = 0; i < 100; i++){
        for(int j = 0; j < length; j++){
            int index = dist(gen);
            results.at(i) += bases.at(index);
        }
        fout << results.at(i) << endl;
    }
     fout.close();
}