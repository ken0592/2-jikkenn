#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
using namespace std;

double freq_A = 7519429.0 / (7519429 * 2 + 4637676 * 2);
double freq_C = 4637676.0 / (7519429 * 2 + 4637676 * 2);
double freq_G = 4637676.0 / (7519429 * 2 + 4637676 * 2);
double freq_T = 7519429.0 / (7519429 * 2 + 4637676 * 2);

vector <double> q = {freq_A, freq_C, freq_G, freq_T};

vector <string> motif_names = {"MATa1", "MATalpha2", "MCM1", "MIG1", "PHO4", "RCS1", "ROX1", "TAF"};

struct motif {
    string motif_name;
    int motif_size;
    vector<vector<double>> p;
    vector<vector<double>> s;

    motif (string name, int size){
        motif_name = name;
        motif_size = size;
        p.resize(4, vector<double>(size, 1.0));
        s.resize(4, vector<double>(size, 0.0));
    }
};

vector <motif> motifs;

int main() {
    for(int i = 0; i < motif_names.size(); i++){
        ifstream file("data/motif/" + motif_names.at(i));

        string seq;
        int len_seq;
        int num_seq = 0;

        while (getline(file, seq)) {
            if(num_seq == 0){
                len_seq = seq.size();
                motifs.push_back(motif(motif_names.at(i), len_seq));
            }
            num_seq++;
            for(int j = 0; j < len_seq; j++){
                if(seq.at(j) == 'A') motifs.at(i).p.at(0).at(j)++;
                else if(seq.at(j) == 'C') motifs.at(i).p.at(1).at(j)++;
                else if(seq.at(j) == 'G') motifs.at(i).p.at(2).at(j)++;
                else if(seq.at(j) == 'T') motifs.at(i).p.at(3).at(j)++;
            }
        }

        for(int j = 0; j < 4; j++){
            for(int k = 0; k < len_seq; k++){
                motifs.at(i).p.at(j).at(k) = motifs.at(i).p.at(j).at(k) / (num_seq+4);
                motifs.at(i).s.at(j).at(k) = log(motifs.at(i).p.at(j).at(k) / q.at(j));
            }
        }
        file.close();
    }

    for(int j = 0; j < 4; j++){
        for(int k = 0; k < 10; k++){
            cout << motifs.at(0).s.at(j).at(k) << " ";
        }
        cout << endl;
    }
}