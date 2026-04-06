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

//モチーフじゃなくて転写因子にした方がよかった
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

struct fasta {
    string name;
    string seq;
};

int main() {
    //各転写因子の対数オッズスコアの計算
    vector <motif> motifs;

    for(int i = 0; i < motif_names.size(); i++){
        ifstream file("data/motif/" + motif_names.at(i));

        string seq;
        int len_seq;
        int num_seq = 0;
        //各コンセンサス配列についてwhile
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
        //対数オッズスコア計算完了
        file.close();
    }

    //各プロモータ配列の取得
    vector <fasta> promoter;
    string line;
    fasta current;

    ifstream file ("data/seq/promoters");

    while(getline(file,line)){
        if(line[0] == '>'){
            if(!current.name.empty()) promoter.push_back(current);
            current.name = line.substr(1);
            current.seq = "";
        }
        else current.seq += line;
    }
    promoter.push_back(current);
    file.close();
    //各プロモータ配列の取得終了
    
    //各モチーフについて全プロモータ配列を回す
    for(int i = 0; i < motifs.size(); i++){
        cout << "Motif:" << motifs.at(i).motif_name << endl << endl;
    }
}