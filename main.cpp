#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
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

double p = 7.50341;

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
    ofstream fout ("output.txt");
    for(int i = 0; i < motifs.size(); i++){
        fout << "Motif:" << motifs.at(i).motif_name << endl << endl;
        for(int j = 0; j < promoter.size(); j++){
            for(int k = 0; k < promoter.at(j).seq.size() - motifs.at(i).motif_size + 1; k++){
                double score = 0.0;
                for(int l = 0; l < motifs.at(i).motif_size; l++){
                    if(promoter.at(j).seq.at(k+l) == 'A') score += motifs.at(i).s.at(0).at(l);
                    else if(promoter.at(j).seq.at(k+l) == 'C') score += motifs.at(i).s.at(1).at(l);
                    else if(promoter.at(j).seq.at(k+l) == 'G') score += motifs.at(i).s.at(2).at(l);
                    else if(promoter.at(j).seq.at(k+l) == 'T') score += motifs.at(i).s.at(3).at(l);
                }
                if(score > p){
                    fout << "Promoter:" << promoter.at(j).name << endl;
                    fout << "Pos:" << k+1 << "~" << k+motifs.at(i).motif_size << endl;
                    fout << "hit(" << promoter.at(j).seq.substr(k, motifs.at(i).motif_size) << ") = " << score << endl << endl;
                }
            }
        }
        fout << "---------------------------------" << endl;
    }
    fout.close();

    //閾値決定
    vector<double> scores;
    ifstream file2("random_seq.txt");
    for(int i = 0; i < motifs.size(); i++){
        while(getline(file2, line)){
                for(int j = 0; j < 500 - motifs.at(i).motif_size + 1; j++){
                    double score = 0.0;
                    for(int k = 0; k < motifs.at(i).motif_size; k++){
                        if(line.at(j+k) == 'A') score += motifs.at(i).s.at(0).at(k);
                        else if(line.at(j+k) == 'C') score += motifs.at(i).s.at(1).at(k);
                        else if(line.at(j+k) == 'G') score += motifs.at(i).s.at(2).at(k);
                        else if(line.at(j+k) == 'T') score += motifs.at(i).s.at(3).at(k);
                    }
                    scores.push_back(score);
                }
            }
    }
    file2.close();

    sort(scores.begin(), scores.end());
    ofstream fout2 ("border.txt");
    fout2 << scores.at(scores.size() - scores.size()/10000) << endl;
    fout2 << scores.size();
    fout2.close();
}