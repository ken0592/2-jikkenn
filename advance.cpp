#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <cmath>

using namespace std; 

double freq_A = 7519429.0 / (7519429 * 2 + 4637676 * 2);
double freq_C = 4637676.0 / (7519429 * 2 + 4637676 * 2);
double freq_G = 4637676.0 / (7519429 * 2 + 4637676 * 2);
double freq_T = 7519429.0 / (7519429 * 2 + 4637676 * 2);

vector <double> q = {freq_A, freq_C, freq_G, freq_T};

struct promoter{
    string name;
    string seq;
    int pos;
};

int L = 7;

vector <vector<double>> p(4, vector<double>(L,1.0));
vector <vector<double>> s(4, vector<double>(L,0.0));

int main(){
    vector <promoter> promoters;
    string line;
    promoter current;

    ifstream file("data/REB1.promoters");
    while(getline(file, line)){
        if(line.at(0) == '>'){
            if(!current.name.empty()) promoters.push_back(current);
            current.seq = "";
            current.name = line.substr(1);
        }
        else current.seq += line;
    }
    promoters.push_back(current);
    file.close();

    random_device rd;
    mt19937 gen(rd());

    for(int i = 0; i < promoters.size(); i++){
        uniform_int_distribution<int> dist(0,promoters.at(i).seq.size()-L);
        promoters.at(i).pos = dist(gen);
    }
    int count = 0;
    while(true){
        count++;
        uniform_int_distribution<int> dist(0,promoters.size()-1);
        int N = dist(gen); 

        for(int i = 0; i < 4; i++){
            for(int j = 0; j < L; j++){
                p.at(i).at(j) = 1.0;
            }
        }

        for(int i = 0; i < promoters.size(); i++){
            if(i == N) continue;
            for(int j = 0; j < L; j++){
                if(promoters.at(i).seq.at(promoters.at(i).pos+j) == 'A') p.at(0).at(j)++;
                else if(promoters.at(i).seq.at(promoters.at(i).pos+j) == 'C') p.at(1).at(j)++;
                else if(promoters.at(i).seq.at(promoters.at(i).pos+j) == 'G') p.at(2).at(j)++;
                else if(promoters.at(i).seq.at(promoters.at(i).pos+j) == 'T') p.at(3).at(j)++;
            }
        }

        for(int i = 0; i < 4; i++){
            for(double& k : p.at(i)){
                k /= (promoters.size()-1+4);
            }
        }

        double s_dif = 0;
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < L; j++){
                s_dif += abs(s.at(i).at(j) - log(p.at(i).at(j)/q.at(i)));
                s.at(i).at(j) = log(p.at(i).at(j)/q.at(i));
            }
        }
        //cout << "s_dif" <<  s_dif << endl;
        //if(s_dif < 0.0001) break;
        if(count == 100000
        ) break;
        
        vector <double> scores(promoters.at(N).seq.size()-L+1,0.0);
        for(int i = 0; i < promoters.at(N).seq.size()-L+1; i++){
            for(int j = 0; j < L; j++){
                if(promoters.at(N).seq.at(i+j) == 'A') scores.at(i) += s.at(0).at(j);
                else if(promoters.at(N).seq.at(i+j) == 'C') scores.at(i) += s.at(1).at(j);
                else if(promoters.at(N).seq.at(i+j) == 'G') scores.at(i) += s.at(2).at(j);
                else if(promoters.at(N).seq.at(i+j) == 'T') scores.at(i) += s.at(3).at(j);
            }
        }
        for(auto &k: scores){
            k = exp(k);
        }
        discrete_distribution<int> dist2(scores.begin(), scores.end());
        promoters.at(N).pos = dist2(gen);
    }

    cout << count << endl;

    ofstream fout("adv_output.txt");
    /*for(int i = 0; i < promoters.size(); i++){
        fout << promoters.at(i).name << endl;
        fout << promoters.at(i).seq << endl;
        fout << promoters.at(i).pos << endl << endl;
    }*/

    fout << "--------対数オッズ表-------" << endl;
    for(int i = 0; i < 4; i++){
        for(double& k : s.at(i)){
            fout << k << " ";
        }
        fout << endl;
    }
    fout << "--------------------------" << endl;

    double score;
    for(int i = 0; i < promoters.size(); i++){
        fout << "promoter:" << promoters.at(i).name << endl << endl;
        for(int j = 0; j < promoters.at(i).seq.size()-L+1; j++){
            score = 0.0;
            for(int k = 0; k < L; k++){
                if(promoters.at(i).seq.at(j+k) == 'A') score += s.at(0).at(k);
                else if(promoters.at(i).seq.at(j+k) == 'C') score += s.at(1).at(k);
                else if(promoters.at(i).seq.at(j+k) == 'G') score += s.at(2).at(k);
                else if(promoters.at(i).seq.at(j+k) == 'T') score += s.at(3).at(k);
            }
            if(score > 4.0){
                fout << "Pos:" << j+1 << "~" << j+L << endl;
                fout << "hit(" << promoters.at(i).seq.substr(j,L) << ") = " << score << endl;
            }
        }
        fout << "--------------------------" << endl;
    }
    fout.close();
}