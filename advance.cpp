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

vector <vector<double>> p(4, vector<double>(L,1));
vector <vector<double>> s(4, vector<double>(L,0));



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

    uniform_int_distribution<int> dist(0,promoters.size()-1);
    int N = dist(gen); 

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
        for(auto& k : p.at(i)){
            k /= (promoters.size()-1+4);
        }
    }

    double s_dif = 0;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < L; j++){
            s_dif += abs(s.at(i).at(j) - log(p.at(i).at(j)/q.at(i)););
            s.at(i).at(j) = log(p.at(i).at(j)/q.at(i));
        }
    }
    
    vector<double> scores;
    for(int i = 0; i < L; i++){
        if(promoters.at(N).seq.at())
    }


    ofstream fout("adv_output.txt");
    /*for(int i = 0; i < promoters.size(); i++){
        fout << promoters.at(i).name << endl;
        fout << promoters.at(i).seq << endl;
        fout << promoters.at(i).pos << endl << endl;
    }*/

    /*for(int i = 0; i < 4; i++){
        for(double& k : s.at(i)){
            fout << k << " ";
        }
        fout << endl;
    }*/
    fout.close();
}