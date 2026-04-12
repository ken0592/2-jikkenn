#include <iostream>
#include <vector>
#include <fstream>

using namespace std; 

struct promoter{
    string name;
    string seq;
    int pos;
};

int L = 7;

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

    ofstream fout("adv_output.txt");
    for(int i = 0; i < promoters.size(); i++){
        fout << promoters.at(i).name << endl;
        fout << promoters.at(i).seq << endl << endl;
    }
}