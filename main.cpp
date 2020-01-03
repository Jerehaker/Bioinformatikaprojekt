//
//  SW_vol1.cpp
//  BioinformaticsProject
//
//  Created by Ema Puljak on 05/12/2019.
//  Copyright © 2019 Ema Puljak. All rights reserved.
//

#include "SW_vol1.hpp"
#include <iostream>
#include "gfagraph.h"
#include "fastqloader.h"
#include <fstream>
#include <vector>       // std::vector
#include <algorithm>    // std::min_element, std::max_element
#include "chrono"
#include <climits>
using namespace std;

vector<int> current_row;
vector<int> next_row;
vector<NodePos> parents;

bool comp(int a, int b)
{
    return (a < b);
}

int node_state_func(string current_node, string query_char, int node_id, int* arguments){
    //int match = 0;
    //int mis = 1;
    //int indel = 0;
    
    int match = arguments[0];
    int mis = arguments[1];
    int indel = arguments[2];
    //vector<int> C_v;
    
    switch(query_char.compare(current_node)){
        case 0:
            switch(parents.size()){
                case 0:
                    return match + current_row[node_id - 1];
                    break;
                default:
                    int min_num = INT_MAX;
                    for (NodePos parent : parents){
                        min_num = min(min_num, current_row[parent.id] + match);
                        //C_v.push_back(current_vec[parent.id] + match);
                    }
                    //return *min_element(C_v.begin(), C_v.end());
                    return min_num;
            }
            break;
        default:
            int min_num;
            switch(parents.size()){
                case 0:
                    /*C_v.push_back(mis + current_vec[node_id - 1]);
                    C_v.push_back(current_vec[node_id] + indel);
                    C_v.push_back(next_vec[node_id - 1] + indel);
                    return *min_element(C_v.begin(), C_v.end());*/
                    return min({mis + current_row[node_id - 1], current_row[node_id] + indel, next_row[node_id - 1] + indel}, comp);
                    break;
                default:
                    //C_v.push_back(current_vec[node_id] + indel);
                    min_num = current_row[node_id] + indel;
                    // mismatch
                    for(NodePos parent : parents){
                        min_num = min({min_num, (current_row[parent.id] + mis), (next_row[parent.id] + indel)}, comp);
                        //C_v.push_back(current_vec[parent.id] + mis);
                        //C_v.push_back(next_vec[parent.id] + indel);
                    }
                    //return *min_element(C_v.begin(), C_v.end());
                    return min_num;
            }
    }
}

int Navarov (string A, GfaGraph B, int A_n, int B_n, int* arguments){
    string query_char;
    string current_node;

    unordered_map<int, std::string> nodes = B.nodes;
    current_row.resize(B_n + 1);

    for(int i=0; i < A_n; i++){
        query_char = A[i];
        if (i==0){
            // init prvi red na 0
            fill (current_row.begin(),current_row.end(),0);
            /*for(int j=0; j <= B_n; j++){
                current.push_back(0);
            }*/
        }
        next_row.push_back(i+1);
        for(int j=1; j < B_n; j++){
            current_node = B.nodes.at(j);
            if(B.parents.find(j) != B.parents.end()){ parents = B.parents.at(j);}
            next_row.push_back(node_state_func(current_node, query_char, j, arguments));
        }
        
        //nadi min vrijednost u nextu
        //vector_min_price_index_row.push_back(getMin(next).first);
        //zamini next i current; next init
        current_row = next_row;
        // obrisi next
        next_row.erase(next_row.begin(), next_row.end());
    }
    return *min_element(next_row.begin(), next_row.end());
}

int main (int argc, char** argv){
    
    //ulazni agrumenti
    int arguments [3];
    for (int i = 0; i < 3; ++i){
        arguments[i] = atoi(argv[i+1]);
    }
    
    // ucitavanje sekvenci iz fastQ file-a
    std::string filename_fastq = argv[4];
    cout << argv[4];
    std::vector<FastQ> seq = loadFastqFromFile(filename_fastq);
    
    // ucitavanje grafa iz gfa file-a
    std::string filename = argv[5];
    GfaGraph graph = GfaGraph::LoadFromFile(filename);
    
    //podaci o grafu
    int num_nodes = (int)graph.nodes.size();
    int num_edges = (int)graph.edges.size();
    int B_n = num_nodes;
    
    auto start = chrono::high_resolution_clock::now();
    //for petlja po sekvencama
    for(int i=0; i < seq.size(); ++i){
        string A = seq[i].sequence;
        int A_n = A.size();
        auto start1 = chrono::high_resolution_clock::now();
        // pokretanje algoritma
        int VMPR = Navarov(A, graph, A_n, B_n, arguments);
        auto stop1 = chrono::high_resolution_clock::now();
        auto duration1 = chrono::duration_cast<chrono::microseconds>(stop1-start1);
        //cout<<"broj sekvence: "<< i <<" duljina: "<< A.size() << " udaljenost je: " << VMPR << " vrijeme: " << duration1.count() << " microsec"<< endl;
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop-start);
    cout << duration.count() << " microseconds";
    std::cout << std::endl;
    
    long duration_time = duration.count();
    
    // zapisi u file
    ofstream outputFile;
    string outputname = filename.substr(0, filename.size()-4) + "_result.txt";
    outputFile.open(outputname);
    outputFile << "Number of nodes: " << num_nodes << endl;
    outputFile << "Number of edges: " << num_edges << endl;
    outputFile << "Duration of algoritm for 74 sequences: " << duration_time << "us" << endl;
    outputFile.close();

    return 0;
}
