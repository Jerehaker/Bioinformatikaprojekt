#include <stdio.h>
#include <iostream>
#include "gfagraph.h"
#include "fastqloader.h"
#include <fstream>
#include <vector>       // std::vector
#include <algorithm>    // std::min_element, std::max_element
#include "chrono"
using namespace std;

vector<int> current_row;
vector<int> next_row;
vector<NodePos> parents;

bool comp(int a, int b)
{
    return (a < b);
}

int node_state_func(string current_node, string query_char, int node_id, int* arguments){
    int match = arguments[0];
    int mis = arguments[1];
    int indel = arguments[2];
    
    switch(query_char.compare(current_node)){
        case 0:
            switch(parents.size()){
                case 0:
                    return match + current_row[node_id - 1];
                    break;
                default:
                    int min_num = INT_MAX;
                    for (NodePos parent : parents){
                        if(parent.id < node_id){
                            min_num = min(min_num, current_row[parent.id] + match);
                        }
                    }
                    return min_num;
            }
            break;
        default:
            int min_num;
            switch(parents.size()){
                case 0:
                    return min({mis + current_row[node_id - 1], current_row[node_id] + indel, next_row[node_id - 1] + indel}, comp);
                    break;
                default:
                    min_num = current_row[node_id] + indel;
                    // mismatch
                    for(NodePos parent : parents){
                        min_num = min({min_num, (current_row[parent.id] + mis), (next_row[parent.id] + indel)}, comp);
                    }
                    return min_num;
            }
    }
}

int Navarov (string A, GfaGraph B, int A_n, int B_n, int* arguments){
    string query_char;
    string current_node;

    unordered_map<int, std::string> nodes = B.nodes;
    current_row.resize(B_n + 1);
    cout << "New sequence\n";
    for(int i=0; i < A_n; i++){
        query_char = A[i];
        if (i==0){
            // init prvi red na 0
            fill (current_row.begin(),current_row.end(),0);
        }
        next_row.push_back(i+1);
        for(int j=1; j < B_n; j++){
            current_node = B.nodes.at(j);
            if(B.parents.find(j) != B.parents.end()){ parents = B.parents.at(j);}
            next_row.push_back(node_state_func(current_node, query_char, j, arguments));
        }

        current_row = next_row;
        next_row.erase(next_row.begin(), next_row.end());
        parents.erase(parents.begin(), parents.end());
    }
    return *min_element(current_row.begin(), current_row.end());
}

int main (int argc, char** argv){
    
    //ulazni agrumenti
    int arguments [3];
    for (int i = 0; i < 3; ++i){
        arguments[i] = atoi(argv[i+1]);
    }
    
    // ucitavanje sekvenci iz fastQ file-a
    std::string filename_fastq = argv[4];
    std::vector<FastQ> seq = loadFastqFromFile(filename_fastq);
    
    // ucitavanje grafa iz gfa file-a
    std::string filename = argv[5];
    GfaGraph graph = GfaGraph::LoadFromFile(filename);
    
    /*std::cout << "Graph nodes contains:";
    for ( auto it = graph.nodes.begin(); it != graph.nodes.end(); ++it )
      std::cout << " " << it->first << ":" << it->second;
    std::cout << std::endl;*/
    
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
        cout<<"Number of sequence: "<< i <<", size of sequence: "<< A.size() << ", distance of alignment: " << VMPR << ", time of execution: " << duration1.count() << " microsec"<< endl;
    }
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop-start);
    cout << duration.count() << " microseconds";
    std::cout << std::endl;
    
    long duration_time = duration.count();
    
    // zapisi u file
    ofstream outputFile;
    outputFile.open("results_twopath_graph.txt");
    outputFile << "Number of nodes: " << num_nodes << endl;
    outputFile << "Number of edges: " << num_edges << endl;
    outputFile << "Duration of algoritm for 74 sequences: " << duration_time << " microseconds" << endl;
    outputFile.close();

    return 0;
}
