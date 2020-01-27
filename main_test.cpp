#include <stdio.h>
#include <iostream>
#include "gfagraph.h"
#include "fastqloader.h"
#include <fstream>
#include <vector>       // std::vector
#include <algorithm>    // std::min_element, std::max_element
#include "chrono"
#include<climits>
using namespace std;

vector<int> current_row;
vector<int> next_row;
vector<NodePos> parents;

inline const char * const BoolToString(bool b)
{
  return b ? "true" : "false";
}

std::ostream &operator<<(std::ostream &os, NodePos const &m) {
    return os << m.id << '/' << BoolToString(m.end);
}

//comparator we send to the *min_element function
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

int Navaro (string A, GfaGraph B, int A_n, int B_n, int* arguments){
    string query_char;
    string current_node;
    unordered_map<int, std::string> nodes = B.nodes;
    //add the room for the j at the beginning of every row
    current_row.resize(B_n + 1);
    for(int i=0; i < A_n; i++){
        query_char = A[i];
        if (i==0){
            // init the first row 0
            fill(current_row.begin(),current_row.end(),0);
        }
        next_row.push_back(i+1);
        for(int j=1; j < B_n + 1; j++){
            current_node = B.nodes.at(j);
            if(B.parents.find(j) != B.parents.end()){ parents = B.parents.at(j);}
            //node_state_func calculates the result for every node in the matrix
            next_row.push_back(node_state_func(current_node, query_char, j, arguments));
        }
        //"swap" the rows in memory and go next
        current_row = next_row;
        next_row.erase(next_row.begin(), next_row.end());
        parents.erase(parents.begin(), parents.end());
    }
    return *min_element(current_row.begin(), current_row.end());
}
int main (int argc, char** argv){
    //The first 3 arguments are costs for match / mis / indel
    int arguments [3];
    for (int i = 0; i < 3; ++i){
        arguments[i] = atoi(argv[i+1]);
    }
    //Reading sequences from the 4th argument
    string sequence = argv[4];
    cout << sequence;

    // reading the graph from the gfa file
    std::string filename = argv[5];
    GfaGraph graph = GfaGraph::LoadFromFile(filename);
    
    
    std::cout << "Graph nodes contains:";
    for ( auto it = graph.nodes.begin(); it != graph.nodes.end(); ++it )
      std::cout << " " << it->first << ":" << it->second;
    std::cout << std::endl;
    
    std::cout << "Graph edges contains:";
    std::cout << graph.edges.size();
    for ( auto it = graph.edges.begin(); it != graph.edges.end(); ++it ){
        std::cout << " " << it->first << "  Nodovi s kojima sam povezan:";
        for(int i=0; i<it->second.size(); ++i){
            std::cout << " " << it->second[i].id << ' ';
        }
      std::cout << std::endl;
    }
    std::cout << std::endl;
    
    

    //The display variable determines whether you print the distance and time required
    //at the end of every sequences calculation or you write it in the _result.txt file
    bool Display = std::stoi(argv[6]);

    //Graph information
    int num_nodes = (int)graph.nodes.size();
    int num_edges = (int)graph.edges.size();
    //We store the results in this variable
    string results = "";
    int VMPR;
    auto start = chrono::high_resolution_clock::now();
    //runs the algorithm over test sequence
    
    auto start1 = chrono::high_resolution_clock::now();
    // Calculate the distance between the graph and the sequence
    VMPR = Navaro(sequence, graph, sequence.size(), num_nodes, arguments);
    //time measuring variables
    auto stop1 = chrono::high_resolution_clock::now();
    auto duration1 = chrono::duration_cast<chrono::microseconds>(stop1-start1);
    string sequence_result = "Sequence size: " + to_string(sequence.size()) + ", score: " + to_string(VMPR) + ", execution time: " + to_string(duration1.count()) + " us\n";
    if(Display == true) cout << sequence_result;
    else results.append(sequence_result);
    auto stop = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(stop-start);
    cout << "Duration for all sequences: " << duration.count() << " microseconds" << endl;
    
    long duration_time = duration.count();
    // Writing results to a _results.txt file in the same directory
    ofstream outputFile;
    string outputname = filename.substr(0, filename.size()-4) + "_result.txt";
    outputFile.open(outputname);
    outputFile << "Number of nodes: " << num_nodes << endl;
    outputFile << "Number of edges: " << num_edges << endl;
    outputFile << "Duration of algoritm for test sequence: " << duration_time << " microseconds" << endl;
    if(Display == false)
        outputFile << results << endl;
    outputFile.close();
    return 0;
}


