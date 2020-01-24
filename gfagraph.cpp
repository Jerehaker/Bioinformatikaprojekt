#include <stdio.h>
#include <fstream>
#include <sstream>
#include "assert.h"
#include "gfagraph.h"
#include <iostream>
#include <cstring>
#include "vector"
#include <algorithm>
#include <unordered_map>

bool hasEnding (std::string const &fullString, std::string const &ending) {
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}


// returns reverse complement of the given string
std::string getReverseComplement (std::string const &sequence) {
    std::string reverse_complement = "";
    for (int i = sequence.size() - 1; i >= 0; i--) {

        char character = sequence.c_str()[i];
        char new_character;
        switch (character) {
            case 'A':
                new_character = 'T';
                break;
            case 'T':
                new_character = 'A';
                break;
            case 'C':
                new_character = 'G';
                break;
            case 'G':
                new_character = 'C';
                break;
        }
        reverse_complement.push_back(new_character);
    }
    return std::string(reverse_complement);
}


NodePos::NodePos() :
id(0),
end(false)
{
}

NodePos::NodePos(int id, bool end) :
id(id),
end(end)
{
}

bool NodePos::operator==(const NodePos& other) const
{
    return id == other.id && end == other.end;
}

bool NodePos::operator!=(const NodePos& other) const
{
    return !(*this == other);
}

NodePos NodePos::Reverse() const
{
    return NodePos { id, !end };
}

GfaGraph::GfaGraph() :
nodes(),
edges(),
edgeOverlap(-1)
{
}

GfaGraph GfaGraph::LoadFromFile(std::string filename)
{
    GfaGraph result;
    std::ifstream file {filename};
    while (file.good())
    {
        std::string line;
        std::getline(file, line);
        if (!file.good()) break;
        if (line.size() == 0) continue;
        if (line[0] != 'S' && line[0] != 'L') continue;
        if (line[0] == 'S')
        {
            std::stringstream sstr {line};
            int id;
            std::string dummy;
            std::string seq;
            sstr >> dummy;
            assert(dummy == "S");
            sstr >> id;
            sstr >> seq;
            result.nodes[id] = seq;
        }
        if (line[0] == 'L')
        {
            std::stringstream sstr {line};
            int from;
            int to;
            std::string fromstart;
            std::string toend;
            std::string dummy;
            int overlap;
            sstr >> dummy;
            assert(dummy == "L");
            sstr >> from;
            sstr >> fromstart;
            sstr >> to;
            sstr >> toend;
            sstr >> overlap;
            assert(overlap >= 0);
            assert(result.edgeOverlap == -1 || overlap == result.edgeOverlap);
            result.edgeOverlap = overlap;
            NodePos frompos {from, fromstart == "+"};
            NodePos topos {to, toend == "+"};
            result.edges[frompos].push_back(topos);
            result.parents[topos.id].push_back(frompos);
        }
    }

    if (hasEnding(filename, "snp.gfa")){
        GfaGraph new_graph;
        new_graph.edgeOverlap = result.edgeOverlap;

        // id of nodes in new_graph
        int id = 1;

        // map contains <key, value> pairs
        // key represents id of result node
        // value represents id of new_graph node
        std::unordered_map<int, int> id_map;
        for (int i = 0; i < result.nodes.size(); i++){
            // temp is string of graph node
            std::vector<NodePos> parents = result.parents[i];
            std::string temp = result.nodes[i];
            for (int j = 0; j < strlen(temp.c_str()); j++){
                // create new node for one char and save
                // connection to old graph node id
                new_graph.nodes[id] = temp.c_str()[j];
                id_map[i] = id;

                // for first character make connections to previous nodes
                // nodes which represent the parents in old graph and last
                // character node made from that exact parent
                if (j == 0){
                    for (auto iter = parents.begin(); iter != parents.end(); ++iter){
                        NodePos topos {id, true};
                        NodePos fromPos {id_map[iter->id], true};
                        new_graph.edges[fromPos].push_back(topos);
                        new_graph.parents[topos.id].push_back(fromPos);
                    }
                }

                // for all others make connection to character (j - 1) that
                // appeared before j character
                else {
                    NodePos topos {id, true};
                    NodePos fromPos {id - 1, true};
                    new_graph.edges[fromPos].push_back(topos);
                    new_graph.parents[topos.id].push_back(fromPos);
                }
                id += 1;
            }
        }
        return new_graph;

    }

    if (hasEnding(filename, "tangle.gfa")){
        GfaGraph new_graph;
        new_graph.edgeOverlap = result.edgeOverlap;

        int id = 1;

        // map contains <key, value> pairs
        // key represents id of result node
        // value represents id of new_graph node
        std::unordered_map<int, int> id_map;
        for (int i = 0; i < result.nodes.size(); i++){
            // temp is string of graph node
            std::vector<NodePos> parents = result.parents[i];
            std::string temp = result.nodes[i];

            for (int j = 0; j < strlen(temp.c_str()); j++){
                // create new node for one char and save
                // connection to old graph node_id times two
                new_graph.nodes[id] = temp.c_str()[j];
                id_map[2 * i] = id;

                // for all but first make connection to character (j - 1) that
                // appeared before j character
                if (j != 0) {
                    NodePos topos {id, true};
                    NodePos fromPos {id - 1, true};
                    new_graph.edges[fromPos].push_back(topos);
                    new_graph.parents[topos.id].push_back(fromPos);
                }
                id += 1;
            }

            std::string reverse_complement = getReverseComplement(temp);
            int length = strlen(reverse_complement.c_str());

            for (int j = 0; j < length; j++){
                // create new node for one char and save connection
                // to old graph node_id times two plus one
                new_graph.nodes[id] = reverse_complement.c_str()[j];
                id_map[2 * i + 1] = id;

                if (j != 0) {
                    NodePos topos {id, true};
                    NodePos fromPos {id - 1, true};
                    new_graph.edges[fromPos].push_back(topos);
                    new_graph.parents[topos.id].push_back(fromPos);
                }
                id += 1;
            }
        }

        id = 1;
        for (int i = 0; i < result.nodes.size(); i++){
            std::vector<NodePos> parents = result.parents[i];
            std::string temp = result.nodes[i];
            for (int j = 0; j < strlen(temp.c_str()); j++){
                // for character after overlap make connections to previous nodes
                // nodes which represent the parents in old graph and last
                // character node made from that exact parent
                if (j == new_graph.edgeOverlap){
                    for (int z = 0; z < parents.size(); z++){
                        std::vector<NodePos> edges_of_parent = result.edges[parents[z]];

                        // look if parent has edge to old_graph temporary node
                        if (std::find(edges_of_parent.begin(), edges_of_parent.end(), NodePos {i, true}) != edges_of_parent.end()) {

                            // if parent end is true make connection
                            // to parent_id times two
                            if (parents[z].end == true ) {
                                NodePos topos{id, true};
                                NodePos fromPos{id_map[parents[z].id * 2], true};
                                new_graph.edges[fromPos].push_back(topos);
                                new_graph.parents[topos.id].push_back(fromPos);
                            }

                            // make connection to parent_id times two plus one
                            else {
                                NodePos topos {id, true};
                                NodePos fromPos {id_map[parents[z].id * 2 + 1], true};
                                new_graph.edges[fromPos].push_back(topos);
                                new_graph.parents[topos.id].push_back(fromPos);
                            }
                        }
                    }
                }
                id += 1;
            }

            std::string reverse_complement = getReverseComplement(temp);
            int length = strlen(reverse_complement.c_str());

            for (int j = 0; j < length; j++){
                // for character after overlap make connections to previous nodes
                // nodes which represent the parents in old graph and last
                // character node made from that exact parent
                if (j == new_graph.edgeOverlap){
                    for (int z = 0; z < parents.size(); z++){
                        std::vector<NodePos> edges_of_parent = result.edges[parents[z]];

                        // look if parent has edge to old_graph temporary node
                        if (std::find(edges_of_parent.begin(), edges_of_parent.end(), NodePos {i, false}) != edges_of_parent.end()) {
                            // if parent end is true make connection
                            // to parent_id times two
                            if (parents[z].end == true ) {
                                NodePos topos{id, true};
                                NodePos fromPos{id_map[parents[z].id * 2], true};
                                new_graph.edges[fromPos].push_back(topos);
                                new_graph.parents[topos.id].push_back(fromPos);
                            }
                            // else make connection to parent_id times two plus one
                            else {
                                NodePos topos {id, true};
                                NodePos fromPos {id_map[parents[z].id * 2 + 1], true};
                                new_graph.edges[fromPos].push_back(topos);
                                new_graph.parents[topos.id].push_back(fromPos);
                            }
                        }
                    }
                }
                id += 1;
            }
        }

        return new_graph;
    }

    return result;
}
