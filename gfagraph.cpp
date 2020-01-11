//
//  gfagraph.cpp
//  BioinformaticsProject
//
//  Created by Ema Puljak on 05/12/2019.
//  Copyright © 2019 Ema Puljak. All rights reserved.
//
#include<assert.h>
#include <stdio.h>
#include <fstream>
#include <sstream>
#include "gfagraph.h"
#include<iostream>

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


GfaGraph GfaGraph::GetSubgraph(const std::unordered_set<int>& ids) const
{
    GfaGraph result;
    result.edgeOverlap = edgeOverlap;
    for (auto node : ids)
    {
        if (nodes.count(node) == 0) continue;
        result.nodes[node] = nodes.at(node);
        NodePos end {node, true};
        if (edges.count(end) == 1)
        {
            for (auto target : edges.at(end))
            {
                if (ids.count(target.id) == 0) continue;
                result.edges[end].push_back(target);
            }
        }
        NodePos start {node, false};
        if (edges.count(start) == 1)
        {
            for (auto target : edges.at(start))
            {
                if (ids.count(target.id) == 0) continue;
                result.edges[start].push_back(target);
            }
        }
    }
    return result;
}


GfaGraph GfaGraph::GetSubgraph(const std::unordered_set<int>& nodeids, const std::unordered_set<std::pair<NodePos, NodePos>>& selectedEdges) const
{
    GfaGraph result;
    result.edgeOverlap = edgeOverlap;
    for (auto node : nodeids)
    {
        if (nodes.count(node) == 0) continue;
        result.nodes[node] = nodes.at(node);
        NodePos end {node, true};
        if (edges.count(end) == 1)
        {
            for (auto target : edges.at(end))
            {
                if (nodeids.count(target.id) == 0) continue;
                if (selectedEdges.count({end, target}) == 1 || selectedEdges.count({target, end}) == 1) result.edges[end].push_back(target);
            }
        }
        NodePos start {node, false};
        if (edges.count(start) == 1)
        {
            for (auto target : edges.at(start))
            {
                if (nodeids.count(target.id) == 0) continue;
                if (selectedEdges.count({start, target}) == 1 || selectedEdges.count({target, start}) == 1) result.edges[start].push_back(target);
            }
        }
    }
    return result;
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
        //std::cout<<line[0]<<" mi je prvi znak u lineu"<<std::endl;
        if (line[0] == 'S')
        {
            //std::cout<<"ok prvo slovo mi je s hhhhhhh"<<std::endl;
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
        int id = 0;
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
                if (j == 0){
                    for (auto iter = parents.begin(); iter != parents.end(); ++iter){
                        NodePos topos {id, true};
                        NodePos fromPos {id_map[iter->id], true};
                        new_graph.edges[fromPos].push_back(topos);
                        new_graph.parents[topos.id].push_back(fromPos);
                    }
                }
                
            }
        }
    }
    
    return result;
}


