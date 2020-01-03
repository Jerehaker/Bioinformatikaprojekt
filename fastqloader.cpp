//
//  fastqloader.cpp
//  BioinformaticsProject
//
//  Created by Ema Puljak on 05/12/2019.
//  Copyright © 2019 Ema Puljak. All rights reserved.
//

#include <stdio.h>
#include <algorithm>
#include <fstream>
#include "fastqloader.h"

std::vector<FastQ> loadFastqFastqFromFile(std::string filename)
{
    std::ifstream file {filename};
    std::vector<FastQ> result;
    do
    {
        std::string line;
        std::getline(file, line);
        if (line[0] != '@') continue;
        FastQ newread;
        if (line.back() == '\r') line.pop_back();
        newread.seq_id = line.substr(1);
        std::getline(file, line);
        if (line.back() == '\r') line.pop_back();
        newread.sequence = line;
        std::getline(file, line);
        std::getline(file, line);
        if (line.back() == '\r') line.pop_back();
        newread.quality = line;
        result.push_back(newread);
    } while (file.good());
    return result;
}

std::vector<FastQ> loadFastqFromFile(std::string filename)
{
    if (filename.substr(filename.size()-6) == ".fastq") return loadFastqFastqFromFile(filename);
    return std::vector<FastQ>{};
}


