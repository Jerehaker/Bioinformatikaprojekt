//
//  fastqloader.h
//  BioinformaticsProject
//
//  Created by Ema Puljak on 05/12/2019.
//  Copyright © 2019 Ema Puljak. All rights reserved.
//

#ifndef fastqloader_h
#define fastqloader_h

#include <string>
#include <vector>

class FastQ {
public:
    FastQ reverseComplement() const;
    std::string seq_id;
    std::string sequence;
    std::string quality;
};

std::vector<FastQ> loadFastqFromFile(std::string filename);

#endif /* fastqloader_h */
