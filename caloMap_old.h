#ifndef CALO_MAP_OLD_H
#define CALO_MAP_OLD_H

#include <map>
#include <vector>
#include <tuple>
#include <iostream>
using namespace std;

inline map<pair<int,int>, vector<int>> GetCaloChMapOld() {
    map<pair<int,int>, vector<int>> chMap;

    // Mapping {(MID, ch)} → (left/right, col, row)
    std::map<std::pair<int,int>, std::tuple<int,int,int>> channelMap;

    std::vector<std::tuple<int,int,int>> mid41_map = {
        {0,7,0}, {1,7,0}, {0,7,1}, {1,7,1}, {1,4,0}, {1,4,1}, {1,4,2}, {1,4,3},
        {1,5,0}, {1,5,1}, {1,5,2}, {1,5,3}, {1,6,0}, {1,6,1}, {1,6,2}, {1,6,3},
        {0,7,2}, {1,7,2}, {0,7,3}, {1,7,3}, {0,4,0}, {0,4,1}, {0,4,2}, {0,4,3},
        {0,5,0}, {0,5,1}, {0,5,2}, {0,5,3}, {0,6,0}, {0,6,1}, {0,6,2}, {0,6,3}
    };
    for (int i = 0; i < 32; ++i) channelMap[{41, i+1}] = mid41_map[i];

    std::vector<std::tuple<int,int,int>> mid42_map = {
        {1,0,0}, {1,0,1}, {1,0,2}, {1,0,3}, {1,1,0}, {1,1,1}, {1,1,2}, {1,1,3},
        {1,2,0}, {1,2,1}, {1,2,2}, {1,2,3}, {1,3,0}, {1,3,1}, {1,3,2}, {1,3,3},
        {0,0,0}, {0,0,1}, {0,0,2}, {0,0,3}, {0,1,0}, {0,1,1}, {0,1,2}, {0,1,3},
        {0,2,0}, {0,2,1}, {0,2,2}, {0,2,3}, {0,3,0}, {0,3,1}, {0,3,2}, {0,3,3}
    };
    for (int i = 0; i < 32; ++i) channelMap[{42, i+1}] = mid42_map[i];

    // Convert to the format expected by the code: {lr, mod, col, layer}
    // For old mapping: lr = tuple[0], col = tuple[1], layer = tuple[2]
    for (const auto& kv : channelMap) {
        int mid = kv.first.first;
        int ch = kv.first.second;
        int lr = std::get<0>(kv.second);
        int col = std::get<1>(kv.second);
        int layer = std::get<2>(kv.second);
        
        // For old mapping, we need to map to 1-32 GeomID range
        // MID 41: col 4-7, layer 0-3 -> GeomID 1-16
        // MID 42: col 0-3, layer 0-3 -> GeomID 17-32
        int geomID = layer * 8 + col + 1;
        
        // Convert to the format: {lr, mod, col, layer}
        // For old mapping, use the calculated geomID directly
        chMap[{mid, ch}] = {lr, geomID, col, layer};
        
        // Debug output
        std::cout << "MID=" << mid << ", CH=" << ch 
                  << " -> {lr=" << lr << ", mod=" << geomID 
                  << ", col=" << col << ", layer=" << layer 
                  << "} -> GeomID=" << geomID << std::endl;
    }

    return chMap;
}

#endif // CALO_MAP_OLD_H 