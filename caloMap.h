// caloMap.h
// CALO channel geometrical mapping header-only implementation
// Simply include this in any ROOT macro, QA script, or C++ analysis code to access the unified map
// 수정 시 이 파일만 바꾸면 모든 코드에서 자동 적용됨

#ifndef CALO_MAP_H
#define CALO_MAP_H

#include <map>
#include <vector>
using namespace std;

inline map<pair<int,int>, vector<int>> GetCaloChMap() {
    map<pair<int,int>, vector<int>> chMap;

    // MID 41
    chMap[{41, 3}] = {0,  6, 7, 0};
    chMap[{41, 4}] = {1,  6, 7, 0};
    chMap[{41, 6}] = {1, 26, 4, 0};
    chMap[{41, 7}] = {1, 30, 4, 1};
    chMap[{41, 8}] = {1, 29, 4, 2};
    chMap[{41,10}] = {1, 18, 5, 0};
    chMap[{41,11}] = {1, 17, 5, 1};
    chMap[{41,12}] = {1, 31, 5, 2};
    chMap[{41,14}] = {1, 27, 6, 0};
    chMap[{41,15}] = {1, 22, 6, 1};
    chMap[{41,16}] = {1, 28, 6, 2};
    chMap[{41,18}] = {1,  7, 7, 1};
    chMap[{41,20}] = {1, 15, 7, 2};
    chMap[{41,22}] = {0, 26, 4, 0};
    chMap[{41,23}] = {0, 30, 4, 1};
    chMap[{41,24}] = {0, 29, 4, 2};
    chMap[{41,26}] = {0, 18, 5, 0};
    chMap[{41,27}] = {0, 17, 5, 1};
    chMap[{41,28}] = {0, 31, 5, 2};
    chMap[{41,30}] = {0, 27, 6, 0};
    chMap[{41,31}] = {0, 22, 6, 1};
    chMap[{41,32}] = {0, 28, 6, 2};

    // MID 42
    chMap[{42, 5}] = {1, 14, 1, 0};
    chMap[{42, 6}] = {1, 13, 1, 1};
    chMap[{42, 7}] = {1, 21, 1, 2};
    chMap[{42, 8}] = {1, 10, 1, 3};
    chMap[{42, 9}] = {1,  4, 2, 0};
    chMap[{42,10}] = {1, 25, 2, 1};
    chMap[{42,11}] = {1, 23, 2, 2};
    chMap[{42,12}] = {1, 32, 2, 3};
    chMap[{42,13}] = {1, 11, 3, 0};
    chMap[{42,14}] = {1, 33, 3, 1};
    chMap[{42,15}] = {1, 20, 3, 2};
    chMap[{42,16}] = {1, 19, 3, 3};
    chMap[{42,21}] = {0, 14, 1, 0};
    chMap[{42,22}] = {0, 13, 1, 1};
    chMap[{42,23}] = {0, 21, 1, 2};
    chMap[{42,24}] = {0, 10, 1, 3};
    chMap[{42,25}] = {0,  4, 2, 0};
    chMap[{42,26}] = {0, 25, 2, 1};
    chMap[{42,27}] = {0, 23, 2, 2};
    chMap[{42,28}] = {0, 32, 2, 3};
    chMap[{42,29}] = {0, 11, 3, 0};
    chMap[{42,30}] = {0, 33, 3, 1};
    chMap[{42,31}] = {0, 20, 3, 2};
    chMap[{42,32}] = {0, 19, 3, 3};

    return chMap;
}

#endif // CALO_MAP_H
