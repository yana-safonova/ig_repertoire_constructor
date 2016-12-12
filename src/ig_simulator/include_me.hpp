#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include "stdlib.h"
#include "stdio.h"
#include <memory>
#include "utils/sequence_tools.hpp"

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <time.h>

size_t RandomInt(size_t max) {
    return rand() % max;
}

size_t RandomIndex(size_t to, size_t from = 0) {
    return rand() % (to - from) + from;
}

char GetRandomNucleotide() {
    int rand_int = rand() % 4;
    if(rand_int == 0)
        return 'A';
    if(rand_int == 1)
        return 'C';
    if(rand_int == 2)
        return 'G';
    return 'T';
}

char GetAnotherRandomNucleotide(const string &seq, size_t pos) {
    char rand_nucl = GetRandomNucleotide();
    while(seq[pos] == rand_nucl)
        rand_nucl = GetRandomNucleotide();
    return rand_nucl;
}

template<typename T>
T StringToType(string str) {
    stringstream ss(str);
    T t;
    ss >> t;
    return t;
}