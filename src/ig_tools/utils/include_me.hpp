#pragma once

#include "iostream"
#include "fstream"
#include "sstream"
#include "string"
#include "vector"
#include "map"
#include "set"
#include "assert.h"
#include <algorithm>
#include "math.h"

using namespace std;

inline size_t abs_diff(size_t a, size_t b) {
	if(a > b)
		return a - b;
	return b - a;
}
