#pragma once

#include "include_me.hpp"

vector<string> split(const string &s, char delim);

pair<string, string> split_by_dots(string str);

template<typename T>
T string_to_number(string str);

pair<size_t, size_t> string_range_to_number(pair<string, string> range_str);

template<typename T>
string number_to_string(T n);

string delete_spaces(string str);

vector<string> split(const string &s, const string &separator);

void StringToKHashes(string &s, size_t k, vector<size_t> &ans);
