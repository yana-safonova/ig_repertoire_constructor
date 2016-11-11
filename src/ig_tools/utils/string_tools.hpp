#pragma once

#include "include_me.hpp"
#include <seqan/seq_io.h>

inline vector<string> split(const string &s, char delim) {
    vector<string> elems;
    stringstream ss(s);
    string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

inline pair<string, string> split_by_dots(string str) {
    size_t del_pos = str.find("..");
    assert(del_pos != string::npos);
    return make_pair(str.substr(0, del_pos),
                     str.substr(del_pos + 2, str.size() - del_pos - 2));
}

template<typename T>
T string_to_number(string str) {
    stringstream ss;
    ss << str;
    T num;
    ss >> num;
    return num;
}

inline pair<size_t, size_t> string_range_to_number(pair<string, string> range_str) {
    return make_pair(string_to_number<size_t>(range_str.first), string_to_number<size_t>(range_str.second));
}

template<typename T>
string number_to_string(T n) {
    stringstream ss;
    ss << n;
    return ss.str();
}

// TODO: rename
inline string delete_spaces(string str) {
    for (size_t i = 0; i < str.size(); i++)
        if (str[i] == ' ')
            str[i] = '_';
    return str;
}

inline vector<string> split(const string &s, const string &separator) {
    vector<string> ans;
    size_t sep_len = separator.length();
    size_t found = s.find(separator);
    size_t previous = 0;
    while (found != string::npos) {
        ans.push_back(s.substr(previous, found - previous));
        previous = found + sep_len;
        found = s.find(separator, found + sep_len);
    }
    ans.push_back(s.substr(previous));
    return ans;
}

inline void string_to_khashes(string &s, size_t k, vector<size_t> &ans) {
    size_t a = 239;
    size_t p = 1;
    size_t s_length = s.length();
    if (k > s_length)
        return;
    ans.resize(0);
    ans.reserve(s_length - k + 1);
    size_t khash = 0;
    for (size_t i = 0; i < k; i++) {
        khash *= a;
        khash += s[i];
        p *= a;
    }
    ans.push_back(khash);
    for (size_t i = k; i < s_length; i++) {
        khash *= a;

        khash -= s[i - k] * p;
        khash += s[i];
        ans.push_back(khash);
    }
    return;
}

template<typename T>
string seqan_string_to_string(const T &s) {
    char buffer[length(s) + 1];
    strcpy(buffer, seqan::String<char, seqan::CStyle>(s));
    return string((char *) buffer);
}