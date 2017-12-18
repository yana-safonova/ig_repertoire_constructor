#pragma once

#include <verify.hpp>
#include <seqan/seq_io.h>

inline std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

inline std::pair<std::string, std::string> split_by_dots(std::string str) {
    size_t del_pos = str.find("..");
    VERIFY(del_pos != std::string::npos);
    return make_pair(str.substr(0, del_pos),
                     str.substr(del_pos + 2, str.size() - del_pos - 2));
}

template<typename T>
T string_to_number(std::string str) {
    std::stringstream ss;
    ss << str;
    T num;
    ss >> num;
    return num;
}

inline std::pair<size_t, size_t> string_range_to_number(std::pair<std::string, std::string> range_str) {
    return std::make_pair(string_to_number<size_t>(range_str.first), string_to_number<size_t>(range_str.second));
}

template<typename T>
std::string number_to_string(T n) {
    std::stringstream ss;
    ss << n;
    return ss.str();
}

// TODO: rename
inline std::string delete_whitespaces(std::string str) {
    for (size_t i = 0; i < str.size(); i++)
        if (str[i] <= ' ')
            str[i] = '_';
    return str;
}

inline std::vector<std::string> split(const std::string &s, const std::string &separator) {
    std::vector<std::string> ans;
    size_t sep_len = separator.length();
    size_t found = s.find(separator);
    size_t previous = 0;
    while (found != std::string::npos) {
        ans.push_back(s.substr(previous, found - previous));
        previous = found + sep_len;
        found = s.find(separator, found + sep_len);
    }
    ans.push_back(s.substr(previous));
    return ans;
}

inline void string_to_khashes(std::string &s, size_t k, std::vector<size_t> &ans) {
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
std::string seqan_string_to_string(const T &s) {
    char buffer[length(s) + 1];
    strcpy(buffer, seqan::String<char, seqan::CStyle>(s));
    return std::string((char *) buffer);
}