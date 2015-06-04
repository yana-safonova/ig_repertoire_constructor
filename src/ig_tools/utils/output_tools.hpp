#include "include_me.hpp"

const size_t HEADER_LEN = 80;
const size_t BEF_HEADER_LEN = 20;

string GetHeader(string header, char del) {
	stringstream ss;
	for(size_t i = 0; i < BEF_HEADER_LEN; i++)
		ss << del;
	ss << " " << header << " ";
	size_t aft_header_len = HEADER_LEN - BEF_HEADER_LEN - header.size() - 2;
	for(size_t i = 0; i < aft_header_len; i++)
		ss << del;
	return ss.str();
}

string GetH1(string header) {
	return GetHeader(header, '=');
}

string GetH2(string header) {
	return GetHeader(header, '-');
}

string GetEnd() {
	stringstream ss;
	for(size_t i = 0; i < HEADER_LEN; i++)
		ss << "=";
	return ss.str();
}
