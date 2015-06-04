#pragma once

#include "include_me.hpp"

char get_complementary(char nucl){
	if(nucl == 'A' || nucl == 'a')
		return 'T';
	if(nucl == 'T' || nucl == 't')
		return 'A';
	if(nucl == 'C' || nucl == 'c')
		return 'G';
	if(nucl == 'G' || nucl == 'g')
		return 'C';
	if(nucl == 'N' || nucl == 'n')
		return 'A';
	cout << "Char " << nucl << " is not from nucleotide alphabet" << endl;
	assert(false);
	return 'A';
}

string reverse_complementary(string seq) {
	string rc_seq = seq;
	for(size_t i = 0; i < seq.size(); i++)
		rc_seq[seq.size() - i - 1] = get_complementary(seq[i]);
	return rc_seq;
}

size_t HammingDistance(string s1, string s2) {
	assert(s1.size() == s2.size());
	size_t dist = 0;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			dist++;
	return dist;
}

set<size_t> DifferentPositions(string s1, string s2) {
	assert(s1.size() == s2.size());
	set<size_t> pos;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			pos.insert(i);
	return pos;
}

string random_correction(string s1, string s2) {
	assert(s1.size() == s2.size());
	string s = s1;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			s[i] = s2[i];
	return s;
}

char GetRangomNucleotide() {
	static char nucls[4] = { 'A', 'C', 'G', 'T'};
	return nucls[rand() % 4];
}

string GetPalindrom(size_t half_size) {
	string str;
	for(size_t i = 0; i < half_size; i++)
		str = str + GetRangomNucleotide();
	for(size_t i = 0; i < half_size; i++)
		str = str + get_complementary(str[half_size - i - 1]);
	return str;
}

char GetAnotherRandomNucleotide(char nucl) {
	char res;
	do
		res = GetRangomNucleotide();
	while(res == nucl);
	return res;
}
