#pragma once

#include <iostream>
#include <set>

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
    std::cout << "Char " << nucl << " is not from nucleotide alphabet" << std::endl;
	VERIFY(false);
	return 'A';
}

std::string reverse_complementary(std::string seq) {
	std::string rc_seq = seq;
	for(size_t i = 0; i < seq.size(); i++)
		rc_seq[seq.size() - i - 1] = get_complementary(seq[i]);
	return rc_seq;
}

size_t HammingDistance(std::string s1, std::string s2) {
	VERIFY(s1.size() == s2.size());
	size_t dist = 0;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			dist++;
	return dist;
}

std::set<size_t> DifferentPositions(std::string s1, std::string s2) {
	VERIFY(s1.size() == s2.size());
    std::set<size_t> pos;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			pos.insert(i);
	return pos;
}

std::string random_correction(std::string s1, std::string s2) {
	VERIFY(s1.size() == s2.size());
	std::string s = s1;
	for(size_t i = 0; i < s1.size(); i++)
		if(s1[i] != s2[i])
			s[i] = s2[i];
	return s;
}

char GetRangomNucleotide() {
	static char nucls[4] = { 'A', 'C', 'G', 'T'};
	return nucls[rand() % 4];
}

std::string GetPalindrom(size_t half_size) {
	std::string str;
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
