#pragma once

class IgGeneProbabilityModel {

};

class NongenomicInsertionModel {

};

class PalindromeDeletionModel {

};

class HCProbabilityRecombinationModel {
    IgGeneProbabilityModel v_gene_probality_model_;
    IgGeneProbabilityModel d_gene_probality_model_;
    IgGeneProbabilityModel j_gene_probality_model_;
    NongenomicInsertionModel nongenomic_insertion_model_;
    PalindromeDeletionModel palindrome_deletion_model_;

public:
    void ExtractFromFile();

};