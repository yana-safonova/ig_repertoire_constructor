1. Task --- statistical analysis of sequences. Somatic mutations etc. 
2. Motivation:
  * Simulating sequences (and hence antibodies repertoires).
  * Clonal trees.
  * Comparison of repertoires.
  * Improvements of IgRepertoireConstructor.
3. Article about T-cells.
  * VDJ recombinations.
  * Model with over9k params.
  * Questions: is this model really okay? No overfitting? Are the results statistically significant? Are the similar facts true about B-cells.
4. Several mutations: ... including clavage. Sometimes we get palindromes.
5. Palindromes. Simplest H0. Biological vs Random. 2 Illustrations. Similarity with T-cell paper.
6. The closest goal is to clasterize gens. 
  * Attempt --- Probability of Clavage. 
  * We want to understand that those probabilities are consistant on different data. Age datasets. 
  * ScatterPlots for V gens. 
  * Because lack of effective sample size we dont use z-test on proportions and only calculate Pearson correlations. 
    Those correlations are significant according on Permutation test.
7. The further goals: find the factor to clusterize V gens effectively. The same about D and J gens. To try to find out the correlations between and distributions of different events: pallindromes, clavage etc.
