from genomic_kmers.genomic_kmers import get_genomic_kmers_all_chains
from chains.chains import Chains

kmers_fr, kmers_cdr = get_genomic_kmers_all_chains()
print(len(kmers_fr[Chains.IGL]))
print(len(kmers_cdr[Chains.IGL]))
