**IgReC** (Ig Repertoire Constructor) is a tool for accurate reconstruction of full-length adaptive immune repertoires from Rep-seq data. IgReC corrects sequencing and PCR errors using sophisticated clustering technique and reports a repertoire as a set of reads corresponding to the identical BCRs/TCRs. Also IgReC performs alignment against immunoglobulin and TCR germline segments and report alignment details.


App Inputs
====
IgReC takes overlapped paired-end or single Illumina Rep-seq reads as an input.   

App Outputs
====
IgReC filters contaminated reads, corrects sequencing and amplification errors and joins reads corresponding to the identical BCRs/TCRs in clusters. The constructed repertoire is a set of clusters characterized by the sequence and the multiplicity. IgReC provides the following information about the constructed repertoire:


* **final_repertoire.fa:** sequences of the constructed clusters (cluster multiplicities are encoded in header lines);
* **final_repertoire.rcm:** correspondence between original reads and constructed clusters.

### Additional output
* **final_repertoire_large.fa:** the sequences corresponding to the large clusters (by default, clusters with multiplicity >= 5);
* **super_reads.fa:** overrepresented original reads;
* **vjfinder:** the directory containing information about alignment against germline segments: cleaned reads, filtered reads and alignment details.

Stages
====
IgReC pipeline consists of the following steps:
* **VJ Finder:** cleaning the _input reads_ using the alignment against the germline genes;
* **HG Constructor:** construction of the Hamming graph on the cleaned reads;
* **Dense Subgraph Finder:** finding the dense subgraphs in the constructed Hamming graph;
* **Consensus Constructor**: construction of the _clusters_ based on the graph decomposition computed at the previous stage.


Known App Limitations
====
Please note that IgReC is designed for the full-length repertoire construction and it expects that the input reads cover whole V(D)J region (~360--400 nt).
Thus, the input reads should be at least 250 nt x2.

**This application is for research use only.**


Link to a Demo Dataset
====
* https://basespace.illumina.com/s/TaAXOszSJglh


Support
====
For support email to [igtools_support@googlegroups.com](mailto:igtools_support@googlegroups.com)


References
====
* **Citations:**
  * Yana Safonova, Stefano Bonissone, Eugene Kurpilyansky, Ekaterina Starostina, Alla Lapidus, Jeremy Stinson, Laura DePalatis, Wendy Sandoval, Jennie Lill and Pavel A. Pevzner. _IgRepertoireConstructor: a novel algorithm for antibody repertoire construction and immunoproteogenomics analysis._ Bioinformatics. 2015 Jun 15; 31(12): i53-i61
http://bioinformatics.oxfordjournals.org/content/31/12/i53.long
  * Alexander Shlemov, Sergey Bankevich, Andrey Bzikadze, Yana Safonova. _New algorithmic challenges of adaptive immune repertoire construction._  RECOMB-SEQ 2016. Los Angeles, USA
http://arxiv.org/abs/1604.02193
* **Web page:** http://yana-safonova.github.io/ig_repertoire_constructor
* **GitHub repository:** https://github.com/yana-safonova/ig_repertoire_constructor
