### draw_barcode.py
Draw and colorize Hamming (or Levenshtein) graph for FASTQ file data
Example:
> src/ashtools/draw_barcode.py --sgraphs age_ig_s3_R12_raw_cropped/ -i s3_bcodes/CGACTATGAAGGCT.fastq -F 1.pdf -G 1.dot --tau 100000 --span --delta 100

### read_barcode_splitter.py
Split barcoded datatset usign barcodes
Example:
> src/ashtools/read_barcode_splitter.py -s test_dataset/age_ig_s3_R12_raw.cropped.fastq -o s3_bcodes -v -f

### align.py
Align reads in file
Example:
> src/ashtools/align.py -v -d 20 s3_bcodes/CGACTATGAAGGCT.fastq /dev/null
