
### find_sequence_occurrence.py 


```
usage: find_sequence_occurrence.py [-h] [-e EVALUE] [--seq_db SEQ_DB]
                                   [--blast_path BLAST_PATH] [-o OUTPUT]
                                   input

DESCRIPTION: Run BLAST and select positive results from a large set of fastq
or fasta sequences to find existence of specific sequences

positional arguments:
  input                 the input fastq/fasta file containing reads sequences

optional arguments:
  -h, --help            show this help message and exit
  -e EVALUE, --evalue EVALUE
                        cutoff evalues, only choose alignment with evalue less
                        than this cutoffs (default: 1e-5).
  --seq_db SEQ_DB       BLAST database of specific sequences, or fasta
                        sequences
  --blast_path BLAST_PATH
                        path for the local blast program
  -o OUTPUT, --output OUTPUT
                        output file contains sequences with both linker
                        alignment
None
```
