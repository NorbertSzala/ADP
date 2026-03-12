--------------------------------------------------------------------------------------
                            SEQUENCE vs. STRUCTURE SEARCH
--------------------------------------------------------------------------------------

Motivation:
Structure is [three to ten times] more conserved than sequence
For more details read: 
https://doi.org/10.1002/prot.22458

--------------------------------------------------------------------------------------

Our main goal is to check how many Mycoplasmoides genitalium G37 homologs can be found 
in CP014940 

--------------------------------------------------------------------------------------

I. SEQUENCE SEARCH
To make structure comparisions you can use many programs like:
- BLAST (very slow & gold standard),
- hhblits (more accurate, slower & gold standard),
- jackhmmer (more accurate, slower)
- Diamond (more accurate, slower),
- usearch (more accurate, slower),
...

Today, we will use classic BLAST

First, copy the relevant fasta files from previous lab and filter sequences below 401 aa and rename accordingly:
python3 ../lab2/fasta_splitter.py -i CP014940.fas -o CP014940_400aa_ -m 400 -s 500
Should create: CP014940_400aa_000-296.fas

python3 ../lab2/fasta_splitter.py -i Mycgen.fas -o Mycgen_400aa_ -m 400 -s 500
Should create: Mycgen_400aa_000-432.fas

For BLAST, we can use:
a) website (non-bioinformatic way)
a) API using biopython [https://biopython.org/docs/latest/Tutorial/index.html], 
b) NCBI [https://blast.ncbi.nlm.nih.gov/doc/blast-help/developerinfo.html],
d) install it locally

(b) & (c) are extremly slow and not very flexible, thus they are not recommended, we will go for option (d).

Download latest version of BLAST+
https://blast.ncbi.nlm.nih.gov/doc/blast-help/downloadblastdata.html
https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/

Or from our mirror (faster and recommended):
https://www.mimuw.edu.pl/~lukaskoz/teaching/adp/labs/lab3/ncbi-blast-2.17.0+-x64-linux.tar.gz

1) Create BLAST formated db from CP014940_400aa_000-296.fas

./ncbi-blast-2.17.0+/bin/makeblastdb -in CP014940_400aa_000-296.fas -dbtype prot -out CP014940

2) Run 'blastp' Mycgen_400aa_000-432.fas against CP014940 database with tabular output (-outfmt 6 and -outfmt 7)

https://www.ncbi.nlm.nih.gov/books/NBK279684/table/appendices.T.options_common_to_all_blast/
https://www.metagenomics.wiki/tools/blast/blastn-output-format-6

The command:
./ncbi-blast-2.17.0+/bin/blastp -query Mycgen_400aa_000-432.fas -db CP014940 -out Mycgen2CP014940_mapping_fmt6.txt -outfmt 6 -evalue 1e-10 -num_threads 4

As you can see, frequently there are multiple hits for the same query sequence (expected), thus use some more verbose format (-outfmt -7)

./ncbi-blast-2.17.0+/bin/blastp -query Mycgen_400aa_000-432.fas -db CP014940 -out Mycgen2CP014940_mapping_fmt7.txt -outfmt 7 -evalue 1e-10 -num_threads 4

We can easily estimate how many queries could not be found at all in CP014940. 
The easiest is to count occurence of '# 0 hits found' in Mycgen2CP014940_mapping_fmt7.txt

Mycgen_400aa_000-432.fas --> 432 sequences
Mycgen2CP014940_mapping_fmt7.txt --> yield 432-253=179 good hits

Obviously, in real life analysis we would need to investigate in detail alignemnts 
(e.g. -outfmt 5 and parse XML format), but for simplicity we will stop here for now.

--------------------------------------------------------------------------------------

II. STRUCTURE SEARCH

To make structure comparisions you can use many programs like:
- DALI (very slow & gold standard), 
- TM-align (slow & gold standard),
- GTalign (fast),
- Foldseek (fast),
- Reseek (fast),
- mTM-align2 (fast),
...

Today, we will use TM-align
https://aideepmed.com/TM-align/

1) create folders 'CP014940' and 'Mycgen' and copy all ESM structures done in previous lab

Note: we already filtered out sequences above 400aa, thus we can start search

TMalign require list of pdb structures, you can easily do it with:
ls ./Mycgen/*.pdb | xargs -n 1 basename -s .pdb > list_mycgen.txt
ls ./CP014940/*.pdb | xargs -n 1 basename -s .pdb > list_cp.txt

And then finally run the TMalign:
time ./TMalign -dir1 ./Mycgen/ list_mycgen.txt -dir2 ./CP014940/ list_cp.txt -suffix .pdb -outfmt 2 -fast >  Mycgen2CP014940_TMalign.tsv

Now we need to do some filtering according TM-score:
awk '$3 > 0.5 || $4 > 0.5' Mycgen2CP014940_TMalign.tsv | cut -f1 | sort | uniq -c | sort -nr

awk '$3 > 0.5 || $4 > 0.5' Mycgen2CP014940_TMalign.tsv | cut -f1 | sort | uniq -c | sort -nr|wc

Obviously, in real life analysis we would need to investigate in detail structural alignements, but for simplicity we will stop here for now.

2) Task: Do the comparision using one (up to 80% grade) or two (up to 100% grade) fast method(s): Foldseek, Reseek, GTalign, mTM-align2.

Compare the runtime and the number of reliable hits. Map how many proteins were detected
by both methods. Write what was not found by blast, but found by structure method. 

Make alignment and 3D superposition of two such cases.

--------------------------------------------------------------------------------------

Homework: Send summary of the runs with Foldseek, Reseek, GTalign, or/and mTM-align2.
Here you can attach tabular, summary files (if available). Describe how many proteins were detected in both. 

What are unique hits for BLAST and structure method(s) (are they false positives or false negatives).

Compare the number of hits with sequence and structure methods.

The homework should be sent by Monday, 16.03.2026, via email to [email here]
with the subject: 'ADP26_lab03_hw_Surname_Name'
and the attachment: 'ADP26_lab03_hw_Surname_Name.7z' 

--------------------------------------------------------------------------------------
