# Bash
Helpful [one-liners]((https://github.com/stephenturner/oneliners)) by Stephen Turner

## Test
Print read lengths and their counts in .fastq file. Originally from [Biostars](https://www.biostars.org/p/72433/)
```Bash
awk 'if (NR%4 == 2) {len[len($0)]++} END {for (l in len) {print l, len[l]}}' [input.fastq]
```
Parse & print first column, followed by the remaining columns. Used to modify Salmon output for DESeq2
```Bash
awk '{split($1,a,"|"); print a[1], substr($0,length($1)+1);}' quant.sf | tr -d " "
```
Sort all but header (Below is for when header consists of only 1 line). Originally from [StackOverflow](https://stackoverflow.com/questions/14562423/is-there-a-way-to-ignore-header-lines-in-a-unix-sort)
```Bash
# Parenthesis collects stdout from commands within it
(head -n1 [filename] && tail -n+2 [filename] | sort) > [outname]
```
# Bioinformatics
* [Coordinate Systems - 1- vs. 0-based](https://www.biostars.org/p/84686/)
* [Variant call format (VCF)](https://samtools.github.io/hts-specs/)
   * Be sure to check the matching specification consistent with your `.vcf` file
* [Mutation annotation format (MAF)](https://docs.gdc.cancer.gov/Data/File_Formats/MAF_Format/#introduction)