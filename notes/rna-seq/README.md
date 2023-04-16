# RSeQC
A suite of tools for working with RNA-Seq data. Useful when little is known about the sequencing data (e.g., quality control, strandedness, etc...).

Check their [SourceForge](https://sourceforge.net/projects/rseqc/files/BED/Human_Homo_sapiens/) page for required dependencies.

- Check for rRNA contamination by aligning reads to rRNA sequences using `split_bam.py`
   - [split_bam.py](https://rseqc.sourceforge.net/#spilt-bam-py)
- "Guess" whether RNA-Seq data is stranded or not by using `infer_experiment.py`
   - [infer_experiment.py](http://rseqc.sourceforge.net/#infer-experiment-py)