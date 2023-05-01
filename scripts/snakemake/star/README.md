# Example Snakemake script for running STAR

Directory structure:
```
star/
|-- env
|-- input_list
|-- reference
|-- README.md
```

## How to run

`snakemake --config curr_file=1 -s gtex_star_snakefile.py`

* `env` - Contains the configuration setting
* `input_list` - List of files to be processed
* `reference` - Reference files required to run STAR (not included in the repo)
    * E.g., `GRCh38.d1.vd1.fa` and `star-2.7.5c_GRCh38.d1.vd1_gencode.v36`

## Note
May need to specify number of jobs and dynamically set required memory as needed.