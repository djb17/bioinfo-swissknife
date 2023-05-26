# Sample GATK workflows

Directory structure:
```
gatk/
|-- mutect2_paired_gdc_snakefile.py
|-- mutect2_tonly_gdc_snakefile.py
|-- README.md
```

Adapted from Genomic Data Commons (GDC) best practices 
* `mutect2_paired_gdc_snakefile.py` - Somatic calling using MuTect2
* `mutect2_tonly_gdc_snakefile.py` - Tumor-only Somatic calling (per chromosome) using MuTect2 using panel of normals