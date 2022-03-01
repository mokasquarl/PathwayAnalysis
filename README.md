# Definitive Desired Diagnosis (d^3)

```
           _/  _/_/_/
      _/_/_/        _/
   _/    _/    _/_/
  _/    _/        _/
   _/_/_/  _/_/_/
```
from KEGG mapped NGS variants ~go mokas!

## All encompassing NGS pipeline in Python go mokas!~

```
      _/    _/  _/_/_/_/  _/_/_/_/_/  _/    _/   
     _/  _/    _/            _/      _/    _/    
    _/_/      _/_/_/        _/      _/    _/     
   _/  _/    _/            _/      _/    _/      
  _/    _/  _/_/_/_/      _/        _/_/  
```

[1] Genome Alignment - From FASTQ short-reads to BAM files
  - Using either BWA or Bowtie algortihms

[2] Process/Recalibrate - Adjust assembeled BAM files according to best practices
  - Header adjustments, mate-pair fixes, etc

[3] Variant Discvoery, Annotation & Filtration

Variant-related Moldules:
        
        
        [1] Discvoery
        [2] Annotation
        [3] Filtration
        [4] Disease Liklihood
        
Ketu Can utilize two different Alignment Algorithms v0.6.28, BWA & Bowtie.
Then preform:

        [1] Remove Duplicate Reads
        [2] Local Realignment
        [3] Base Quality Recalibration
        
Ketu can Detect Variants via UnifiedGenotyper & mPileup, followed by Variant Annotations & Functional Annotations.
And finally help create custom filters for Genomic Variants with 4 levels of stringencies. Before moving onto pathway analysis.

<img width="1680" src="https://user-images.githubusercontent.com/61995876/156232790-25d57a67-010b-4892-8224-96c480bf8cb4.JPG">
