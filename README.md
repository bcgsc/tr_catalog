# A tandem repeat (TR) catalog generated from high-quality long-read human genome assemblies
This repository keeps the analysis scripts that were used to generated the TR catalog from public diploid long-read human genome assemblies from the following data soucres:
1. [Human Pangenome Reference Consortium (HPRC)](https://humanpangenome.org/)
2. [Human Genome Structural Variation Consortium (HGSVC2)](https://www.internationalgenome.org/human-genome-structural-variation-consortium/)
3. [1000G ONT Sequencing Consortium](https://millerlaboratory.com/1000G-ONT.html)

## Workflow
![workflow](1b_300_cropped.png)

## Catalog
[v1](https://zenodo.org/records/11522276)
- haplotype names separated by semi-colons are shown in first header line preceded by '#'
- column descriptions:

| Column | Description |
| ------ | ----------- |
| chrom | chromosome |
| start | start coordinate |
| end | end coordinate |
| motif | consensus repeat motif |
| copy_numbers | copy numbers in haplotypes separated by semi-colons ('-' for missing genotypes) |
| sizes | sizes (bp) in haplotypes separated by semi-colons  ('-' for missing genotypes) |
| motifs | motifs in haplotypes separated by semi-colons  ('-' for missing genotypes) |
| max_change | maximum change (of all haplotypes) in size (bp) substracted from reference genome size |
| num_samples | number of samples with genotype |
| num_calls | number of haplotypes with genotype |
| motif_frequency | number of haplotypes associated with each motif observed e.g. CAG(10);CAA(2) |
| feature | gene element overlapped. Format: gene\|transcript\|<element>, where <element> = exon#\|intron#\|utr5\|utr3\|cds\|promoter\|exon_bound (exon boundary) |
