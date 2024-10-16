# SNP Distance Pipeline

This pipeline is designed to compute SNP (Single Nucleotide Polymorphism) distances between samples in a provided genome alignment and a set of novel assemblies. It utilizes several bioinformatics tools such as Snippy, Bcftools, BEDTools, and a custom Python script to process the assemblies, create consensus sequences, apply masks, and ultimately compute SNP distances.

Install with snk ()

```bash
snk install wytamma/snp-dist-pipe
```

# Usage

```bash
❯ snp-dist-pipe run -h
╭─ Workflow Configuration ────────────────────────────────────────────────────────────╮
│ *  --assemblies        PATH  [default: None] [required]                             │
│ *  --alignment         PATH  [default: None] [required]                             │
│    --extension         TEXT  [default: .fa]                                         │
│    --mask              PATH  [default: resources/masks/marin-conservative.bed]      │
│    --reference         PATH  [default: resources/H37Rv.fasta]                       │
╰─────────────────────────────────────────────────────────────────────────────────────╯
```

```bash
snp-dist-pipe run \
    --cores 16 \
    --alignment path/to/alignment.fa \
    --assemblies path/to/assemblies/dir
```
