# Automated PCR Primer Design

This is a [Marimo](https://marimo.io) notebook, that automates the process of designing PCR primers
for a given SNP. See the [corresponding blog post](https://harry.vangberg.name/posts/automated-pcr-primer-design).

 The process is as follows:

1. You enter a SNP ID.
2. The sequence flanking the SNP is downloaded from [Ensembl](https://www.ensembl.org).
3. [Primer3](https://primer3.org) is used to design primers for the SNP.
4. The primers are BLASTed against the human genome.

## Prerequisites

You need to have the following installed:

- [Marimo](https://marimo.io)
- [Primer3](https://primer3.org)
- [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

## Usage

1. Download the human reference genome in the same directory as the notebook – this might take a while:

   ```bash
   update_blastdb.pl --decompress human_genome
   ```

2. Open the notebook in Marimo:

   ```bash
   marimo edit notebook.py
   ```
