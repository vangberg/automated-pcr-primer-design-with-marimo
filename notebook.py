# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "biopython==1.85",
#     "marimo",
#     "openai==1.81.0",
#     "pandas==2.2.3",
#     "primer3-py==2.2.0",
#     "primers==0.5.10",
# ]
# ///

import marimo

__generated_with = "0.14.7"
app = marimo.App(width="medium", app_title="PCR Primer Design")

with app.setup:
    import marimo as mo
    from Bio import Entrez, SeqIO, Blast
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import Alignment
    import primer3
    import subprocess

    # Set the email for NCBI Entrez
    Entrez.email = "harry@vangberg.name"


@app.cell
def _():
    mo.md(
        r"""
    ## SNP

    Enter the RS number for the SNP you want to design primers for, e.g. "429358" (APOE) or "334" (sickle-cell anemia):
    """
    )
    return


@app.cell
def _():
    snp_id = mo.ui.text(label="SNP RS number:", value="429358")
    snp_id
    return (snp_id,)


@app.cell
def _():
    mo.md(
        r"""The sequence flanking the SNP (+/- 1000 bp) is fetched from the NCBI Nucleotide database:"""
    )
    return


@app.cell
def _(snp_id):
    mo.stop(snp_id.value == "", "Enter a SNP RS number to fetch the sequence")

    # Fetch the record from the NCBI SNP database.
    _handle = Entrez.esummary(db="snp", id=snp_id.value, retmode="xml")
    _record = Entrez.read(_handle)

    # Extract the accession number
    acc = _record["DocumentSummarySet"]["DocumentSummary"][0]["ACC"]

    # Extract the chromosome position
    chrpos = int(
        _record["DocumentSummarySet"]["DocumentSummary"][0]["CHRPOS_SORT"]
    )

    # Fetch the sequence around the SNP position
    _handle = Entrez.efetch(
        db="nucleotide",
        id=acc,
        seq_start=chrpos - 1000,
        seq_stop=chrpos + 1000,
        rettype="fasta",
        retmode="text",
    )
    seqrec = SeqIO.read(_handle, "fasta")
    mo.md(str(seqrec.seq)).style({"word-break": "break-all"})
    return (seqrec,)


@app.cell
def _():
    mo.md(
        r"""
    ## Primers

    [Primer3](https://primer3.org/) is used to design primer pairs:
    """
    )
    return


@app.cell
def _():
    product_size_range = mo.ui.range_slider(
        label="Product size range",
        show_value=True,
        start=0,
        stop=1000,
        value=[250, 500],
    )
    product_size_range
    return (product_size_range,)


@app.cell
def _(product_size_range, seqrec):
    # See the Primer3 manual for a full list of options ("tags" in Primer3 parlance): https://primer3.org/manual.html
    primers = primer3.bindings.design_primers(
        seq_args={
            "SEQUENCE_ID": "PRIMER",
            "SEQUENCE_TEMPLATE": str(seqrec.seq),
            "SEQUENCE_TARGET": [1000, 1],
        },
        global_args={
            "PRIMER_OPT_SIZE": 20,
            "PRIMER_MIN_SIZE": 18,
            "PRIMER_MAX_SIZE": 25,
            "PRIMER_OPT_TM": 60.0,
            "PRIMER_MIN_TM": 57.0,
            "PRIMER_MAX_TM": 63.0,
            "PRIMER_MIN_GC": 20.0,
            "PRIMER_MAX_GC": 80.0,
            "PRIMER_MAX_POLY_X": 100,
            "PRIMER_PRODUCT_SIZE_RANGE": [product_size_range.value],
        },
    )
    return (primers,)


@app.cell
def _(primers):
    import pandas as pd

    # Create a DataFrame from the primers dictionary
    _data = [
        {
            "Index": _i,
            "Pair": f"Pair {_i + 1}",
            "Length (bp)": primers["PRIMER_PAIR"][_i]["PRODUCT_SIZE"],
            "Left TM": round(primers["PRIMER_RIGHT"][_i]["TM"], 1),
            "Left GC%": primers["PRIMER_LEFT"][_i]["GC_PERCENT"],
            "Left Length": primers["PRIMER_LEFT"][_i]["COORDS"][1],
            "Right TM": round(primers["PRIMER_RIGHT"][_i]["TM"], 1),
            "Right GC%": primers["PRIMER_RIGHT"][_i]["GC_PERCENT"],
            "Right Length": primers["PRIMER_RIGHT"][_i]["COORDS"][1],
        }
        for _i in range(5)
    ]

    primers_table = mo.ui.table(
        data=_data,
        selection="single",
        initial_selection=[0],
    )
    primers_table
    return (primers_table,)


@app.cell
def _():
    mo.md(r"""### Primer sequence""")
    return


@app.cell
def _(primers, primers_table):
    mo.stop(len(primers_table.value) == 0, "Select a primer pair to show sequence")
    _idx = primers_table.value[0]["Index"]
    _left_primer = primers["PRIMER_LEFT"][_idx]["SEQUENCE"]
    _right_primer = primers["PRIMER_RIGHT"][_idx]["SEQUENCE"]
    mo.md(
        f"""
        ```
        Left:  {_left_primer}
        Right: {_right_primer}
        ```
        """
    )
    return


@app.cell
def _():
    mo.md(r"""## BLAST""")
    return


@app.cell
def _(primers):
    # Primer sequences are stored on disk, so `blastn` can read them.
    primer_seqs = []

    for _i in range(5):
        left_seq = SeqRecord(
            Seq(primers["PRIMER_LEFT"][_i]["SEQUENCE"]),
            id=f"Primer Pair {_i + 1} Left",
        )
        right_seq = SeqRecord(
            Seq(primers["PRIMER_RIGHT"][_i]["SEQUENCE"]),
            id=f"Primer Pair {_i + 1} Right",
        )

        primer_seqs.append(left_seq)
        primer_seqs.append(right_seq)

    # This is a Marimo hack. We need to run `blastn` whenever we save the `fasta` file. When we store the file name
    # in a variable, Marimo will automatically re-run all cells that read the `primers_fasta_file` variable whenever
    # this cell runs.
    primers_fasta_file = "primer_seqs.fasta"
    SeqIO.write(primer_seqs, primers_fasta_file, "fasta")
    None
    return (primers_fasta_file,)


@app.cell
def _(primers_fasta_file):
    blast_results_file = "blast_results.xml"

    with mo.status.spinner(title="Running blastn…"):
        _cmd = [
            "blastn",
            "-query",
            primers_fasta_file,
            "-db",
            "GCF_000001405.39_top_level",
            "-outfmt",
            "5",
            "-out",
            blast_results_file,
            # "-outfmt",
            # "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-task",
            "blastn-short",
            "-perc_identity",
            "95",
            "-strand",
            "both",
            "-max_target_seqs",
            "10",
            "-num_threads",
            "4",  # Remember to convert to string if you use a variable here
        ]

        subprocess.run(
            _cmd,
            check=True,
        )
    return (blast_results_file,)


@app.cell
def _(blast_results_file, primers_fasta_file):
    " ".join(
        [
            "blastn",
            "-query",
            primers_fasta_file,
            "-db",
            "GCF_000001405.39_top_level",
            "-outfmt",
            "5",
            "-out",
            blast_results_file,
            # "-outfmt",
            # "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            "-task",
            "blastn-short",
            "-perc_identity",
            "95",
            "-strand",
            "both",
            "-max_target_seqs",
            "10",
            "-num_threads",
            "4",  # Remember to convert to string if you use a variable here
        ]
    )
    return


@app.cell
def _(blast_results_file):
    # Parse BLAST results
    _records = Blast.parse(blast_results_file)
    blast_records = []
    for _record in _records:
        blast_records.append(_record)
    return (blast_records,)


@app.cell
def _():
    mo.md(r"""Select a query to show hits:""")
    return


@app.cell
def _(blast_records):
    _data = [
        {
            "Index": _idx,
            "Primer Pair": _record.query.description,
            # "Identity (%)": r[0].identity,
        }
        for _idx, _record in enumerate(blast_records)
    ]

    blast_table = mo.ui.table(
        data=_data, selection="single", initial_selection=[0]
    )
    blast_table
    return (blast_table,)


@app.cell
def _():
    mo.md(
        r"""
    ### Hits

    Select a hit to show alignments:
    """
    )
    return


@app.cell
def _(blast_records, blast_table):
    mo.stop(len(blast_table.value) == 0, "Select a primer to show hits")
    blast_record = blast_records[blast_table.value[0]["Index"]]
    _data = [
        {
            "Index": _idx,
            "Target": _hit.target.description,
            "Best e-value": _hit[0].annotations["evalue"],
        }
        for _idx, _hit in enumerate(blast_record)
    ]
    hits_table = mo.ui.table(_data, selection="single", initial_selection=[0])
    hits_table
    return blast_record, hits_table


@app.cell
def _():
    mo.md(r"""### Alignments""")
    return


@app.cell
def _(blast_record, hits_table):
    mo.stop(len(hits_table.value) == 0, "Select a hit to show alignments")
    hit = blast_record[hits_table.value[0]["Index"]]
    hit
    return


if __name__ == "__main__":
    app.run()
