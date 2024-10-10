"""
This part of the workflow prepares sequences for constructing the phylogenetic tree.

REQUIRED INPUTS:

    metadata    = data/metadata.tsv
    sequences   = data/sequences.fasta
    reference   = ../shared/reference.fasta

OUTPUTS:

    prepared_sequences = results/prepared_sequences.fasta

This part of the workflow usually includes the following steps:

    - augur index
    - augur filter
    - augur align
    - augur mask

See Augur's usage docs for these commands for more details.
"""

if not config['ingest_url_prefix'].endswith("/"):
    config['ingest_url_prefix']+="/"

rule download_metadata:
    output:
        metadata = "data/metadata.tsv"
    params:
        address = config['ingest_url_prefix'] + 'metadata.tsv.zst'
    shell:
        r"""
        curl -fsSL --compressed {params.address:q} |
        zstd -d -c > {output.metadata}
        """

rule download_sequences_for_segment:
    output:
        sequences = "data/{segment}/sequences.fasta"
    params:
        address = lambda w: f"{config['ingest_url_prefix']}{w.segment}/sequences.fasta.zst"
    shell:
        r"""
        curl -fsSL --compressed {params.address:q} |
        zstd -d -c > {output.sequences}
        """


rule filter:
    """
    Filtering to
      - {params.sequences_per_group} sequence(s) per {params.group_by!s}
      - excluding strains in {input.exclude}
      - excluding strains matching {params.exclude_where}
    """
    input:
        sequences = "data/{segment}/sequences.fasta",
        metadata = "data/metadata.tsv",
        exclude = config['filter']['exclude']
    output:
        sequences = "results/{segment}/filtered.fasta"
    params:
        strain_id_field = config["strain_id_field"],
        min_length = lambda w: config['filter']['min_length'][w.segment],
        exclude = config['filter']['exclude'],
        exclude_where = config['filter']['exclude_where']
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output {output.sequences} \
            --min-length {params.min_length} \
            --exclude {input.exclude} \
            --exclude-where "{params.exclude_where}"
        """


rule align:
    """
    Aligning sequences to {input.reference}
      - filling gaps with N
    """
    input:
        sequences = "results/{segment}/filtered.fasta",
        reference = "defaults/oropouche_{segment}.gb"
    output:
        alignment = "results/{segment}/aligned.fasta"
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --fill-gaps \
            --remove-reference
        """
