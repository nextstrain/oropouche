"""
This part of the workflow creates additonal annotations for the phylogenetic tree.

REQUIRED INPUTS:

    metadata            = data/metadata.tsv
    prepared_sequences  = results/prepared_sequences.fasta
    tree                = results/tree.nwk

OUTPUTS:

    node_data = results/*.json

    There are no required outputs for this part of the workflow as it depends
    on which annotations are created. All outputs are expected to be node data
    JSON files that can be fed into `augur export`.

    See Nextstrain's data format docs for more details on node data JSONs:
    https://docs.nextstrain.org/page/reference/data-formats.html

This part of the workflow usually includes the following steps:

    - augur traits
    - augur ancestral
    - augur translate
    - augur clades

See Augur's usage docs for these commands for more details.

Custom node data files can also be produced by build-specific scripts in addition
to the ones produced by Augur commands.
"""
rule ancestral:
    """Reconstructing ancestral sequences and mutations"""
    input:
        tree = "results/{segment}/tree.nwk",
        alignment = "results/{segment}/aligned.fasta",
    output:
        node_data = "results/{segment}/nt_muts.json"
    params:
        inference = "joint"
    log:
        "logs/{segment}/ancestral.txt",
    benchmark:
        "benchmarks/{segment}/ancestral.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference}
        """

rule translate:
    """Translating amino acid sequences"""
    input:
        tree = "results/{segment}/tree.nwk",
        node_data = "results/{segment}/nt_muts.json",
        reference = "defaults/oropouche_{segment}.gb"
    output:
        node_data = "results/{segment}/aa_muts.json"
    log:
        "logs/{segment}/translate.txt",
    benchmark:
        "benchmarks/{segment}/translate.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output-node-data {output.node_data}
        """

rule traits:
    """Inferring ancestral traits for {params.columns!s}"""
    input:
        tree = "results/{segment}/tree.nwk",
        metadata = "data/metadata.tsv",
    output:
        node_data = "results/{segment}/traits.json",
    params:
        strain_id_field = config["strain_id_field"],
        columns = config['traits']['columns']
    log:
        "logs/{segment}/traits.txt",
    benchmark:
        "benchmarks/{segment}/traits.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-node-data {output.node_data} \
            --columns {params.columns} \
            --confidence
        """
