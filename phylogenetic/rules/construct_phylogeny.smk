"""
This part of the workflow constructs the phylogenetic tree.

REQUIRED INPUTS:

    metadata            = data/metadata.tsv
    prepared_sequences  = results/prepared_sequences.fasta

OUTPUTS:

    tree            = results/tree.nwk
    branch_lengths  = results/branch_lengths.json

This part of the workflow usually includes the following steps:

    - augur tree
    - augur refine

See Augur's usage docs for these commands for more details.
"""
rule tree:
    """Building tree"""
    input:
        alignment = "results/{segment}/aligned.fasta"
    output:
        tree = "results/{segment}/tree_raw.nwk"
    params:
        method = "iqtree"
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method}
        """

rule refine:
    """
    Refining tree
      - estimate timetree
      - use {params.coalescent} coalescent timescale
      - estimate {params.date_inference} node dates
      - fix clock rate at {params.clock_rate}
    """
    input:
        tree = "results/{segment}/tree_raw.nwk",
        alignment = "results/{segment}/aligned.fasta",
        metadata = "data/{segment}/metadata.tsv",
    output:
        tree = "results/{segment}/tree.nwk",
        node_data = "results/{segment}/branch_lengths.json"
    params:
        strain_id_field = config["strain_id_field"],
        coalescent = config['refine']['coalescent'],
        date_inference = config['refine']['date_inference'],
        # clock_rate = config['refine']['clock_rate'],
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --coalescent {params.coalescent} \
            --date-confidence \
            --date-inference {params.date_inference}
        """