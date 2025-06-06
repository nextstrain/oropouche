"""
This part of the workflow collects the phylogenetic tree and annotations to
export a Nextstrain dataset.

REQUIRED INPUTS:

    metadata        = data/metadata.tsv
    tree            = results/tree.nwk
    branch_lengths  = results/branch_lengths.json
    node_data       = results/*.json

OUTPUTS:

    auspice_json = auspice/${build_name}.json

    There are optional sidecar JSON files that can be exported as part of the dataset.
    See Nextstrain's data format docs for more details on sidecar files:
    https://docs.nextstrain.org/page/reference/data-formats.html

This part of the workflow usually includes the following steps:

    - augur export v2
    - augur frequencies

See Augur's usage docs for these commands for more details.
"""
rule colors:
    input:
        color_schemes = "defaults/color_schemes.tsv",
        color_orderings = "defaults/color_orderings.tsv",
        metadata = "data/metadata.tsv",
    output:
        colors = "results/{segment}/colors.tsv"
    log:
        "logs/{segment}/colors.txt",
    benchmark:
        "benchmarks/{segment}/colors.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/assign-colors.py \
            --color-schemes {input.color_schemes} \
            --ordering {input.color_orderings} \
            --metadata {input.metadata} \
            --output {output.colors}
        """

rule export:
    """Exporting data files for for auspice"""
    input:
        tree = "results/{segment}/tree.nwk",
        metadata = "data/metadata.tsv",
        branch_lengths = "results/{segment}/branch_lengths.json",
        traits = "results/{segment}/traits.json",
        nt_muts = "results/{segment}/nt_muts.json",
        aa_muts = "results/{segment}/aa_muts.json",
        colors = "results/{segment}/colors.tsv",
        description = config['export']['description'],
        auspice_config = config['export']['auspice_config'],
    output:
        auspice = "auspice/oropouche_{segment}.json",
    params:
        strain_id_field = config["strain_id_field"],
        metadata_columns = lambda w: [name.format(segment=w.segment) for name in config["export"]["segment_metadata_columns"]]
    log:
        "logs/{segment}/export.txt",
    benchmark:
        "benchmarks/{segment}/export.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id_field} \
            --node-data {input.branch_lengths} {input.traits} {input.nt_muts} {input.aa_muts} \
            --description {input.description} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --metadata-columns {params.metadata_columns} \
            --output {output.auspice} \
            --include-root-sequence-inline
        """
