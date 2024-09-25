
rule group_segments:
    input:
        metadata="data/metadata_merged.tsv",
        resolutions=config["grouping"]["resolutions"],
    output:
        metadata="results/metadata.tsv"
    params:
        common_strain_fields = config["grouping"]["common_strain_fields"],
        segments = segments,
    shell:
        r"""
        python3 scripts/group_segments.py \
            --metadata {input.metadata} \
            --common-strain-fields {params.common_strain_fields} \
            --segments {params.segments} \
            --resolutions {input.resolutions} \
            --output-metadata {output.metadata}
        """

rule subset_sequences_by_segment:
    input:
        metadata = "results/metadata.tsv",
        sequences = "data/sequences.fasta",
    output:
        kv_map = temp("data/kv-map_{segment}.tsv"),
        sequences = "results/{segment}/sequences.fasta",
    params:
        columns = lambda w: f"accession_{w.segment},strain",
        filter_exp = lambda w: f"len($accession_{w.segment})>0",
        drop_key = "__DROP__",
    shell:
        r"""
        cat results/metadata.tsv \
            | csvtk cut -t -f {params.columns} \
            | csvtk filter2 -t -U -f {params.filter_exp:q} \
            > {output.kv_map} && \
        seqkit replace \
            -p "(.*)" --replacement "{{kv}}" --kv-file {output.kv_map} -m {params.drop_key} \
            {input.sequences} \
            | seqkit grep -v -r -p '^{params.drop_key}$' \
            > {output.sequences}
        """
