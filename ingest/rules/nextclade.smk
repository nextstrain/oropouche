"""
This part of the workflow handles running Nextclade on the curated metadata
and sequences to split the sequences into L, M, and S segments.

REQUIRED INPUTS:

    metadata     = data/subset_metadata.tsv
    all_metadata = results/all/metadata.tsv
    sequences    = results/all/sequences.fasta

OUTPUTS:

    metadata        = results/{segment}/metadata.tsv
    sequences       = results/{segment}/sequences.fasta

See Nextclade docs for more details on usage, inputs, and outputs if you would
like to customize the rules:
https://docs.nextstrain.org/projects/nextclade/page/user/nextclade-cli.html
"""

rule run_nextclade_to_identify_segment:
    input:
        sequences = "data/sequences.fasta",
        segment_reference = config["nextclade"]["segment_reference"],
    output:
        nextclade = temp("data/nextclade_{segment}.tsv"),
    params:
        min_seed_cover = config["nextclade"]["min_seed_cover"],
    shell:
        r"""
        nextclade run \
            --input-ref {input.segment_reference} \
            --output-tsv {output.nextclade} \
            --min-seed-cover {params.min_seed_cover} \
            --silent \
            {input.sequences}
        """

rule parse_nextclade_tsv:
    input:
        nextclade = "data/nextclade_{segment}.tsv",
    output:
        summary = "data/nextclade_{segment}_summary.tsv",
    params:
        nextclade_cols = 'seqName,qc.overallStatus',
        new_cols = lambda w: f'accession,qc_{w.segment}',
        mutate_exp = lambda w: f'len($qc_{w.segment})>0 ? "1" : "0"',
        segment_col = lambda w: f'segment_{w.segment}',
    shell:
        r"""
        csvtk cut -t -H -f {params.nextclade_cols:q} {input.nextclade:q} \
            | csvtk rename -t -f {params.nextclade_cols:q} -n {params.new_cols:q} \
            | csvtk mutate2 -t -n {params.segment_col:q} --at 2 -e {params.mutate_exp:q} \
            > {output.summary:q}

        echo "Nextclade aligned $(( $(cat {output.summary} | csvtk grep -t -f {params.segment_col} -p '1' -U | wc -l) ))/$(( $(wc -l < {input.nextclade}) -1 )) sequences to segment {wildcards.segment}"
        """


rule merge_metadata:
    input:
        strain="data/strain-names.tsv",
        main="data/metadata_subset.tsv",
        segments=expand("data/nextclade_{segment}_summary.tsv", segment=segments),
    output:
        metadata="data/metadata_merged.tsv",
    params:
        # augur merge requires NAME=FILEPATH argments, so we transform the inputs here:
        segments = lambda w,input: " ".join([f"s_{idx}={s}" for idx,s in enumerate(input.segments)])
    shell:
        r"""
        augur merge \
            --metadata strains={input.strain} main={input.main} {params.segments} \
            --metadata-id-columns accession \
            --output-metadata {output.metadata}
        """
