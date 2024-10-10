rule copy_example_sequences:
    input:
        sequences="example_data/sequences_{segment}.fasta",
    output:
        sequences="data/{segment}/sequences.fasta",
    shell:
        """
        cp -f {input.sequences} {output.sequences}
        """

rule copy_example_metadata:
    input:
        metadata="example_data/metadata.tsv",
    output:
        metadata="data/metadata.tsv",
    shell:
        """
        cp -f {input.metadata} {output.metadata}
        """

# Add a Snakemake ruleorder directive here if you need to resolve ambiguous rules
# that have the same output as the copy_example_data rule.

ruleorder: copy_example_sequences > download_sequences_for_segment
ruleorder: copy_example_metadata > download_metadata
