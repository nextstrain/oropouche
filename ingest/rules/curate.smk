"""
This part of the workflow handles the curation of data from NCBI

REQUIRED INPUTS:

    ndjson      = data/ncbi.ndjson

OUTPUTS:

    metadata    = data/subset_metadata.tsv
    sequences   = results/all/sequences.fasta

"""


def format_field_map(field_map: dict[str, str]) -> list[str]:
    """
    Format entries to the format expected by `augur curate --field-map`.

    When used in a Snakemake shell block, the list is automatically expanded and
    spaces are handled by quoted interpolation.
    """
    return  [f'{key}={value}' for key, value in field_map.items()]


# This curate pipeline is based on existing pipelines for pathogen repos using NCBI data.
# You may want to add and/or remove steps from the pipeline for custom metadata
# curation for your pathogen. Note that the curate pipeline is streaming NDJSON
# records between scripts, so any custom scripts added to the pipeline should expect
# the input as NDJSON records from stdin and output NDJSON records to stdout.
# The final step of the pipeline should convert the NDJSON records to two
# separate files: a metadata TSV and a sequences FASTA.
rule curate:
    input:
        sequences_ndjson="data/ncbi.ndjson",
        local_geolocation_rules=config["curate"]["local_geolocation_rules"],
        annotations=config["curate"]["annotations"],
    output:
        metadata="data/metadata_curated.tsv",
        sequences="data/sequences.fasta",
    params:
        field_map=format_field_map(config["curate"]["field_map"]),
        strain_regex=config["curate"]["strain_regex"],
        strain_backup_fields=config["curate"]["strain_backup_fields"],
        date_fields=config["curate"]["date_fields"],
        expected_date_formats=config["curate"]["expected_date_formats"],
        genbank_location_field=config["curate"]["genbank_location_field"],
        articles=config["curate"]["titlecase"]["articles"],
        abbreviations=config["curate"]["titlecase"]["abbreviations"],
        titlecase_fields=config["curate"]["titlecase"]["fields"],
        authors_field=config["curate"]["authors_field"],
        authors_default_value=config["curate"]["authors_default_value"],
        abbr_authors_field=config["curate"]["abbr_authors_field"],
        annotations_id=config["curate"]["annotations_id"],
        id_field=config["curate"]["output_id_field"],
        sequence_field=config["curate"]["output_sequence_field"],
    log:
        "logs/curate.txt",
    benchmark:
        "benchmarks/curate.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        cat {input.sequences_ndjson} \
            | augur curate rename \
                --field-map {params.field_map:q} \
            | augur curate normalize-strings \
            | augur curate transform-strain-name \
                --strain-regex {params.strain_regex} \
                --backup-fields {params.strain_backup_fields} \
            | augur curate format-dates \
                --date-fields {params.date_fields} \
                --expected-date-formats {params.expected_date_formats} \
            | augur curate parse-genbank-location \
                --location-field {params.genbank_location_field} \
            | augur curate titlecase \
                --titlecase-fields {params.titlecase_fields} \
                --articles {params.articles} \
                --abbreviations {params.abbreviations} \
            | augur curate abbreviate-authors \
                --authors-field {params.authors_field} \
                --default-value {params.authors_default_value} \
                --abbr-authors-field {params.abbr_authors_field} \
            | augur curate apply-geolocation-rules \
                --geolocation-rules {input.local_geolocation_rules} \
            | augur curate apply-record-annotations \
                --annotations {input.annotations} \
                --id-field {params.annotations_id} \
                --output-metadata {output.metadata} \
                --output-fasta {output.sequences} \
                --output-id-field {params.id_field} \
                --output-seq-field {params.sequence_field}
        """

rule subset_curated_metadata_columns:
    input:
        metadata="data/metadata_curated.tsv",
    output:
        metadata="data/metadata_subset.tsv",
    params:
        metadata_fields=",".join(config["curate"]["metadata_columns"]),
    log:
        "logs/subset_curated_metadata_columns.txt",
    benchmark:
        "benchmarks/subset_curated_metadata_columns.txt"
    shell:
        r"""
        exec &> >(tee {log:q})

        csvtk cut -t -f {params.metadata_fields} \
          {input.metadata} \
        > {output.metadata}
        """
