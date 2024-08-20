import pandas as pd
import json
import argparse
from augur.io import read_metadata

def replace_name_recursive(node, lookup, saveoldcolumn):
    if node["name"] in lookup:
        if saveoldcolumn == "accession":
            node["node_attrs"][saveoldcolumn] = node["name"]
        elif saveoldcolumn == "genbank_accession":
            node["node_attrs"][saveoldcolumn] = {}
            node["node_attrs"][saveoldcolumn]["value"] = node["name"]
        else:
            node["node_attrs"][saveoldcolumn] = node["name"]

    if node["name"] in lookup:
        node["name"] = lookup[node["name"]]

    if "children" in node:
        for child in node["children"]:
            replace_name_recursive(child, lookup, saveoldcolumn)

if __name__=="__main__":
    parser = argparse.ArgumentParser(
        description="Swaps out the strain names in the Auspice JSON with the final strain name",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--input-auspice-json', type=str, required=True, help="input auspice_json")
    parser.add_argument('--metadata', type=str, required=True, help="input data")
    parser.add_argument('--metadata-id-columns', nargs="+", help="names of possible metadata columns containing identifier information, ordered by priority. Only one ID column will be inferred.")
    parser.add_argument('--display-strain-name', type=str, required=True, help="field to use as strain name in auspice")
    parser.add_argument('--output', type=str, metavar="JSON", required=True, help="output Auspice JSON")
    args = parser.parse_args()

    metadata = read_metadata(args.metadata, id_columns=args.metadata_id_columns)
    name_lookup = {}
    for ri, row in metadata.iterrows():
        strain_id = row.name
        name_lookup[strain_id] = args.display_strain_name if pd.isna(row[args.display_strain_name]) else row[args.display_strain_name]

    with open(args.input_auspice_json, 'r') as fh:
        data = json.load(fh)

    replace_name_recursive(data['tree'], name_lookup, "genbank_accession")

    with open(args.output, 'w') as fh:
        json.dump(data, fh)
