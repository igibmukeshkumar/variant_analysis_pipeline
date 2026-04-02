#!/bin/bash

# Check if input file and data-type are given
if [ -z "$1" ] || [ "$2" != "--data-type" ] || [ -z "$3" ]; then
  echo "Usage: $0 <input_json_file> --data-type <joint|exome|genome>"
  exit 1
fi

input_file="$1"
data_type="$3"
prefix="$data_type"

# Strip .json safely and define output
output_file="$(basename "$input_file" .json)_sorted.txt"

# Extract variant and population-level information using jq
jq -r --arg dt "$data_type" --arg pfx "$prefix" '
# Extract variants and map to flat structure
[
  .data.gene.variants[] |
  {
    variant_id,
    gene_symbol,
    pos,
    rsids: (.rsids | join(",") // "NA"),
    transcript_id,
    transcript_version,
    hgvs,
    hgvsc,
    hgvsp,
    consequence,

    # Main cohort-level counts
    ("\($pfx)_ac"): (.[$dt].ac // 0),
    ("\($pfx)_ac_hemizygote_count"): (.[$dt].ac_hemizygote_count // 0),
    ("\($pfx)_ac_homozygote_count"): (.[$dt].ac_homozygote_count // 0),
    ("\($pfx)_an"): (.[$dt].an // 0),

    # Subpopulation-level counts
    ("\($pfx)_afr_ac"): (.[$dt].populations[]? | select(.id == "afr") | .ac // 0),
    ("\($pfx)_afr_an"): (.[$dt].populations[]? | select(.id == "afr") | .an // 0),
    ("\($pfx)_afr_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "afr") | .ac_hemizygote_count // 0),
    ("\($pfx)_afr_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "afr") | .ac_homozygote_count // 0),

    ("\($pfx)_amr_ac"): (.[$dt].populations[]? | select(.id == "amr") | .ac // 0),
    ("\($pfx)_amr_an"): (.[$dt].populations[]? | select(.id == "amr") | .an // 0),
    ("\($pfx)_amr_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "amr") | .ac_hemizygote_count // 0),
    ("\($pfx)_amr_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "amr") | .ac_homozygote_count // 0),

    ("\($pfx)_asj_ac"): (.[$dt].populations[]? | select(.id == "asj") | .ac // 0),
    ("\($pfx)_asj_an"): (.[$dt].populations[]? | select(.id == "asj") | .an // 0),
    ("\($pfx)_asj_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "asj") | .ac_hemizygote_count // 0),
    ("\($pfx)_asj_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "asj") | .ac_homozygote_count // 0),

    ("\($pfx)_eas_ac"): (.[$dt].populations[]? | select(.id == "eas") | .ac // 0),
    ("\($pfx)_eas_an"): (.[$dt].populations[]? | select(.id == "eas") | .an // 0),
    ("\($pfx)_eas_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "eas") | .ac_hemizygote_count // 0),
    ("\($pfx)_eas_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "eas") | .ac_homozygote_count // 0),

    ("\($pfx)_fin_ac"): (.[$dt].populations[]? | select(.id == "fin") | .ac // 0),
    ("\($pfx)_fin_an"): (.[$dt].populations[]? | select(.id == "fin") | .an // 0),
    ("\($pfx)_fin_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "fin") | .ac_hemizygote_count // 0),
    ("\($pfx)_fin_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "fin") | .ac_homozygote_count // 0),

    ("\($pfx)_mid_ac"): (.[$dt].populations[]? | select(.id == "mid") | .ac // 0),
    ("\($pfx)_mid_an"): (.[$dt].populations[]? | select(.id == "mid") | .an // 0),
    ("\($pfx)_mid_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "mid") | .ac_hemizygote_count // 0),
    ("\($pfx)_mid_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "mid") | .ac_homozygote_count // 0),

    ("\($pfx)_nfe_ac"): (.[$dt].populations[]? | select(.id == "nfe") | .ac // 0),
    ("\($pfx)_nfe_an"): (.[$dt].populations[]? | select(.id == "nfe") | .an // 0),
    ("\($pfx)_nfe_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "nfe") | .ac_hemizygote_count // 0),
    ("\($pfx)_nfe_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "nfe") | .ac_homozygote_count // 0),

    ("\($pfx)_sas_ac"): (.[$dt].populations[]? | select(.id == "sas") | .ac // 0),
    ("\($pfx)_sas_an"): (.[$dt].populations[]? | select(.id == "sas") | .an // 0),
    ("\($pfx)_sas_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "sas") | .ac_hemizygote_count // 0),
    ("\($pfx)_sas_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "sas") | .ac_homozygote_count // 0),

    ("\($pfx)_remaining_ac"): (.[$dt].populations[]? | select(.id == "remaining") | .ac // 0),
    ("\($pfx)_remaining_an"): (.[$dt].populations[]? | select(.id == "remaining") | .an // 0),
    ("\($pfx)_remaining_ac_hemizygote_count"): (.[$dt].populations[]? | select(.id == "remaining") | .ac_hemizygote_count // 0),
    ("\($pfx)_remaining_ac_homozygote_count"): (.[$dt].populations[]? | select(.id == "remaining") | .ac_homozygote_count // 0)
  }
] 
| (.[0] | keys_unsorted) as $keys 
| $keys, map([.[ $keys[] ]])[] 
| @tsv
' "$input_file" > "$output_file"

echo "✅ Output written to: $output_file"
