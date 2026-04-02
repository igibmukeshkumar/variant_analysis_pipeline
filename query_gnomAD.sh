#!/bin/bash

# Function to display usage
usage() {
  echo "Usage: $0 --gl <gene_list.txt> --data-type <exome|genome|joint>"
  echo "  --gl        : Path to the gene list file (e.g., genes.txt)"
  echo "  --data-type : Choose data type (exome, genome, joint)"
  echo "  --help      : Display this help message"
  exit 1
}

# Parse command line options
while [[ "$1" != "" ]]; do
  case "$1" in
    --gl)
      shift
      GENE_LIST_FILE="$1"
      ;;
    --data-type)
      shift
      DATA_TYPE="$1"
      ;;
    --help)
      usage
      ;;
    *)
      echo "Unknown option: $1"
      usage
      ;;
  esac
  shift
done

# Validate inputs
if [ -z "$GENE_LIST_FILE" ]; then
  echo "Error: Gene list file is required (--gl)"
  usage
fi

if [ -z "$DATA_TYPE" ] || [[ "$DATA_TYPE" != "exome" && "$DATA_TYPE" != "genome" && "$DATA_TYPE" != "joint" ]]; then
  echo "Error: Invalid data type. Please specify exome, genome, or joint (--data-type)"
  usage
fi

# Check if gene list file exists
if [ ! -f "$GENE_LIST_FILE" ]; then
  echo "File not found: $GENE_LIST_FILE"
  exit 1
fi

# Loop through each gene in the file
while read -r GENE; do
  # Skip empty lines
  [ -z "$GENE" ] && continue

  echo "Querying gnomAD for gene: $GENE"
  OUTFILE="${GENE}_gnomad_${DATA_TYPE}.json"

  # Set FIELDS depending on DATA_TYPE
  if [ "$DATA_TYPE" == "joint" ]; then
    FIELDS=$(cat <<EOF
      ${DATA_TYPE} {
        ac
        an
        hemizygote_count
        homozygote_count
        populations {
          id
          ac
          an
          hemizygote_count
          homozygote_count
        }
      }
EOF
)
  else
    FIELDS=$(cat <<EOF
      ${DATA_TYPE} {
        ac
        an
        af
        ac_hemi
        ac_hom
        populations {
          id
          ac
          an
          ac_hemi
          ac_hom
        }
      }
EOF
)
  fi

  # Send GraphQL query to gnomAD
  curl https://gnomad.broadinstitute.org/api \
    -H 'Content-Type: application/graphql; charset=utf-8' \
    --data-binary @- <<EOF > "$OUTFILE"
query VariantsInGene {
  gene(gene_symbol: "$GENE", reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      gene_symbol
      pos
      rsids
      transcript_id
      transcript_version
      hgvs
      hgvsc
      hgvsp
      consequence
      $FIELDS
    }
  }
}
EOF

  echo "Saved output to $OUTFILE"

done < "$GENE_LIST_FILE"
