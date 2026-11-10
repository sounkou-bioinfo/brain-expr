#!/usr/bin/env bash
# Downloads the inputs of extract/mouse_dev.R into <raw_dir> (about 21 GB, wild-type samples only).
set -euo pipefail
cd "${1:?usage: fetch_mouse_dev.sh <raw_dir>}"
geo() { # accession file -> geo/<accession>/<file>
  local acc=$1 f=$2 p
  p=$(echo "$acc" | sed -E 's/[0-9]{3}$/nnn/')
  mkdir -p "geo/$acc"
  [ -s "geo/$acc/$f" ] || curl -sSfL --retry 5 -o "geo/$acc/$f" "https://ftp.ncbi.nlm.nih.gov/geo/series/$p/$acc/suppl/$f"
}
gsm() { # sample accession file -> geo/<series>/<file>
  local series=$1 gsm=$2 f=$3 p
  p=$(echo "$gsm" | sed -E 's/[0-9]{3}$/nnn/')
  mkdir -p "geo/$series"
  [ -s "geo/$series/$f" ] || curl -sSfL --retry 5 -o "geo/$series/$f" "https://ftp.ncbi.nlm.nih.gov/geo/samples/$p/$gsm/suppl/$f"
}

# Allaway 2021 (GSE164570): Dlx lineage, E13 MGE -> P28 cortex; Mef2c cKO excluded.
for s in GSM5014302_MGE_Dlx6pos GSM5014303_MGE_Dlx6neg GSM5014304_P2_Cortex_Dlxpos \
         GSM5014305_P10_ALM_Dlxpos GSM5014306_P10_V1_Dlxpos GSM5014307_P28_ALM_Dlxpos GSM5014308_P28_V1_Dlxpos; do
  for k in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do gsm GSE164570 "${s%%_*}" "${s}_RNA_${k}"; done
done

# Lee 2022 (GSE167013): E12.5 cortex, MGE, CGE, LGE.
for s in GSM5090774_CTX GSM5090775_MGE GSM5090776_CGE GSM5090777_LGE; do
  for k in barcodes.tsv.gz features.tsv.gz matrix.mtx.gz; do gsm GSE167013 "${s%%_*}" "${s}_${k}"; done
done

# Di Bella 2021 (GSE153162): wild-type cortex E10.5 -> P4; Fezf2 mutants excluded.
for f in GSM5277843_E10_v1_filtered_feature_bc_matrix.h5 GSM4635072_E11_5_filtered_gene_bc_matrices_h5.h5 \
         GSM4635073_E12_5_filtered_gene_bc_matrices_h5.h5 GSM4635074_E13_5_filtered_gene_bc_matrices_h5.h5 \
         GSM4635075_E14_5_filtered_gene_bc_matrices_h5.h5 GSM4635076_E15_5_S1_filtered_gene_bc_matrices_h5.h5 \
         GSM4635077_E16_filtered_gene_bc_matrices_h5.h5 GSM5277844_E17_5_filtered_feature_bc_matrix.h5 \
         GSM4635078_E18_5_S1_filtered_gene_bc_matrices_h5.h5 GSM4635079_E18_S3_filtered_gene_bc_matrices_h5.h5 \
         GSM4635080_P1_S1_filtered_gene_bc_matrices_h5.h5 GSM4635081_P1_S2_filtered_gene_bc_matrices_h5.h5 \
         GSM5277845_P4_filtered_feature_bc_matrix.h5; do
  gsm GSE153162 "${f%%_*}" "$f"
done

# Mayer 2018 (GSE104156): 10x, E13.5 MGE Lhx6+, CGE, E18.5 and P10 cortex and subcortex.
geo GSE104156 GSE104156_digital_expression.csv.gz
# Hochgerner 2018 (GSE104323): dentate gyrus E16.5 -> P132, published cell types.
geo GSE104323 GSE104323_10X_expression_data_V2.tab.gz
geo GSE104323 GSE104323_metadata_barcodes_24185cells.txt.gz

du -sh geo/*

# Yao 2023, Allen Brain Cell Atlas: adult striatum/amygdala (STR) and hippocampus (HPF), with cell metadata.
abc=https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com
mkdir -p abc
for f in WMB-10Xv3-STR WMB-10Xv3-HPF; do
  [ -s "abc/$f-raw.h5ad" ] || curl -sSfL --retry 5 -o "abc/$f-raw.h5ad" "$abc/expression_matrices/WMB-10Xv3/20230630/$f-raw.h5ad"
done
[ -s abc/cell_metadata_with_cluster_annotation.csv ] || curl -sSfL --retry 5 -o abc/cell_metadata_with_cluster_annotation.csv \
  "$abc/metadata/WMB-10X/20231215/views/cell_metadata_with_cluster_annotation.csv"
du -sh abc
