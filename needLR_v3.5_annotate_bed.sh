#!/bin/bash

##Run from inside needLR_local

##Make a bed file with 10 columns (create empty columns for less than 10)
#Chr\tStart_Pos\tEnd_Pos\tCol_4\tCol_5\tCol_6\tCol_7\tCol_8\tCol_9\tCol_10

##Load environment with 
  #bedtools v2.31.1 


while getopts g:b: flag
do
    case "${flag}" in
        g) REF_GENOME=${OPTARG};;
        b) QUERY_BEDS=${OPTARG};;
    esac
done

##This starts the loop for all of the files in QUERY_BEDS
while IFS= read -r QUERY_FILE_PATH; do

##This will be the basename for the output files
QUERY_BED_NAME=$(basename "$QUERY_FILE_PATH" .bed)_needLR_v3.5_bed
QUERY_SAMPLE_ID=$(basename "$QUERY_FILE_PATH" .bed)

##This is the output directory specific to each query VCF
UNSOLVED_DIR="needLR_output/$QUERY_BED_NAME"
mkdir $UNSOLVED_DIR

##These are the bed files used for annotating the SVs
GENES="backend_files/bed_files/PROTEIN_CODING_GENE_gencode.v45.annotation_5kb_slop.bed"
UTR="backend_files/bed_files/ENSEMBL_CANONICAL_PROTEIN_CODING_UTR_gencode.v45.annotation.bed"
CDS="backend_files/bed_files/ENSEMBL_CANONICAL_EXON_in_PROTEIN_CODING_GENE_NO_UTR_gencode.v45.annotation.bed"
OMIM_GENES="backend_files/bed_files/20251213_OMIM_genemap2_cut.bed"
GENCC="backend_files/bed_files/gencc-submissions_20250510_cleaned.tsv"
HPO="backend_files/bed_files/HPO_genes_to_phenotype_cut.txt"
PLI="backend_files/bed_files/gnomAD_v4.1_constraint_CURATED.txt"
OREGANNO="backend_files/bed_files/ORegAnno.bed"
CENTROMERES="backend_files/bed_files/hg38_centromeres_endtoend.bed"
PERICENTROMERES="backend_files/bed_files/hg38_pericentromeres_5Mb.bed"
TELOMERES="backend_files/bed_files/hg38_telomeres_5Mb.bed"
VAMOS="backend_files/bed_files/from_SG_STR_VNTR_original_motifs.set148_cut.bed"
SEGDUPS="backend_files/bed_files/SEGDUPS_GIAB_v3.3.bed"
HOMOPOLYMERS="backend_files/bed_files/gus_homemade_homopolymers_ge50_nodecoy.bed"
REPEAT_MASKER="backend_files/bed_files/UCSC_hg38_Repeats_RepeatMasker.bed"
DEFRABB_HICONF="backend_files/bed_files/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed"
GAPS="backend_files/bed_files/UCSC_hg38_Mapping_and_Sequencing_Gap.bed"
GENOME_FILE="backend_files/bed_files/hg38.chrom.sizes.txt"


########################################
# Annotate BED (10-column version)
########################################

awk '($5 > -10000000 && $5 < 10000000) { print }' $QUERY_FILE_PATH > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp0.txt

# Annotate with gene info (+/-1kb)
bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp0.txt -b $GENES > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp1.txt

# Sort by gene name
sort -k14,14 "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp1.txt > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp2.txt

# Join with OMIM annotations
join -t $'\t' -a 1 -1 14 -2 4 -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.14,2.5' \
   "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp2.txt "$OMIM_GENES" > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp3.txt

# Join with GenCC annotations
join -t $'\t' -a 1 -1 11 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2' \
   "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp3.txt "$GENCC" > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp4.txt

# Join with HPO annotations
join -t $'\t' -a 1 -1 11 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2,2.3' \
   "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp4.txt "$HPO" > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp5.txt

# Join with pLI scores
join -t $'\t' -a 1 -1 11 -2 1 -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.2' \
   "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp5.txt "$PLI" > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp6.txt

# Consolidate values in columns 11–16
awk '
BEGIN { FS = OFS = "\t" }
{
   key = $1
   for (i = 2; i <= 10; i++) key = key FS $i

   a[key][$11] = 1
   b[key][$12] = 1
   c[key][$13] = 1
   d[key][$14] = 1
   e[key][$15] = 1
   f[key][$16] = 1
}
END {
   for (k in a) {
       split(k, fields, FS)
       for (i = 1; i <= 10; i++)
           printf "%s%s", fields[i], (i < 10 ? OFS : "")

       a_vals = ""; for (v in a[k]) a_vals = a_vals ? a_vals ", " v : v
       b_vals = ""; for (v in b[k]) b_vals = b_vals ? b_vals ", " v : v
       c_vals = ""; for (v in c[k]) c_vals = c_vals ? c_vals ", " v : v
       d_vals = ""; for (v in d[k]) d_vals = d_vals ? d_vals ", " v : v
       e_vals = ""; for (v in e[k]) e_vals = e_vals ? e_vals ", " v : v
       f_vals = ""; for (v in f[k]) f_vals = f_vals ? f_vals ", " v : v

       printf "%s%s%s%s%s%s\n", OFS a_vals, OFS b_vals, OFS c_vals, OFS d_vals, OFS e_vals, OFS f_vals
   }
}' "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp6.txt > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp7.txt

# Replace empty fields with "."
awk 'BEGIN {FS=OFS="\t"} {for (i=1; i<=NF; i++) if ($i=="") $i="."} 1' \
   "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp7.txt > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp8.txt


# Annotate
bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp8.txt \
  -b "$UTR" \
| cut -f1-16,18 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($17 == -1 || $17 == ".") $17=".";
         else $17="UTR";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp9.txt


bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp9.txt \
  -b "$CDS" \
| cut -f1-17,19 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($18 == -1 || $18 == ".") $18=".";
         else $18="CDS";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp10.txt

bedtools sort -i "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp10.txt > "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp10.1.txt

bedtools map \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp10.1.txt \
  -b "$OREGANNO" \
  -c 4 \
  -o distinct \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp11.txt

bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp11.txt \
  -b "$CENTROMERES" \
| cut -f1-19,21 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($20 == -1 || $20 == ".") $20=".";
         else $20="centromeric";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp12.txt

bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp12.txt \
  -b "$PERICENTROMERES" \
| cut -f1-20,22 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($21 == -1 || $21 == ".") $21=".";
         else $21="pericentromeric";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp13.txt

bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp13.txt \
  -b "$TELOMERES" \
| cut -f1-21,23 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($22 == -1 || $22 == ".") $22=".";
         else $22="telomeric";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp14.txt

###Vamos

awk 'BEGIN{FS=OFS="\t"}
     $4!="INS" && $4!="DEL" && $4!="DUP" {
       print $0, ".";
     }' \
  "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp14.txt \
> "$UNSOLVED_DIR"/novamos.temp

bedtools slop -b 50 -i "$VAMOS" -g "$GENOME_FILE" \
| bedtools intersect -wa -wb -loj \
    -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp14.txt \
    -b - \
| awk -F'\t' '$4=="INS"' \
> "$UNSOLVED_DIR"/ins.temp

bedtools slop -b 50 -i "$VAMOS" -g "$GENOME_FILE" \
| bedtools intersect -wa -wb -loj \
    -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp14.txt \
    -b - \
| awk -F'\t' '$4=="DUP"' \
> "$UNSOLVED_DIR"/dup.temp

bedtools intersect -wa -wb -loj -f 0.5 \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp14.txt \
  -b "$VAMOS" \
| awk -F'\t' '$4=="DEL"' \
> "$UNSOLVED_DIR"/del.temp

cat "$UNSOLVED_DIR"/ins.temp \
    "$UNSOLVED_DIR"/del.temp \
    "$UNSOLVED_DIR"/dup.temp \
| cut -f1-22,24 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($23 == -1 || $23 == ".") $23=".";
         else $23="vamos";
         print;
       }' \
> "$UNSOLVED_DIR"/vamos_annot.temp


cat "$UNSOLVED_DIR"/novamos.temp \
    "$UNSOLVED_DIR"/vamos_annot.temp \
| sort | uniq \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp15.txt

###


bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp15.txt \
  -b "$SEGDUPS" \
| cut -f1-23,25 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($24 == -1 || $24 == ".") $24=".";
         else $24="segdups";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp16.txt


bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp16.txt \
  -b "$REPEAT_MASKER" \
| cut -f1-24,26 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($25 == -1 || $25 == ".") $25=".";
         else $25="repeat";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp17.txt


bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp17.txt \
  -b "$GAPS" \
| cut -f1-25,27 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($26 == -1 || $26 == ".") $26=".";
         else $26="Gap";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp18.txt


bedtools intersect -wa -wb -loj -f 1.0 \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp18.txt \
  -b "$HOMOPOLYMERS" \
| cut -f1-26,28 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($27 == -1 || $27 == ".") $27=".";
         else $27="HP>=50bp";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp19.txt


bedtools intersect -wa -wb -loj -f 1.0 \
  -a "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp19.txt \
  -b "$DEFRABB_HICONF" \
| cut -f1-27,29 \
| awk 'BEGIN{FS=OFS="\t"}
       {
         if ($28 == -1 || $28 == ".") $28=".";
         else $28="hiconf";
         print;
       }' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_temp20.txt


# Add headers
echo -e "Chr\tStart_Pos\tEnd_Pos\tSV_type\tSV_length\tCol_6\tCol_7\tCol_8\tCol_9\tCol_10\tGenes\tOMIM\tGenCC\tHPO\tHPO\tpLI\tUTR\tCDS\tORegAnno\tCentromeric\tPericentromeric\tTelomeric\tVAMOS\tSegdup\tRepeat\tGap\tHomopolymer\tHiConf" \
> "$UNSOLVED_DIR"/"$QUERY_BED_NAME"_RESULTS.txt

# Clean up and finalize
sed -i 's/?//g' "${UNSOLVED_DIR}/${QUERY_BED_NAME}_temp20.txt"
sed -i -e :a -e '/^[[:space:]]*$/{$d;N;ba}' "${UNSOLVED_DIR}/${QUERY_BED_NAME}_temp20.txt"
cat "${UNSOLVED_DIR}/${QUERY_BED_NAME}_temp20.txt" >> "${UNSOLVED_DIR}/${QUERY_BED_NAME}_RESULTS.txt"

# Remove temporary files
rm ${UNSOLVED_DIR}/*temp*

##This facilitates and ends the loop for the input VCF list
done < "$QUERY_BEDS"
