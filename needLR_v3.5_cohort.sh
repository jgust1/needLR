#!/bin/bash

##Run from inside needLR_local

##Load environment with 
  #truvari v4.2.2 
  #bedtools v2.31.1 
  #bcftools v1.19 

##needLR default Truvari parameters
truvari_k="common"
truvari_s=50
truvari_S=10000000
truvari_r=2000
truvari_p=0
truvari_P=0.2
truvari_O=0.2

##Usage:
##  needLR_v3.5_cohort.sh -g REF_GENOME -c COHORT_VCF -n COHORT_SAMPLE_NAMES -s COHORT_SAMPLES
##
##  -g  Path to reference genome (fasta)
##  -v  Pre-Truvari-merged multi-sample cohort VCF (.vcf.gz)
##  -a  File listing cohort sample names (one per line, in VCF column order)
##  -d  Number of samples in the cohort VCF

while getopts g:v:a:d: flag
do
    case "${flag}" in
        g) REF_GENOME=${OPTARG};;
        v) COHORT_VCF=${OPTARG};;
        a) COHORT_SAMPLE_NAMES=${OPTARG};;
        d) COHORT_SAMPLES=${OPTARG};;
    esac
done

COHORT_PREFIX=$(basename "$COHORT_VCF" .vcf.gz)
QUERY_FILE_NAME="${COHORT_PREFIX}_needLR_v3.5_cohort"

UNSOLVED_DIR="needLR_output/$QUERY_FILE_NAME"
mkdir -p $UNSOLVED_DIR

##Hardcoded 1KGP control set (500 samples)
CTRL_MERGED="backend_files/UWONT_500_sniffles_2.6.2_preproc2_T31_v4.2.2.vcf.gz"
CTRL_SAMPLE_NAMES="backend_files/SAMPLE_order.txt"
CTRL_SAMPLES=500
AFR=163
AMR=67
EAS=95
EUR=72
SAS=103

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

##Precompute haplotype counts
COHORT_HAPS=$(( COHORT_SAMPLES * 2 ))
CTRL_HAPS=$(( CTRL_SAMPLES * 2 ))
AFRHAPS=$(( AFR * 2 ))
AMRHAPS=$(( AMR * 2 ))
EASHAPS=$(( EAS * 2 ))
EURHAPS=$(( EUR * 2 ))
SASHAPS=$(( SAS * 2 ))

##---------------------------------------------------------------------------
## STEP 1: Preprocess cohort VCF and merge with 1KGP control set
##---------------------------------------------------------------------------

bcftools view \
  -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49) && GT!="0/0" && GT!="0|0"' \
  -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM \
  -o "$UNSOLVED_DIR"/preprocessed_temp_"$COHORT_PREFIX".vcf \
  "$COHORT_VCF"

bgzip "$UNSOLVED_DIR"/preprocessed_temp_"$COHORT_PREFIX".vcf
tabix "$UNSOLVED_DIR"/preprocessed_temp_"$COHORT_PREFIX".vcf.gz

##Create merge list (cohort first, then 1KGP — order determines column zones)
echo "$UNSOLVED_DIR"/preprocessed_temp_"$COHORT_PREFIX".vcf.gz > "$UNSOLVED_DIR"/list_temp.txt
echo "$CTRL_MERGED" >> "$UNSOLVED_DIR"/list_temp.txt

bcftools merge -m none --force-samples \
  -l "$UNSOLVED_DIR"/list_temp.txt \
  -Oz -o "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_temp.vcf.gz

tabix "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_temp.vcf.gz

##Truvari collapse (T31)
truvari collapse \
  -i "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_temp.vcf.gz \
  -o "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_temp.vcf \
  -c "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_collapsed_temp.vcf \
  -f $REF_GENOME -k $truvari_k -s $truvari_s -S $truvari_S -r $truvari_r -p $truvari_p -P $truvari_P -O $truvari_O

##---------------------------------------------------------------------------
## STEP 2: Parse VCF
##
## bcftools query column layout (1-indexed):
##   $1-$6  : CHROM POS REF ALT SVLEN SVTYPE
##   $7,$8,$9 .. $(7+(COHORT_SAMPLES-1)*3) : cohort GT/DV/DR triplets
##   $(7+COHORT_SAMPLES*3) ..               : 1KGP GT/DV/DR triplets
##---------------------------------------------------------------------------

bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVLEN\t%INFO/SVTYPE[\t%GT\t%DV\t%DR]\n' \
  "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_temp.vcf \
  > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_temp.txt

##Keep only rows where at least one cohort sample has a non-missing genotype
COHORT_LAST_GT=$(( 7 + (COHORT_SAMPLES - 1) * 3 ))

awk -v first=7 -v last="$COHORT_LAST_GT" -v step=3 'BEGIN{OFS="\t"} {
    for (i = first; i <= last; i += step) {
        if ($i != "./.") { print; next }
    }
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_temp.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_cohort_filtered_temp.txt

##Remove BNDs and SVs >10Mb
grep -v BND "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_cohort_filtered_temp.txt \
| awk '($5 > -10000000 && $5 < 10000000)' \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp1.txt

##---------------------------------------------------------------------------
## STEP 3: Collapse to GT-only columns (drop DV/DR)
##
## temp2 layout (1-indexed):
##   $1-$6                  : CHROM POS REF ALT SVLEN SVTYPE
##   $7 .. $6+COHORT        : cohort GTs (one col per sample)
##   $7+COHORT .. $6+COHORT+CTRL : 1KGP GTs (one col per sample)
##---------------------------------------------------------------------------

COHORT_RAW_GT_START=7
COHORT_RAW_GT_END=$(( 7 + (COHORT_SAMPLES - 1) * 3 ))
CTRL_RAW_GT_START=$(( 7 + COHORT_SAMPLES * 3 ))
CTRL_RAW_GT_END=$(( CTRL_RAW_GT_START + (CTRL_SAMPLES - 1) * 3 ))

awk -v cgs="$COHORT_RAW_GT_START" -v cge="$COHORT_RAW_GT_END" \
    -v kgs="$CTRL_RAW_GT_START"   -v kge="$CTRL_RAW_GT_END" \
'BEGIN{OFS="\t"} {
    printf "%s\t%s\t%s\t%s\t%s\t%s", $1, $2, $3, $4, $5, $6
    for (i = cgs; i <= cge; i += 3) printf "\t%s", $i
    for (i = kgs; i <= kge; i += 3) printf "\t%s", $i
    printf "\n"
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp1.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp2.txt

##Cohort and 1KGP GT column positions in temp2 (reused throughout)
C_GT_START=7
C_GT_END=$(( 6 + COHORT_SAMPLES ))
K_GT_START=$(( 7 + COHORT_SAMPLES ))
K_GT_END=$(( 6 + COHORT_SAMPLES + CTRL_SAMPLES ))

##---------------------------------------------------------------------------
## STEP 4: Build SUPP_VECs (one binary string per zone, appended as new cols)
##---------------------------------------------------------------------------

##Append COHORT_SUPP_VEC
awk -v cs="$C_GT_START" -v ce="$C_GT_END" 'BEGIN{OFS="\t"} {
    sv = ""
    for (i = cs; i <= ce; i++) sv = sv ($i ~ /1/ ? "1" : "0")
    print $0, sv
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp2.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp3.txt

##Append CTRL_SUPP_VEC
awk -v ks="$K_GT_START" -v ke="$K_GT_END" 'BEGIN{OFS="\t"} {
    sv = ""
    for (i = ks; i <= ke; i++) sv = sv ($i ~ /1/ ? "1" : "0")
    print $0, sv
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp3.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp4.txt

COHORT_SUPP_COL=$(( K_GT_END + 1 ))
CTRL_SUPP_COL=$(( K_GT_END + 2 ))

##---------------------------------------------------------------------------
## STEP 5: Rearrange to clean working format
##
## temp5 layout (1-indexed):
##   $1  CHROM
##   $2  POS
##   $3  temp (End_Pos placeholder, filled in later)
##   $4  REF
##   $5  ALT
##   $6  SVLEN
##   $7  SVTYPE
##   $8  COHORT_SUPP_VEC
##   $9  CTRL_SUPP_VEC
##   $10 .. $9+COHORT_SAMPLES       : cohort GTs
##   $10+COHORT_SAMPLES .. $9+COHORT+CTRL : 1KGP GTs
##---------------------------------------------------------------------------

awk -v cs="$C_GT_START" -v ce="$C_GT_END" \
    -v ks="$K_GT_START" -v ke="$K_GT_END" \
    -v csup="$COHORT_SUPP_COL" -v ksup="$CTRL_SUPP_COL" \
'BEGIN{OFS="\t"} {
    printf "%s\t%s\ttemp\t%s\t%s\t%s\t%s\t%s\t%s", \
           $1, $2, $3, $4, $5, $6, $(csup), $(ksup)
    for (i = cs; i <= ce; i++) printf "\t%s", $i
    for (i = ks; i <= ke; i++) printf "\t%s", $i
    printf "\n"
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp4.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp5.txt

##GT column positions in temp5 (used for all frequency calculations)
COH_COL_START=10
COH_COL_END=$(( 9 + COHORT_SAMPLES ))
CTRL_COL_START=$(( 10 + COHORT_SAMPLES ))
CTRL_COL_END=$(( 9 + COHORT_SAMPLES + CTRL_SAMPLES ))

##---------------------------------------------------------------------------
## STEP 6: Calculate cohort-level frequency statistics
##
## Appends 5 columns:
##   Pop_COHORT, Allele_COHORT, GT_coh_homWT, GT_coh_het, GT_coh_homVAR
##---------------------------------------------------------------------------

awk -v cs="$COH_COL_START" -v ce="$COH_COL_END" \
'BEGIN{OFS="\t"} {
    Pop=0; Allele=0; homWT=0; het=0; homVAR=0
    for (i = cs; i <= ce; i++) {
        if ($i ~ /1/) Pop++
        n = split($i, arr, "")
        for (j = 1; j <= n; j++) if (arr[j] == "1") Allele++
        if ($i == "./." || $i == ".|." || $i == "0/0" || $i == "0|0") homWT++
        if ($i == "0/1" || $i == "0|1" || $i == "1/0" || $i == "1|0") het++
        if ($i == "1/1" || $i == "1|1") homVAR++
    }
    print $0, Pop, Allele, homWT, het, homVAR
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp5.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp6.txt

POP_COH_COL=$(( CTRL_COL_END + 1 ))
ALLELE_COH_COL=$(( CTRL_COL_END + 2 ))
GT_COH_HOMWT_COL=$(( CTRL_COL_END + 3 ))
GT_COH_HET_COL=$(( CTRL_COL_END + 4 ))
GT_COH_HOMVAR_COL=$(( CTRL_COL_END + 5 ))

##---------------------------------------------------------------------------
## STEP 7: Calculate 1KGP population-stratified frequency statistics
##
## Population order within CTRL_SUPP_VEC and GT columns:
##   AFR (163) → AMR (67) → EAS (95) → EUR (72) → SAS (103)
##
## Appends 15 columns:
##   Pop_AFR, Pop_AMR, Pop_EAS, Pop_EUR, Pop_SAS, Pop_ALL
##   Allele_AFR, Allele_AMR, Allele_EAS, Allele_EUR, Allele_SAS, Allele_ALL
##   GT_homWT, GT_het, GT_homVAR
##---------------------------------------------------------------------------

##SUPP_VEC substring offsets (1-indexed within the vector string)
TEMPAMR=$(( AFR + 1 ))
TEMPEAS=$(( AFR + AMR + 1 ))
TEMPEUR=$(( AFR + AMR + EAS + 1 ))
TEMPSAS=$(( AFR + AMR + EAS + EUR + 1 ))

##Absolute GT column boundaries for each population
TEMPAFRC=$(( CTRL_COL_START + AFR - 1 ))
TEMPAMRB=$(( CTRL_COL_START + AFR ))
TEMPAMRC=$(( CTRL_COL_START + AFR + AMR - 1 ))
TEMPEASB=$(( CTRL_COL_START + AFR + AMR ))
TEMPEASC=$(( CTRL_COL_START + AFR + AMR + EAS - 1 ))
TEMPEURB=$(( CTRL_COL_START + AFR + AMR + EAS ))
TEMPEURC=$(( CTRL_COL_START + AFR + AMR + EAS + EUR - 1 ))
TEMPSASB=$(( CTRL_COL_START + AFR + AMR + EAS + EUR ))
TEMPSASC=$(( CTRL_COL_START + AFR + AMR + EAS + EUR + SAS - 1 ))

awk -v afr="$AFR" -v amr="$AMR" -v eas="$EAS" -v eur="$EUR" -v sas="$SAS" \
    -v ctrl_samples="$CTRL_SAMPLES" \
    -v ks="$CTRL_COL_START" -v ke="$CTRL_COL_END" \
    -v tempamr="$TEMPAMR" -v tempeas="$TEMPEAS" \
    -v tempeur="$TEMPEUR" -v tempsas="$TEMPSAS" \
    -v tempafrc="$TEMPAFRC" \
    -v tempamrb="$TEMPAMRB" -v tempamrc="$TEMPAMRC" \
    -v tempeasb="$TEMPEASB" -v tempeasc="$TEMPEASC" \
    -v tempeurb="$TEMPEURB" -v tempeurc="$TEMPEURC" \
    -v tempsasb="$TEMPSASB" -v tempsasc="$TEMPSASC" \
'BEGIN{OFS="\t"} {
    ##Pop counts via CTRL_SUPP_VEC (col $9) substrings
    sub_AFR = substr($9, 1, afr);           pop_AFR = gsub("1","",sub_AFR)
    sub_AMR = substr($9, tempamr, amr);     pop_AMR = gsub("1","",sub_AMR)
    sub_EAS = substr($9, tempeas, eas);     pop_EAS = gsub("1","",sub_EAS)
    sub_EUR = substr($9, tempeur, eur);     pop_EUR = gsub("1","",sub_EUR)
    sub_SAS = substr($9, tempsas, sas);     pop_SAS = gsub("1","",sub_SAS)
    sub_ALL = substr($9, 1, ctrl_samples);  pop_ALL = gsub("1","",sub_ALL)

    al_AFR=0; al_AMR=0; al_EAS=0; al_EUR=0; al_SAS=0; al_ALL=0
    GT_homWT=0; GT_het=0; GT_homVAR=0

    for (i=ks;       i<=tempafrc; i++) { n=split($i,a,""); for(j=1;j<=n;j++) if(a[j]=="1") al_AFR++ }
    for (i=tempamrb; i<=tempamrc; i++) { n=split($i,a,""); for(j=1;j<=n;j++) if(a[j]=="1") al_AMR++ }
    for (i=tempeasb; i<=tempeasc; i++) { n=split($i,a,""); for(j=1;j<=n;j++) if(a[j]=="1") al_EAS++ }
    for (i=tempeurb; i<=tempeurc; i++) { n=split($i,a,""); for(j=1;j<=n;j++) if(a[j]=="1") al_EUR++ }
    for (i=tempsasb; i<=tempsasc; i++) { n=split($i,a,""); for(j=1;j<=n;j++) if(a[j]=="1") al_SAS++ }
    for (i=ks;       i<=ke;       i++) { n=split($i,a,""); for(j=1;j<=n;j++) if(a[j]=="1") al_ALL++ }

    for (i=ks; i<=ke; i++) {
        if ($i=="./." || $i==".|." || $i=="0/0" || $i=="0|0") GT_homWT++
        if ($i=="0/1" || $i=="0|1" || $i=="1/0" || $i=="1|0") GT_het++
        if ($i=="1/1" || $i=="1|1") GT_homVAR++
    }

    print $0, pop_AFR, pop_AMR, pop_EAS, pop_EUR, pop_SAS, pop_ALL, \
              al_AFR, al_AMR, al_EAS, al_EUR, al_SAS, al_ALL, \
              GT_homWT, GT_het, GT_homVAR
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp6.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp7.txt

POP_AFR_COL=$(( GT_COH_HOMVAR_COL + 1 ))
POP_AMR_COL=$(( GT_COH_HOMVAR_COL + 2 ))
POP_EAS_COL=$(( GT_COH_HOMVAR_COL + 3 ))
POP_EUR_COL=$(( GT_COH_HOMVAR_COL + 4 ))
POP_SAS_COL=$(( GT_COH_HOMVAR_COL + 5 ))
POP_ALL_COL=$(( GT_COH_HOMVAR_COL + 6 ))
AL_AFR_COL=$(( GT_COH_HOMVAR_COL + 7 ))
AL_AMR_COL=$(( GT_COH_HOMVAR_COL + 8 ))
AL_EAS_COL=$(( GT_COH_HOMVAR_COL + 9 ))
AL_EUR_COL=$(( GT_COH_HOMVAR_COL + 10 ))
AL_SAS_COL=$(( GT_COH_HOMVAR_COL + 11 ))
AL_ALL_COL=$(( GT_COH_HOMVAR_COL + 12 ))
GT_HOMWT_COL=$(( GT_COH_HOMVAR_COL + 13 ))
GT_HET_COL=$(( GT_COH_HOMVAR_COL + 14 ))
GT_HOMVAR_COL=$(( GT_COH_HOMVAR_COL + 15 ))

##---------------------------------------------------------------------------
## STEP 8: Expand SUPP_VEC strings to named sample lists
##---------------------------------------------------------------------------

##Replace $8 (COHORT_SUPP_VEC) with cohort sample names
mapfile -t cohort_names < "$COHORT_SAMPLE_NAMES"

awk -v names="$(printf "%s\n" "${cohort_names[@]}")" \
'BEGIN{ split(names, na, "\n"); OFS="\t" } {
    sv = $8
    if (sv ~ /^0+$/) {
        label = "cohort_unique"
    } else {
        label = ""
        for (i = 1; i <= length(sv); i++)
            if (substr(sv, i, 1) == "1")
                label = (label == "") ? na[i] : label "_" na[i]
    }
    $8 = label; print
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp7.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp8.txt

##Replace $9 (CTRL_SUPP_VEC) with 1KGP sample names
mapfile -t ctrl_names < "$CTRL_SAMPLE_NAMES"

awk -v names="$(printf "%s\n" "${ctrl_names[@]}")" \
'BEGIN{ split(names, na, "\n"); OFS="\t" } {
    sv = $9
    if (sv ~ /^0+$/) {
        label = "query_only"
    } else {
        label = ""
        for (i = 1; i <= length(sv); i++)
            if (substr(sv, i, 1) == "1")
                label = (label == "") ? na[i] : label "_" na[i]
    }
    $9 = label; print
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp8.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp9.txt

##---------------------------------------------------------------------------
## STEP 9: HWE calculation for 1KGP
##---------------------------------------------------------------------------

##p = WT allele frequency
awk -v hwt="$GT_HOMWT_COL" -v hhet="$GT_HET_COL" -v hvar="$GT_HOMVAR_COL" '{
    denom = $(hwt) + $(hhet) + $(hvar)
    p = (denom != 0) ? ((2 * $(hwt)) + $(hhet)) / (2 * denom) : 0
    print $0, p
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp9.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp10.txt

P_COL=$(( GT_HOMVAR_COL + 1 ))

##q = 1-p
awk -v p_col="$P_COL" '{ print $0, (1 - $(p_col)) }' \
"$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp10.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp11.txt

Q_COL=$(( P_COL + 1 ))

##np^2, 2npq, nq^2
awk -v p_col="$P_COL" -v q_col="$Q_COL" \
    -v n="$CTRL_SAMPLES" -v nhaps="$CTRL_HAPS" \
'{ p=$(p_col); q=$(q_col); print $0, (n*p*p), (nhaps*p*q), (n*q*q) }' \
"$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp11.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp12.txt

NP2_COL=$(( Q_COL + 1 ))
TWONPQ_COL=$(( Q_COL + 2 ))
NQ2_COL=$(( Q_COL + 3 ))

##Chi-squared
awk -v hwt="$GT_HOMWT_COL" -v hhet="$GT_HET_COL" -v hvar="$GT_HOMVAR_COL" \
    -v np2="$NP2_COL" -v twonpq="$TWONPQ_COL" -v nq2="$NQ2_COL" \
'{
    if ($(np2) > 0 && $(twonpq) > 0 && $(nq2) > 0) {
        x = ((($(hwt)-$(np2))^2)/$(np2)) + \
            ((($(hhet)-$(twonpq))^2)/$(twonpq)) + \
            ((($(hvar)-$(nq2))^2)/$(nq2))
        print $0, x
    } else {
        print $0, "NULL"
    }
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp12.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp13.txt

CHISQ_COL=$(( NQ2_COL + 1 ))

##HWE label
awk -v c="$CHISQ_COL" 'BEGIN{OFS="\t"} {
    if      ($(c) == "NULL")   $(c) = "HWE_NULL"
    else if ($(c)+0 > 3.84)    $(c) = "HWE_FALSE"
    else                       $(c) = "HWE_TRUE"
    print
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp13.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp14.txt

HWE_COL=$CHISQ_COL

##---------------------------------------------------------------------------
## STEP 10: Condense to 43-column summary (drops raw GT column arrays)
##
## $1-$7   : CHROM POS temp REF ALT SVLEN SVTYPE
## $8       : Cohort_support (named)
## $9       : 1KGP_support (named)
## $10      : Cohort_Pop_Count
## $11      : Cohort_Pop_Freq
## $12      : Cohort_Allele_Count
## $13      : Cohort_Allele_Freq
## $14-$15  : Pop_Count_AFR, Pop_Freq_AFR
## $16-$17  : Pop_Count_AMR, Pop_Freq_AMR
## $18-$19  : Pop_Count_EAS, Pop_Freq_EAS
## $20-$21  : Pop_Count_EUR, Pop_Freq_EUR
## $22-$23  : Pop_Count_SAS, Pop_Freq_SAS
## $24-$25  : Pop_Count_ALL, Pop_Freq_ALL
## $26-$27  : Allele_Count_AFR, Allele_Freq_AFR
## $28-$29  : Allele_Count_AMR, Allele_Freq_AMR
## $30-$31  : Allele_Count_EAS, Allele_Freq_EAS
## $32-$33  : Allele_Count_EUR, Allele_Freq_EUR
## $34-$35  : Allele_Count_SAS, Allele_Freq_SAS
## $36-$37  : Allele_Count_ALL, Allele_Freq_ALL
## $38      : GT_homWT (1KGP)
## $39      : GT_het (1KGP)
## $40      : GT_homVAR (1KGP)
## $41      : p
## $42      : q
## $43      : HWE
##---------------------------------------------------------------------------

awk -v pcc="$POP_COH_COL"   -v acc="$ALLELE_COH_COL" \
    -v coh_n="$COHORT_SAMPLES" -v coh_h="$COHORT_HAPS" \
    -v pafr="$POP_AFR_COL"  -v pamr="$POP_AMR_COL"  \
    -v peas="$POP_EAS_COL"  -v peur="$POP_EUR_COL"  \
    -v psas="$POP_SAS_COL"  -v pall="$POP_ALL_COL"  \
    -v aafr="$AL_AFR_COL"   -v aamr="$AL_AMR_COL"   \
    -v aeas="$AL_EAS_COL"   -v aeur="$AL_EUR_COL"   \
    -v asas="$AL_SAS_COL"   -v aall="$AL_ALL_COL"   \
    -v hwt="$GT_HOMWT_COL"  -v hhet="$GT_HET_COL"   \
    -v hvar="$GT_HOMVAR_COL" \
    -v p_col="$P_COL"       -v q_col="$Q_COL"       \
    -v hwe="$HWE_COL" \
    -v afr="$AFR" -v amr="$AMR" -v eas="$EAS" -v eur="$EUR" -v sas="$SAS" \
    -v ctrl_n="$CTRL_SAMPLES" \
    -v afrhaps="$AFRHAPS" -v amrhaps="$AMRHAPS" -v eashaps="$EASHAPS" \
    -v eurhaps="$EURHAPS" -v sashaps="$SASHAPS" -v ctrlhaps="$CTRL_HAPS" \
'BEGIN{OFS="\t"} {
    print $1,$2,$3,$4,$5,$6,$7,$8,$9,
          $(pcc), ($(pcc)/coh_n),
          $(acc), ($(acc)/coh_h),
          $(pafr), ($(pafr)/afr),
          $(pamr), ($(pamr)/amr),
          $(peas), ($(peas)/eas),
          $(peur), ($(peur)/eur),
          $(psas), ($(psas)/sas),
          $(pall), ($(pall)/ctrl_n),
          $(aafr), ($(aafr)/afrhaps),
          $(aamr), ($(aamr)/amrhaps),
          $(aeas), ($(aeas)/eashaps),
          $(aeur), ($(aeur)/eurhaps),
          $(asas), ($(asas)/sashaps),
          $(aall), ($(aall)/ctrlhaps),
          $(hwt), $(hhet), $(hvar),
          $(p_col), $(q_col), $(hwe)
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp14.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp15.txt

##---------------------------------------------------------------------------
## STEP 11: Add End_Pos (replaces "temp" placeholder in $3)
##---------------------------------------------------------------------------

awk 'BEGIN{OFS="\t"} {
    if ($7 == "INS") $3 = $2 + 1
    else             $3 = $2 + ($6 < 0 ? -1*$6 : $6)
    print
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp15.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp16.txt

##---------------------------------------------------------------------------
## STEP 12: Gene annotation
##
## bedtools intersect with GENES appends 4 bed cols ($44-$47); gene name = $47
## sort -k47,47 then sequential joins with OMIM, GenCC, HPO, pLI
##
## Post-join layout (48 cols):
##   $44=Genes $45=OMIM $46=GenCC $47=HPO $48=pLI
##---------------------------------------------------------------------------

bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp16.txt \
  -b $GENES \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp17.txt

sort -k47,47 "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp17.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp18.txt

##Join OMIM (gene=col47, OMIM gene key=col4, value=col5)
join -t $'\t' -a 1 -1 47 -2 4 \
  -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.47,2.5' \
  "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp18.txt "$OMIM_GENES" \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp19.txt

##Join GenCC (gene=col44, GenCC gene key=col1, value=col2)
join -t $'\t' -a 1 -1 44 -2 1 \
  -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,2.2' \
  "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp19.txt "$GENCC" \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp20.txt

##Join HPO (gene=col44, HPO gene key=col1, value=col3)
join -t $'\t' -a 1 -1 44 -2 1 \
  -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,2.3' \
  "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp20.txt "$HPO" \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp21.txt

##Join pLI (gene=col44, pLI gene key=col1, value=col2)
join -t $'\t' -a 1 -1 44 -2 1 \
  -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.45,1.46,1.47,2.2' \
  "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp21.txt "$PLI" \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp22.txt

##Consolidate multi-gene rows: key=cols 1-43, annotation cols 44-48
awk 'BEGIN{FS=OFS="\t"} {
    key = $1; for (i=2; i<=43; i++) key = key FS $i
    a[key][$44]=1; b[key][$45]=1; c[key][$46]=1; d[key][$47]=1; e[key][$48]=1
}
END {
    for (k in a) {
        split(k, f, FS)
        for (i=1; i<=43; i++) printf "%s%s", f[i], (i<43 ? OFS : "")
        av=""; for(v in a[k]) av = av ? av ", " v : v
        bv=""; for(v in b[k]) bv = bv ? bv ", " v : v
        cv=""; for(v in c[k]) cv = cv ? cv ", " v : v
        dv=""; for(v in d[k]) dv = dv ? dv ", " v : v
        ev=""; for(v in e[k]) ev = ev ? ev ", " v : v
        printf "%s%s%s%s%s%s%s%s%s%s\n", OFS,av,OFS,bv,OFS,cv,OFS,dv,OFS,ev
    }
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp22.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp23.txt

##Replace empty fields with "."
awk 'BEGIN{FS=OFS="\t"} {for(i=1;i<=NF;i++) if($i=="") $i="."; print}' \
"$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp23.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.txt

##---------------------------------------------------------------------------
## STEP 13: Bed annotation pipeline
##
## Starting from temp24 (48 cols: 43 summary + gene/OMIM/GenCC/HPO/pLI)
## Each step appends one annotation column:
##   $49=UTR  $50=CDS  $51=ORegAnno  $52=Centromeric  $53=Pericentromeric
##   $54=Telomeric  $55=Vamos  $56=Segdup  $57=Repeat  $58=Gap
##   $59=Homopolymer  $60=HiConf
##---------------------------------------------------------------------------

bedtools sort -i "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.1.txt

bedtools map \
  -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.1.txt \
  -b "$UTR" \
  -c 5 \
  -o distinct \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp25.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp25.txt -b "$CDS" \
| cut -f1-49,51 \
| awk 'BEGIN{FS=OFS="\t"} {$50=($50==-1||$50==".")?".":"CDS"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.txt

bedtools sort -i "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.1.txt

bedtools map \
  -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.1.txt \
  -b "$OREGANNO" -c 4 -o distinct \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp27.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp27.txt -b "$CENTROMERES" \
| cut -f1-51,53 \
| awk 'BEGIN{FS=OFS="\t"} {$52=($52==-1||$52==".")?".":"centromeric"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp28.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp28.txt -b "$PERICENTROMERES" \
| cut -f1-52,54 \
| awk 'BEGIN{FS=OFS="\t"} {$53=($53==-1||$53==".")?".":"pericentromeric"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp29.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp29.txt -b "$TELOMERES" \
| cut -f1-53,55 \
| awk 'BEGIN{FS=OFS="\t"} {$54=($54==-1||$54==".")?".":"telomeric"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp30.txt

###VAMOS — type-specific intersection logic
awk 'BEGIN{FS=OFS="\t"} $7!="INS" && $7!="DEL" && $7!="DUP" {print $0, "."}' \
"$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp30.txt > "$UNSOLVED_DIR"/novamos_temp.txt

bedtools slop -b 50 -i "$VAMOS" -g "$GENOME_FILE" \
| bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp30.txt -b - \
| awk -F'\t' '$7=="INS"' > "$UNSOLVED_DIR"/ins_temp.txt

bedtools slop -b 50 -i "$VAMOS" -g "$GENOME_FILE" \
| bedtools intersect -wa -wb -loj \
  -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp30.txt -b - \
| awk -F'\t' '$7=="DUP"' > "$UNSOLVED_DIR"/dup_temp.txt

bedtools intersect -wa -wb -loj -f 0.5 \
  -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp30.txt -b "$VAMOS" \
| awk -F'\t' '$7=="DEL"' > "$UNSOLVED_DIR"/del_temp.txt

cat "$UNSOLVED_DIR"/ins_temp.txt \
    "$UNSOLVED_DIR"/del_temp.txt \
    "$UNSOLVED_DIR"/dup_temp.txt \
| cut -f1-54,56 \
| awk 'BEGIN{FS=OFS="\t"} {$55=($55==-1||$55==".")?".":"vamos"; print}' \
> "$UNSOLVED_DIR"/vamos_annot_temp.txt

cat "$UNSOLVED_DIR"/novamos_temp.txt "$UNSOLVED_DIR"/vamos_annot_temp.txt \
| sort | uniq \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp31.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp31.txt -b "$SEGDUPS" \
| cut -f1-55,57 \
| awk 'BEGIN{FS=OFS="\t"} {$56=($56==-1||$56==".")?".":"segdups"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp32.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp32.txt -b "$REPEAT_MASKER" \
| cut -f1-56,58 \
| awk 'BEGIN{FS=OFS="\t"} {$57=($57==-1||$57==".")?".":"repeat"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp33.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp33.txt -b "$GAPS" \
| cut -f1-57,59 \
| awk 'BEGIN{FS=OFS="\t"} {$58=($58==-1||$58==".")?".":"Gap"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp34.txt

bedtools intersect -wa -wb -loj -f 1.0 \
  -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp34.txt -b "$HOMOPOLYMERS" \
| cut -f1-58,60 \
| awk 'BEGIN{FS=OFS="\t"} {$59=($59==-1||$59==".")?".":"HP>=50bp"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp35.txt

bedtools intersect -wa -wb -loj -f 1.0 \
  -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp35.txt -b "$DEFRABB_HICONF" \
| cut -f1-59,61 \
| awk 'BEGIN{FS=OFS="\t"} {$60=($60==-1||$60==".")?".":"hiconf"; print}' \
| sort -u \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp36.txt

##---------------------------------------------------------------------------
## STEP 14: Reorder to final output column layout
##
## temp36 column map:
##   $1-$7   : CHROM POS End_Pos REF ALT SVLEN SVTYPE
##   $8       : Cohort_support
##   $9       : 1KGP_support
##   $10-$13  : Cohort freq (Pop_Count, Pop_Freq, Allele_Count, Allele_Freq)
##   $14-$25  : 1KGP Pop_Count/Freq pairs (AFR AMR EAS EUR SAS ALL)
##   $26-$37  : 1KGP Allele_Count/Freq pairs (AFR AMR EAS EUR SAS ALL)
##   $38-$40  : GT_homWT, GT_het, GT_homVAR
##   $41-$42  : p, q
##   $43      : HWE
##   $44      : Genes
##   $45      : OMIM
##   $46      : GenCC
##   $47      : HPO
##   $48      : pLI
##   $49      : UTR
##   $50      : CDS
##   $51      : ORegAnno
##   $52      : Centromeric
##   $53      : Pericentromeric
##   $54      : Telomeric
##   $55      : Vamos
##   $56      : Segdup
##   $57      : Repeat
##   $58      : Gap
##   $59      : Homopolymer
##   $60      : HiConf
##
## Final output order:
##   $1-$7   : Chr Start End REF ALT SV_Length SV_Type
##   $8       : Cohort_support
##   $9       : 1KGP_support
##   $10-$13  : Cohort_Pop_Count Cohort_Pop_Freq Cohort_Allele_Count Cohort_Allele_Freq
##   $14      : Allele_Freq_ALL_1KGP           ← $37 from temp36
##   $15-$31  : Annotations (17 cols)
##   $32-$43  : 1KGP Pop_Count/Freq pairs (AFR AMR EAS EUR SAS ALL)
##   $44-$55  : 1KGP Allele_Count/Freq pairs (AFR AMR EAS EUR SAS ALL)
##   $56-$58  : GT_homWT, GT_het, GT_homVAR
##   $59      : HWE
##---------------------------------------------------------------------------

awk 'BEGIN{FS=OFS="\t"} {
    print $1,$2,$3,$4,$5,$6,$7,
          $8,$9,
          $10,$11,$12,$13,
          $37,
          $44,$45,$46,$47,$48,$49,$50,$51,$52,$53,$54,$55,$56,$57,$58,$59,$60,
          $14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,
          $26,$27,$28,$29,$30,$31,$32,$33,$34,$35,$36,$37,
          $38,$39,$40,
          $43
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp36.txt \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp1.txt

##---------------------------------------------------------------------------
## STEP 15: Truncate oversized fields, write header, assemble RESULTS files
##---------------------------------------------------------------------------

for col in 4 5 16; do
    awk -F'\t' -v c="$col" 'BEGIN{OFS="\t"} {
        if (length($c) > 30000) $c = "[too long, see vcf]"
        print
    }' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp1.txt \
    > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp1b.txt
    mv "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp1b.txt \
       "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp1.txt
done

echo -e "Chr\tStart_Pos\tEnd_Pos\tREF\tALT\tSV_Length\tSV_Type\tCohort_support\t1KGP_support\tCohort_Pop_Count\tCohort_Pop_Freq\tCohort_Allele_Count\tCohort_Allele_Freq\tAllele_Freq_ALL_1KGP\tGenes\tOMIM\tGenCC\tHPO\tpLI\tUTR\tCDS\tORegAnno\tCentromeric\tPericentromeric\tTelomeric\tVamos\tSegdup\tRepeat\tGap\tHomopolymer\tHiConf\tPop_Count_AFR\tPop_Freq_AFR\tPop_Count_AMR\tPop_Freq_AMR\tPop_Count_EAS\tPop_Freq_EAS\tPop_Count_EUR\tPop_Freq_EUR\tPop_Count_SAS\tPop_Freq_SAS\tPop_Count_ALL\tPop_Freq_ALL\tAllele_Count_AFR\tAllele_Freq_AFR\tAllele_Count_AMR\tAllele_Freq_AMR\tAllele_Count_EAS\tAllele_Freq_EAS\tAllele_Count_EUR\tAllele_Freq_EUR\tAllele_Count_SAS\tAllele_Freq_SAS\tAllele_Count_ALL\tAllele_Freq_ALL\tGT_homWT\tGT_het\tGT_homVAR\tHWE" \
> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS.txt

sed -i 's/?//g' "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp1.txt"
sed -i -e :a -e '/^[[:space:]]*$/{$d;N;ba}' "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp1.txt"
cat "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp1.txt" >> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.txt"

##Filtered output 1: SVs absent from 1KGP entirely ($14 = Allele_Freq_ALL_1KGP)
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 || $14==0 {print}' \
"${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.txt" \
> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_cohort_unique.txt"

##Filtered output 2: SVs with 1KGP allele freq <=0.01
awk -F'\t' 'BEGIN{OFS="\t"} NR==1 || $14<=0.01 {print}' \
"${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.txt" \
> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_0.01.txt"

##---------------------------------------------------------------------------
## STEP 16: Build output VCF
##---------------------------------------------------------------------------

cat <<EOF > "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp.vcf"
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chr10,length=133797422>
##contig=<ID=chr11,length=135086622>
##contig=<ID=chr12,length=133275309>
##contig=<ID=chr13,length=114364328>
##contig=<ID=chr14,length=107043718>
##contig=<ID=chr15,length=101991189>
##contig=<ID=chr16,length=90338345>
##contig=<ID=chr17,length=83257441>
##contig=<ID=chr18,length=80373285>
##contig=<ID=chr19,length=58617616>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##contig=<ID=chrM,length=16569>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=Cohort_support,Number=.,Type=String,Description="Cohort samples carrying the SV">
##INFO=<ID=1KGP_support,Number=.,Type=String,Description="1KGP control samples (n=500) sharing the SV">
##INFO=<ID=Cohort_Pop_Count,Number=1,Type=Integer,Description="Number of cohort individuals carrying the SV">
##INFO=<ID=Cohort_Pop_Freq,Number=1,Type=Float,Description="Proportion of cohort individuals carrying the SV">
##INFO=<ID=Cohort_Allele_Count,Number=1,Type=Integer,Description="Number of cohort alleles carrying the SV">
##INFO=<ID=Cohort_Allele_Freq,Number=1,Type=Float,Description="Allele frequency of the SV in the cohort">
##INFO=<ID=Allele_Freq_ALL_1KGP,Number=1,Type=Float,Description="Allele frequency of the SV in the 1KGP control cohort (n=500)">
##INFO=<ID=Genes,Number=.,Type=String,Description="Genes overlapped by the SV (gencode v45, +/-5kb)">
##INFO=<ID=OMIM,Number=.,Type=String,Description="OMIM phenotypes for overlapping genes (OMIM 20251213)">
##INFO=<ID=GenCC,Number=.,Type=String,Description="GenCC phenotypes for overlapping genes (GenCC 20250510)">
##INFO=<ID=HPO,Number=.,Type=String,Description="HPO terms for overlapping genes">
##INFO=<ID=pLI,Number=.,Type=String,Description="pLI scores from gnomAD v4.1">
##INFO=<ID=UTR,Number=1,Type=String,Description="SV overlaps a UTR (gencode v45)">
##INFO=<ID=CDS,Number=1,Type=String,Description="SV overlaps a coding exon (gencode v45)">
##INFO=<ID=ORegAnno,Number=1,Type=String,Description="SV overlaps a regulatory element (ORegAnno)">
##INFO=<ID=Centromeric,Number=1,Type=String,Description="SV overlaps a centromeric region (UCSC hg38)">
##INFO=<ID=Pericentromeric,Number=1,Type=String,Description="SV overlaps a pericentromeric region (+/-5Mb)">
##INFO=<ID=Telomeric,Number=1,Type=String,Description="SV overlaps a telomeric region (5Mb from chromosome end)">
##INFO=<ID=Vamos,Number=1,Type=String,Description="SV overlaps a tandem repeat region (VAMOS, n=148)">
##INFO=<ID=Segdup,Number=1,Type=String,Description="SV overlaps a segmental duplication (GIAB v3.3)">
##INFO=<ID=Repeat,Number=1,Type=String,Description="SV overlaps a repeat region (UCSC RepeatMasker)">
##INFO=<ID=Gap,Number=1,Type=String,Description="SV overlaps an hg38 assembly gap (UCSC)">
##INFO=<ID=Homopolymer,Number=1,Type=String,Description="SV overlaps a homopolymeric region >=50bp">
##INFO=<ID=HiConf,Number=1,Type=String,Description="SV fully contained in a high-confidence region (GIAB T2TQ100-V1.0)">
##INFO=<ID=Pop_Count_AFR,Number=1,Type=Integer,Description="1KGP AFR samples carrying the SV (n=163)">
##INFO=<ID=Pop_Freq_AFR,Number=1,Type=Float,Description="Proportion of 1KGP AFR samples carrying the SV">
##INFO=<ID=Pop_Count_AMR,Number=1,Type=Integer,Description="1KGP AMR samples carrying the SV (n=67)">
##INFO=<ID=Pop_Freq_AMR,Number=1,Type=Float,Description="Proportion of 1KGP AMR samples carrying the SV">
##INFO=<ID=Pop_Count_EAS,Number=1,Type=Integer,Description="1KGP EAS samples carrying the SV (n=95)">
##INFO=<ID=Pop_Freq_EAS,Number=1,Type=Float,Description="Proportion of 1KGP EAS samples carrying the SV">
##INFO=<ID=Pop_Count_EUR,Number=1,Type=Integer,Description="1KGP EUR samples carrying the SV (n=72)">
##INFO=<ID=Pop_Freq_EUR,Number=1,Type=Float,Description="Proportion of 1KGP EUR samples carrying the SV">
##INFO=<ID=Pop_Count_SAS,Number=1,Type=Integer,Description="1KGP SAS samples carrying the SV (n=103)">
##INFO=<ID=Pop_Freq_SAS,Number=1,Type=Float,Description="Proportion of 1KGP SAS samples carrying the SV">
##INFO=<ID=Pop_Count_ALL,Number=1,Type=Integer,Description="1KGP samples carrying the SV (all populations, n=500)">
##INFO=<ID=Pop_Freq_ALL,Number=1,Type=Float,Description="Proportion of 1KGP samples carrying the SV (all populations)">
##INFO=<ID=Allele_Count_AFR,Number=1,Type=Integer,Description="1KGP AFR alleles carrying the SV">
##INFO=<ID=Allele_Freq_AFR,Number=1,Type=Float,Description="1KGP AFR allele frequency">
##INFO=<ID=Allele_Count_AMR,Number=1,Type=Integer,Description="1KGP AMR alleles carrying the SV">
##INFO=<ID=Allele_Freq_AMR,Number=1,Type=Float,Description="1KGP AMR allele frequency">
##INFO=<ID=Allele_Count_EAS,Number=1,Type=Integer,Description="1KGP EAS alleles carrying the SV">
##INFO=<ID=Allele_Freq_EAS,Number=1,Type=Float,Description="1KGP EAS allele frequency">
##INFO=<ID=Allele_Count_EUR,Number=1,Type=Integer,Description="1KGP EUR alleles carrying the SV">
##INFO=<ID=Allele_Freq_EUR,Number=1,Type=Float,Description="1KGP EUR allele frequency">
##INFO=<ID=Allele_Count_SAS,Number=1,Type=Integer,Description="1KGP SAS alleles carrying the SV">
##INFO=<ID=Allele_Freq_SAS,Number=1,Type=Float,Description="1KGP SAS allele frequency">
##INFO=<ID=Allele_Count_ALL,Number=1,Type=Integer,Description="1KGP alleles carrying the SV (all populations)">
##INFO=<ID=Allele_Freq_ALL,Number=1,Type=Float,Description="1KGP allele frequency (all populations)">
##INFO=<ID=GT_homWT,Number=1,Type=Integer,Description="Count of hom-ref genotypes in 1KGP">
##INFO=<ID=GT_het,Number=1,Type=Integer,Description="Count of het genotypes in 1KGP">
##INFO=<ID=GT_homVAR,Number=1,Type=Integer,Description="Count of hom-alt genotypes in 1KGP">
##INFO=<ID=HWE,Number=1,Type=String,Description="Hardy-Weinberg Equilibrium result (1KGP)">
EOF

echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' \
>> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp.vcf"

awk -F"\t" 'BEGIN{OFS="\t"} NR>1 {
    info = "END=" $3 \
         ";SVLEN=" $6 ";SVTYPE=" $7 \
         ";Cohort_support=" $8 ";1KGP_support=" $9 \
         ";Cohort_Pop_Count=" $10 ";Cohort_Pop_Freq=" $11 \
         ";Cohort_Allele_Count=" $12 ";Cohort_Allele_Freq=" $13 \
         ";Allele_Freq_ALL_1KGP=" $14 \
         ";Genes=" $15 ";OMIM=" $16 ";GenCC=" $17 ";HPO=" $18 ";pLI=" $19 \
         ";UTR=" $20 ";CDS=" $21 ";ORegAnno=" $22 \
         ";Centromeric=" $23 ";Pericentromeric=" $24 ";Telomeric=" $25 \
         ";Vamos=" $26 ";Segdup=" $27 ";Repeat=" $28 \
         ";Gap=" $29 ";Homopolymer=" $30 ";HiConf=" $31 \
         ";Pop_Count_AFR=" $32 ";Pop_Freq_AFR=" $33 \
         ";Pop_Count_AMR=" $34 ";Pop_Freq_AMR=" $35 \
         ";Pop_Count_EAS=" $36 ";Pop_Freq_EAS=" $37 \
         ";Pop_Count_EUR=" $38 ";Pop_Freq_EUR=" $39 \
         ";Pop_Count_SAS=" $40 ";Pop_Freq_SAS=" $41 \
         ";Pop_Count_ALL=" $42 ";Pop_Freq_ALL=" $43 \
         ";Allele_Count_AFR=" $44 ";Allele_Freq_AFR=" $45 \
         ";Allele_Count_AMR=" $46 ";Allele_Freq_AMR=" $47 \
         ";Allele_Count_EAS=" $48 ";Allele_Freq_EAS=" $49 \
         ";Allele_Count_EUR=" $50 ";Allele_Freq_EUR=" $51 \
         ";Allele_Count_SAS=" $52 ";Allele_Freq_SAS=" $53 \
         ";Allele_Count_ALL=" $54 ";Allele_Freq_ALL=" $55 \
         ";GT_homWT=" $56 ";GT_het=" $57 ";GT_homVAR=" $58 \
         ";HWE=" $59
    print $1, $2, ".", $4, $5, ".", ".", info
}' "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.txt" \
>> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp.vcf"

bcftools sort "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp.vcf" \
> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.vcf"

bgzip "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.vcf"
tabix "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.vcf.gz"

##This removes intermediate temp files
rm ${UNSOLVED_DIR}/*temp*
