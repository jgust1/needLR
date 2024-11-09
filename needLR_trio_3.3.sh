#!/bin/bash

##Run from inside needLR_local

while getopts p:m:f:g: flag
do
    case "${flag}" in
        p) PROBAND_VCF=${OPTARG};;
        m) MATERNAL_VCF=${OPTARG};;
        f) PATERNAL_VCF=${OPTARG};;
        g) REF_GENOME=${OPTARG};;
    esac
done

PROBAND_PREFIX=$(basename "$PROBAND_VCF" .vcf.gz)
MATERNAL_PREFIX=$(basename "$MATERNAL_VCF" .vcf.gz)
PATERNAL_PREFIX=$(basename "$PATERNAL_VCF" .vcf.gz)

QUERY_FILE_NAME=$(basename "$PROBAND_VCF" .vcf.gz)_needLR-3.3_trio_UWONT150

##This is a filepath to an hg38 reference FASTA - USER ADDS
REF_GENOME="/n/users/jgust1/reference_files/hg38_reference_genome/hg38.no_alt.fa"

UNSOLVED_DIR="needLR_output/$QUERY_FILE_NAME"
mkdir $UNSOLVED_DIR

##This is the pre-bcftools merged background vcf
ONTUW150_BCFTOOLS_MERGED="backend_files/UWONT150_truvari.vcf.gz"

##This is the list of the 1KGP sample names in the order they will be merged in Jasmine
KGP_SAMPLE_NAMES="backend_files/ID_order.txt"

##These are the bed files used for annotating the SVs
GENES=""$WD_PREFIX"backend_files/bed_files/PROTEIN_CODING_GENE_gencode.v45.annotation_1kb_slop.bed"
CODING_REGIONS=""$WD_PREFIX"backend_files/bed_files/ENSEMBL_CANONICAL_EXON_in_PROTEIN_CODING_GENE_gencode.v45.annotation.bed"
OMIM_GENES=""$WD_PREFIX"backend_files/bed_files/OMIM_gene_phen_hg38.bed"
CENTROMERES=""$WD_PREFIX"backend_files/bed_files/hg38_centromeres_endtoend.bed"
PERICENTROMERES=""$WD_PREFIX"backend_files/bed_files/hg38_pericentromeres_5Mb.bed"
TELOMERES=""$WD_PREFIX"backend_files/bed_files/hg38_telomeres_5Mb.bed"
DEFRABB_HICONF=""$WD_PREFIX"backend_files/bed_files/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed"
STR=""$WD_PREFIX"backend_files/bed_files/STR_original_motifs.set148.bed"
VNTR=""$WD_PREFIX"backend_files/bed_files/VNTR_original_motifs.set148.bed"
SEGDUPS=""$WD_PREFIX"backend_files/bed_files/SEGDUPS_GIAB_v3.3.bed"
REPEAT_MASKER=""$WD_PREFIX"backend_files/bed_files/UCSC_hg38_Repeats_RepeatMasker.bed"
GAPS=""$WD_PREFIX"backend_files/bed_files/UCSC_hg38_Mapping_and_Sequencing_Gap.bed"

##Copy trio vcfs and indexes into putput directory
scp "$PROBAND_VCF"* $UNSOLVED_DIR/
scp "$MATERNAL_VCF"* $UNSOLVED_DIR/
scp "$PATERNAL_VCF"* $UNSOLVED_DIR/

#This preprocesses the query VCFs to be >=50bp, on full chromosomes, and FILTER=PASS
for vcf in $UNSOLVED_DIR/*vcf.gz; do
preprocessed_vcf=$(basename "$vcf" .gz)
bcftools view -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o "$UNSOLVED_DIR"/preprocessed_"$preprocessed_vcf".vcf "$vcf"
done

bgzip "$UNSOLVED_DIR"/preprocessed*

find "$UNSOLVED_DIR" -type f -name "preprocessed*gz" -exec tabix {} \;

#This creates the input list of VCFs for bcftools merge
echo "$UNSOLVED_DIR"/preprocessed_"$PROBAND_PREFIX"*.vcf.gz > "$UNSOLVED_DIR"/list.txt
echo "$UNSOLVED_DIR"/preprocessed_"$MATERNAL_PREFIX"*.vcf.gz >> "$UNSOLVED_DIR"/list.txt
echo "$UNSOLVED_DIR"/preprocessed_"$PATERNAL_PREFIX"*.vcf.gz >> "$UNSOLVED_DIR"/list.txt
echo "$ONTUW150_BCFTOOLS_MERGED" >> "$UNSOLVED_DIR"/list.txt

#bcftools merge
bcftools merge -m none --force-samples -l "$UNSOLVED_DIR"/list.txt -Oz -o "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_with_UWONT150.vcf.gz

tabix "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_with_UWONT150.vcf.gz

#Truvari collapse
truvari collapse -i "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_with_UWONT150.vcf.gz -o "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged.vcf -c "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_collapsed.vcf -f $REF_GENOME -k first -s 50 -S 10000000 --passonly -B 50 -r 2000 -p 0 -P 0.2 -O 0.2 --chain

##This parses the Truvari VCF output for relevant information
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVLEN\t%INFO/SVTYPE[\t%GT\t%DV\t%DR]\n' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged.vcf > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged.txt

#Only keep rows where column 6 is not "./."
awk 'BEGIN {OFS="\t"} $6 != "./."' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp.txt

##Rearrange the columns so that you keep the GT, DV, DR, (and total read depth) for any SVs in the query, and then order the GT outputs for the 1KGP samples
awk 'BEGIN {OFS="\t"} {print $1, $2, $4, $5, $6, $7, $8, ($7 + $8), $9, $12, $15, $18, $21, $24, $27, $30, $33, $36, $39, $42, $45, $48, $51, $54, $57, $60, $63, $66, $69, $72, $75, $78, $81, $84, $87, $90, $93, $96, $99, $102, $105, $108, $111, $114, $117, $120, $123, $126, $129, $132, $135, $138, $141, $144, $147, $150, $153, $156, $159, $162, $165, $168, $171, $174, $177, $180, $183, $186, $189, $192, $195, $198, $201, $204, $207, $210, $213, $216, $219, $222, $225, $228, $231, $234, $237, $240, $243, $246, $249, $252, $255, $258, $261, $264, $267, $270, $273, $276, $279, $282, $285, $288, $291, $294, $297, $300, $303, $306, $309, $312, $315, $318, $321, $324, $327, $330, $333, $336, $339, $342, $345, $348, $351, $354, $357, $360, $363, $366, $369, $372, $375, $378, $381, $384, $387, $390, $393, $396, $399, $402, $405, $408, $411, $414, $417, $420, $423, $426, $429, $432, $435, $438, $441, $444, $447, $450, $453, $456, $459, $462}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp_2.txt

##This creates a SUPP_VEC (like in Jasmine)
awk 'BEGIN {OFS="\t"} {
    result = "";
    for (i=11; i<=160; i++) {
        if ($i ~ /1/) {
            result = result "1";
        } else {
            result = result "0";
        }
    }
    print $0, result
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp_2.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp_3.txt

##Rearrange columns so SUPP_VEC is $6
awk 'BEGIN {OFS="\t"} {print $1, $2, "temp", $3, $4, $161, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54, $55, $56, $57, $58, $59, $60, $61, $62, $63, $64, $65, $66, $67, $68, $69, $70, $71, $72, $73, $74, $75, $76, $77, $78, $79, $80, $81, $82, $83, $84, $85, $86, $87, $88, $89, $90, $91, $92, $93, $94, $95, $96, $97, $98, $99, $100, $101, $102, $103, $104, $105, $106, $107, $108, $109, $110, $111, $112, $113, $114, $115, $116, $117, $118, $119, $120, $121, $122, $123, $124, $125, $126, $127, $128, $129, $130, $131, $132, $133, $134, $135, $136, $137, $138, $139, $140, $141, $142, $143, $144, $145, $146, $147, $148, $149, $150, $151, $152, $153, $154, $155, $156, $157, $158, $159, $160}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp_3.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp_4.txt

##This removes BNDs from the analysis (too large to efficiently annotate)
cat "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp_4.txt | grep -v BND > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp1.txt

## AFR – 51
## AMR – 18
## EAS – 25
## EUR – 24
## SAS – 32

##This calculates the population frequencies, allele frequencies, and allele types (for HWE)
awk 'BEGIN {OFS="\t"} {
   substring_AFR = substr($6, 1, 51)
   ## 2 = starting position, 51 = length
   substring_AMR = substr($6, 52, 18)
   substring_EAS = substr($6, 70, 25)
   substring_EUR = substr($6, 95, 24)
   substring_SAS = substr($6, 119, 32)
   substring_ALL = substr($6, 1, 150)

    ## Count the number of 1s in the substring
   Pop_AFR = gsub("1", "", substring_AFR)
   Pop_AMR = gsub("1", "", substring_AMR)
   Pop_EAS = gsub("1", "", substring_EAS)
   Pop_EUR = gsub("1", "", substring_EUR)
   Pop_SAS = gsub("1", "", substring_SAS)
   Pop_ALL = gsub("1", "", substring_ALL)

   Allele_AFR = 0;
   Allele_AMR = 0;
   Allele_EAS = 0;
   Allele_EUR = 0;
   Allele_SAS = 0;
   Allele_ALL = 0;
   GT_homWT = 0;
   GT_het = 0;
   GT_homVAR = 0;
   for (i = 11; i <= 61; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_AFR++;
       }
     }
   }
   for (i = 62; i <= 79; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_AMR++;
       }
     }
   }
   for (i = 80; i <= 104; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_EAS++;
       }
     }
   }
   for (i = 105; i <= 128; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_EUR++;
       }
     }
   }
   for (i = 129; i <= 160; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_SAS++;
       }
     }
   }
   for (i = 11; i <= 160; i++) {
      n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_ALL++;
       }
     }
   }
   for (i = 11; i <= 160; i++) {
     if ($i == "./.") {
           GT_homWT++;
       }
     }
   for (i = 11; i <= 160; i++) {
     if ($i == "0/1" || $i == "1/0") {
           GT_het++;
       }
     }
   for (i = 11; i <= 160; i++) {
     if ($i == "1/1") {
          GT_homVAR++;
       }
     }

   print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, Pop_AFR, Pop_AMR, Pop_EAS, Pop_EUR, Pop_SAS, Pop_ALL, Allele_AFR, Allele_AMR, Allele_EAS, Allele_EUR, Allele_SAS, Allele_ALL, GT_homWT, GT_het, GT_homVAR
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp1.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp2.txt


##This matches the SUPP_VEC string to the sample list (and lists which 1KGP samples have the SV)
mapfile -t names < "$KGP_SAMPLE_NAMES"

awk -v names="$(printf "%s\n" "${names[@]}")" '
BEGIN {
    split(names, nameArray, "\n");
    OFS="\t";
}
{
    # Get the binary string from column 6
    binary_string = $6;
    
    # Check if the binary string is all zeros
    if (binary_string ~ /^0+$/) {
        names_string = "query_only";
    } else {
        # Initialize the names string
        names_string = "";
        
        # Loop over each character in the binary string
        for (i = 1; i <= length(binary_string); i++) {
            if (substr(binary_string, i, 1) == "1") {
                if (names_string == "") {
                    names_string = nameArray[i];
                } else {
                    names_string = names_string "_" nameArray[i];
                }
            }
        }
    }
    
    # Replace column 6 with the names string
    $6 = names_string;
    
    # Print the modified row
    print $0;
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp2.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp3.txt


# $27 - GT_homWT 
# $28 - GT_het 
# $29 - GT_homVAR

## This calculates p for HWE (based on 150 samples) - WT alle freq ((2*homWT) + het)/(2* (homWT + het + homVar))
## Checks for demoninator of zero (if SV is counted but genotype is ./. for a 1KGP sample)

awk '{
  if ($25 + $26 + $27 != 0) {
    result = ((2 * $25) + $26) / (2 * ($25 + $26 + $27));
  } else {
    result = 0;
  }
  print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, result
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp3.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp4.txt

##This calculates q for HWE (based on 150 samples) - (1-p)
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, (1-$28)}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp4.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp5.txt

##This calculates np^2, 2npq, and nq^2 for HWE (based on 150 samples)
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, (150*$28*$28), (300*$28*$29), (150*$29*$29)}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp5.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp6.txt

#This calculates chi-squared value for HWE. (Can't divide by zero, so if expected is zero, chi-sq=NULL) 
awk '{
   if ($30 > 0 && $31 > 0 && $32 > 0) {
       result = ((($25 - $30) * ($25 - $30)) / $30) + ((($26 - $31) * ($26 - $31)) / $31) + ((($27 - $32) * ($27 - $32)) / $32);
       print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, result;
   } else {
       print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, "NULL";
   }
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp6.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp7.txt

## $15 - Pop_AFR 
## $16 - Pop_AMR 
## $17 - Pop_EAS
## $18 - Pop_EUR
## $19 - Pop_SAS
## $20 - Pop_ALL
## $21 - Allele_AFR
## $22 - Allele_AMR
## $23 - Allele_EAS
## $24 - Allele_EUR
## $25 - Allele_SAS
## $26 - Allele_ALL
## $27 - GT_homWT
## $28 - T_het
## $29 - GT_homVAR


#This lists all population counts/frequencies, allele counts/frequencies, and HWE info
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, ($13/51), $14, ($14/18), $15, ($15/25), $16, ($16/24), $17, ($17/32), $18, ($18/150), $19, ($19/102), $20, ($20/36), $21, ($21/50), $22, ($22/48), $23, ($23/64), $24, ($24/300), $25, $26, $27, $28, $29, $33}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp7.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp8.txt

#This determines if the SV is in HWE
awk '{
   if ($42 == "NULL") {
       $42 = "HWE_NULL";
   } else if ($42 > 3.84) {
       $42 = "HWE_FALSE";
   } else {
       $42 = "HWE_TRUE";
   }
   print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp8.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp9.txt

##This adds the absolute value of the length of the SV to the start position (DELs are listed as negative lengths in the vcf, but the range is from the start position plus that number of nucleotides)
awk 'BEGIN {OFS="\t"} { if ($5 == "INS") $3 = $2 + 1; else $3 = $2 + ($4 < 0 ? -1*$4 : $4); print }' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp9.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp10.txt

#This annotates each variant with gene info (+/-1kb to try to capture promoters, etc)
bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp10.txt -b $GENES > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp11.txt

#This sorts the file by gene name
sort -k46,46 "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp11.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp12.txt

#This joins the gene name with any OMIM annotations
join -t $'\t' -a 1 -1 46 -2 4 -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.46,2.6' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp12.txt "$OMIM_GENES" > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp13.txt

#This condolidates rows that are identical except for the gene and lists the genes together in the last column (ie. a large deletion that spans 50 genes will be consolidated to one row)
awk 'BEGIN{FS=OFS="\t"} {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8 FS $9 FS $10 FS $11 FS $12 FS $13 FS $14 FS $15 FS $16 FS $17 FS $18 FS $19 FS $20 FS $21 FS $22 FS $23 FS $24 FS $25 FS $26 FS $27 FS $28 FS $29 FS $30 FS $31 FS $32 FS $33 FS $34 FS $35 FS $36 FS $37 FS $38 FS $39 FS $40 FS $41 FS $42} {a[key]=(a[key]!="")?(a[key] ", " $43):(a[key] $43)} {b[key]=(b[key]!="")?(b[key] ", " $44):(b[key] $44)} END {for (k in a) print k, a[k], b[k]}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp13.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp14.txt

#This replaces all empoty fields with a . (important for next command)
awk 'BEGIN {FS = OFS = "\t"} {for (i = 1; i <= NF; i++) if ($i == "") $i = "."} 1' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp14.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp15.txt

#These commands annotates each SV with info (for each bed file), removes unnecessary columns, labels annotated regions as such, and sorts the file (necessary for bedtools intersect)

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp15.txt -b $CODING_REGIONS | \
cut -f1-44,46 | \
awk 'BEGIN {FS=OFS="\t"} $45 == -1 {$45 = "."} $45 != -1 && $45 != "." {$45 = "exonic"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp16.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp16.txt -b $CENTROMERES | \
cut -f1-45,47 | \
awk 'BEGIN {FS=OFS="\t"} $46 == -1 {$46 = "."} $46 != -1 && $46 != "." {$46 = "centromeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp17.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp17.txt -b $PERICENTROMERES | \
cut -f1-46,48 | \
awk 'BEGIN {FS=OFS="\t"} $47 == -1 {$47 = "."} $47 != -1 && $47 != "." {$47 = "pericentromeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp18.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp18.txt -b $TELOMERES | \
cut -f1-47,49 | \
awk 'BEGIN {FS=OFS="\t"} $48 == -1 {$48 = "."} $48 != -1 && $48 != "." {$48 = "telomeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp19.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp19.txt -b $STR | \
cut -f1-48,50 | \
awk 'BEGIN {FS=OFS="\t"} $49 == -1 {$49 = "."} $49 != -1 && $49 != "." {$49 = "STR"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp20.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp20.txt -b $VNTR | \
cut -f1-49,51 | \
awk 'BEGIN {FS=OFS="\t"} $50 == -1 {$50 = "."} $50 != -1 && $50 != "." {$50 = "VNTR"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp21.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp21.txt -b $SEGDUPS | \
cut -f1-50,52 | \
awk 'BEGIN {FS=OFS="\t"} $51 == -1 {$51 = "."} $51 != -1 && $51 != "." {$51 = "segdups"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp22.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp22.txt -b $REPEAT_MASKER | \
cut -f1-51,53 | \
awk 'BEGIN {FS=OFS="\t"} $52 == -1 {$52 = "."} $52 != -1 && $52 != "." {$52 = "repeat"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp23.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp23.txt -b $GAPS | \
cut -f1-52,54 | \
awk 'BEGIN {FS=OFS="\t"} $53 == -1 {$53 = "."} $53 != -1 && $53 != "." {$53 = "gap"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.txt -b $DEFRABB_HICONF  -f 1.0 | \
cut -f1-53,55 | \
awk 'BEGIN {FS=OFS="\t"} $54 == -1 {$54 = "."} $54 != -1 && $54 != "." {$54 = "hiconf"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp25.txt

## $1-$12 - SV info (Chr, Start_Pos, End_Pos, SV_Length, SV_Type, 1KGP_support, Genotype, Alt_reads, Ref_reads, Total_reads, mat, pat)
## $13-$36 - allele/pop freq ($34 = AlleleFreqAll)
## $37-$42 - HWE (GT_homWT, GThet, GT_homVar, p, q, HWE)
## $43-$54 - Annotations (gene, OMIM, exonic, cent, peri, tel, str, vntr, segdup, repeat, gap, hiconf)

#This prints the columns in the correct order
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $36, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp25.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.txt

awk 'BEGIN {FS=OFS="\t"} {
    if ($11 == "./." && $12 == "./.") {
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "de novo", $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54, $55
    } else {
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, "inherited", $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54, $55
    }
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp27.txt


#This creates an output file with the correct column headers
echo -e "Chr\tStart_Pos\tEnd_Pos\tSV_Length\tSV_Type\t1KGP_support\tGenotype\tAlt_reads\tRef_reads\tTotal_reads\tMaternal_genotype\tPaternal_genotype\tInheritance\tAllele_Freq_ALL\tGenes\tOMIM\tExonic\tCentromeric\tPericentromeric\tTelomeric\tSTR\tVNTR\tSegdup\tRepeat\tGap\tHiConf\tPop_Count_AFR\tPop_Freq_AFR\tPop_Count_AMR\tPop_Freq_AMR\tPop_Count_EAS\tPop_Freq_EAS\tPop_Count_EUR\tPop_Freq_EUR\tPop_Count_SAS\tPop_Freq_SAS\tPop_Count_ALL\tPop_Freq_ALL\tAllele_Count_AFR\tAllele_Freq_AFR\tAllele_Count_AMR\tAllele_Freq_AMR\tAllele_Count_EAS\tAllele_Freq_EAS\tAllele_Count_EUR\tAllele_Freq_EUR\tAllele_Count_SAS\tAllele_Freq_SAS\tAllele_Count_ALL\tAllele_Freq_ALL\tGT_homWT\tGT_het\tGT_homVAR\tHWE-p\tHWE-q\tHWE" > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS.txt

#This adds the data to the output file
cat "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp27.txt >> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS.txt

#This removes "?s" from OMIM annotations (hard to sort) - optional
sed -i 's/?//g' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS.txt

##This removes intermediate "temp" files 
find "$UNSOLVED_DIR" -type f -name "*temp*" -exec rm {} \;
find "$UNSOLVED_DIR" -type f -name "*collapsed*" -exec rm {} \;
find "$UNSOLVED_DIR" -type f -name "*bcftools_merged*" -exec rm {} \;
find "$UNSOLVED_DIR" -type f -name "*merged.txt*" -exec rm {} \;

