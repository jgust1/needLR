#!/bin/bash

##Run from inside needLR_local


while getopts f:g:t: flag
do
    case "${flag}" in
        f) QUERY_VCFS=${OPTARG};;
        g) REF_GENOME=${OPTARG};;
        t) THREADS=${OPTARG};;
    esac
done

##This starts the loop for all of the files in QUERY_VCFS
while IFS= read -r QUERY_FILE_PATH; do

##This will be the basename for the output files
QUERY_FILE_NAME=$(basename "$QUERY_FILE_PATH" .vcf.gz)_needLR_3.2

# ##This is the output directory specific to each query VCF
UNSOLVED_DIR="needLR_output/$QUERY_FILE_NAME"
mkdir $UNSOLVED_DIR

#These are the file paths of the 1KGP input VCFs by ancestry (order important for allele frequencies)
AFR_LIST="backend_files/1KGP_VCFs/AFR_list.txt"
AMR_LIST="backend_files/1KGP_VCFs/AMR_list.txt"
EAS_LIST="backend_files/1KGP_VCFs/EAS_list.txt"
EUR_LIST="backend_files/1KGP_VCFs/EUR_list.txt"
SAS_LIST="backend_files/1KGP_VCFs/SAS_list.txt"

##This is the list of the 1KGP sample names in the order they will be merged in Jasmine
KGP_SAMPLE_NAMES="backend_files/1KGP_VCFs/ID_order.txt"

##These are the bed files used for annotating the SVs
GENES="backend_files/bed_files/PROTEIN_CODING_GENE_gencode.v45.annotation_1kb_slop.bed"
CODING_REGIONS="backend_files/bed_files/ENSEMBL_CANONICAL_EXON_in_PROTEIN_CODING_GENE_gencode.v45.annotation.bed"
OMIM_GENES="backend_files/bed_files/OMIM_gene_phen_hg38.bed"
CENTROMERES="backend_files/bed_files/hg38_centromeres_endtoend.bed"
PERICENTROMERES="backend_files/bed_files/hg38_pericentromeres_5Mb.bed"
TELOMERES="backend_files/bed_files/hg38_telomeres_5Mb.bed"
DEFRABB_HICONF="backend_files/bed_files/GRCh38_HG002-T2TQ100-V1.0_stvar.benchmark.bed"
STR="backend_files/bed_files/STR_original_motifs.set148.bed"
VNTR="backend_files/bed_files/VNTR_original_motifs.set148.bed"
SEGDUPS="backend_files/bed_files/SEGDUPS_GIAB_v3.3.bed"
REPEAT_MASKER="backend_files/bed_files/UCSC_hg38_Repeats_RepeatMasker.bed"
GAPS="backend_files/bed_files/UCSC_hg38_Mapping_and_Sequencing_Gap.bed"

#This preprocesses the query VCFs to be >=50bp, on full chromosomes, and FILTER=PASS
bcftools view -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o "$UNSOLVED_DIR"/preprocessed_"$QUERY_FILE_NAME".vcf "$QUERY_FILE_PATH"

#This creates the input list of VCFs for Jasmine 
echo "$UNSOLVED_DIR"/preprocessed_"$QUERY_FILE_NAME".vcf > "$UNSOLVED_DIR"/list.txt
cat "$AFR_LIST" >> "$UNSOLVED_DIR"/list.txt
cat "$AMR_LIST" >> "$UNSOLVED_DIR"/list.txt
cat "$EAS_LIST" >> "$UNSOLVED_DIR"/list.txt
cat "$EUR_LIST" >> "$UNSOLVED_DIR"/list.txt
cat "$SAS_LIST" >> "$UNSOLVED_DIR"/list.txt

#This merges the query VCF and all of the 1KGP VCFs 

#new2
jasmine file_list="$UNSOLVED_DIR"/list.txt out_file="$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_jasmine.vcf threads="$THREADS" --maxdist=500 --nonlinear_dist --centroid_merging --allow_intrasample --output_genotypes spec_reads=3 genome_file="$REF_GENOME"

#This parses the Jasmine output for all SVs seen in the query sample
bcftools view -i 'INFO/SUPP_VEC ~ "^1"' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_jasmine.vcf > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_jasmine_query_only.vcf

##This parses the Jasmine VCF output for relevant information - customizable, but could require updating the rest of the script
bcftools query -f '%CHROM\t%POS\t%ID\t%INFO/SVLEN\t%INFO/SVTYPE\t%INFO/SUPP_VEC[\t%GT\t%DV\t%DR]\n' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_jasmine_query_only.vcf > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_jasmine_query_only.txt

##Rearrange the columns so that you keep the GT, DV, DR, (and total read depth) for any SVs in the query, and then order the GT outputs for the 1KGP samples
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, ($8 + $9), $10, $13, $16, $19, $22, $25, $28, $31, $34, $37, $40, $43, $46, $49, $52, $55, $58, $61, $64, $67, $70, $73, $76, $79, $82, $85, $88, $91, $94, $97, $100, $103, $106, $109, $112, $115, $118, $121, $124, $127, $130, $133, $136, $139, $142, $145, $148, $151, $154, $157, $160, $163, $166, $169, $172, $175, $178, $181, $184, $187, $190, $193, $196, $199, $202, $205, $208, $211, $214, $217, $220, $223, $226, $229, $232, $235, $238, $241, $244, $247, $250, $253, $256, $259, $262, $265, $268, $271, $274, $277, $280, $283, $286, $289, $292, $295, $298, $301, $304, $307, $310, $313, $316, $319, $322, $325, $328, $331, $334, $337, $340, $343, $346, $349, $352, $355, $358, $361, $364, $367, $370, $373, $376, $379, $382, $385, $388, $391, $394, $397, $400, $403, $406, $409, $412, $415, $418, $421, $424, $427, $430, $433, $436, $439, $442, $445, $448, $451, $454, $457}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_jasmine_query_only.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_jasmine_query_only_2.txt

##This removes BNDs from the analysis (too large to efficiently annotate)
cat "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_jasmine_query_only_2.txt | grep -v BND > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp1.txt

## AFR – 51
## AMR – 18
## EAS – 25
## EUR – 24
## SAS – 32

##This calculates the population frequencies, allele frequencies, and allele types (for HWE)
awk '{
   substring_AFR = substr($6, 2, 51)
   ## 2 = starting position, 51 = length
   substring_AMR = substr($6, 53, 18)
   substring_EAS = substr($6, 71, 25)
   substring_EUR = substr($6, 96, 24)
   substring_SAS = substr($6, 120, 32)
   substring_ALL = substr($6, 2, 150)

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

    ## (substr($6, 2) removes the first 1 from the sample string (that is the query sample)

   print $1, $2, $3, $4, $5, substr($6, 2), $7, $8, $9, $10, Pop_AFR, Pop_AMR, Pop_EAS, Pop_EUR, Pop_SAS, Pop_ALL, Allele_AFR, Allele_AMR, Allele_EAS, Allele_EUR, Allele_SAS, Allele_ALL, GT_homWT, GT_het, GT_homVAR
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


# $23 - GT_homWT 
# $24 - GT_het 
# $25 - GT_homVAR

## This calculates p for HWE (based on 150 samples) - WT alle freq ((2*homWT) + het)/(2* (homWT + het + homVar))
## Checks for demoninator of zero (if SV is counted but genotype is ./. for a 1KGP sample)

awk '{
  if ($23 + $24 + $25 != 0) {
    result = ((2 * $23) + $24) / (2 * ($23 + $24 + $25));
  } else {
    result = 0;
  }
  print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, result
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp3.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp4.txt

##This calculates q for HWE (based on 150 samples) - (1-p)
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, (1-$26)}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp4.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp5.txt

##This calculates np^2, 2npq, and nq^2 for HWE (based on 150 samples)
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, (150*$26*$26), (300*$26*$27), (150*$27*$27)}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp5.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp6.txt

#This calculates chi-squared value for HWE. (Can't divide by zero, so if expected is zero, chi-sq=NULL) 
awk '{
   if ($28 > 0 && $29 > 0 && $30 > 0) {
       result = ((($23 - $28) * ($23 - $28)) / $28) + ((($24 - $29) * ($24 - $29)) / $29) + ((($25 - $30) * ($25 - $30)) / $30);
       print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, result;
   } else {
       print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, "NULL";
   }
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp6.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp7.txt

## $11 - Pop_AFR 
## $12 - Pop_AMR 
## $13 - Pop_EAS
## $14 - Pop_EUR
## $15 - Pop_SAS
## $16 - Pop_ALL
## $17 - Allele_AFR
## $18 - Allele_AMR
## $19 - Allele_EAS
## $20 - Allele_EUR
## $21 - Allele_SAS
## $22 - Allele_ALL
## $23 - GT_homWT
## $24 - T_het
## $25 - GT_homVAR


#This lists all population counts/frequencies, allele counts/frequencies, and HWE info
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, ($11/51), $12, ($12/18), $13, ($13/25), $14, ($14/24), $15, ($15/32), $16, ($16/150), $17, ($17/102), $18, ($18/36), $19, ($19/50), $20, ($20/48), $21, ($21/64), $22, ($22/300), $23, $24, $25, $26, $27, $31}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp7.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp8.txt

#This determines if the SV is in HWE
awk '{
   if ($40 == "NULL") {
       $40 = "HWE_NULL";
   } else if ($40 > 3.84) {
       $40 = "HWE_FALSE";
   } else {
       $40 = "HWE_TRUE";
   }
   print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp8.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp9.txt

##This adds the absolute value of the length of the SV to the start position (DELs are listed as negative lengths in the vcf, but the range is from the start position plus that number of nucleotides)
awk 'BEGIN {OFS="\t"} { if ($5 == "INS") $3 = $2 + 1; else $3 = $2 + ($4 < 0 ? -1*$4 : $4); print }' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp9.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp10.txt

#This annotates each variant with gene info (+/-1kb to try to capture promoters, etc)
bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp10.txt -b $GENES > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp11.txt

#This sorts the file by gene name
sort -k44,44 "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp11.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp12.txt

#This joins the gene name with any OMIM annotations
join -t $'\t' -a 1 -1 44 -2 4 -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.44,2.6' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp12.txt "$OMIM_GENES" > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp13.txt

#This condolidates rows that are identical except for the gene and lists the genes together in the last column (ie. a large deletion that spans 50 genes will be consolidated to one row)
awk 'BEGIN{FS=OFS="\t"} {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8 FS $9 FS $10 FS $11 FS $12 FS $13 FS $14 FS $15 FS $16 FS $17 FS $18 FS $19 FS $20 FS $21 FS $22 FS $23 FS $24 FS $25 FS $26 FS $27 FS $28 FS $29 FS $30 FS $31 FS $32 FS $33 FS $34 FS $35 FS $36 FS $37 FS $38 FS $39 FS $40} {a[key]=(a[key]!="")?(a[key] ", " $41):(a[key] $41)} {b[key]=(b[key]!="")?(b[key] ", " $42):(b[key] $42)} END {for (k in a) print k, a[k], b[k]}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp13.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp14.txt

#This replaces all empoty fields with a . (important for next command)
awk 'BEGIN {FS = OFS = "\t"} {for (i = 1; i <= NF; i++) if ($i == "") $i = "."} 1' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp14.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp15.txt

#These commands annotates each SV with info (for each bed file), removes unnecessary columns, labels annotated regions as such, and sorts the file (necessary for bedtools intersect)

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp15.txt -b $CODING_REGIONS | \
cut -f1-42,44 | \
awk 'BEGIN {FS=OFS="\t"} $43 == -1 {$43 = "."} $43 != -1 && $43 != "." {$43 = "exonic"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp16.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp16.txt -b $CENTROMERES | \
cut -f1-43,45 | \
awk 'BEGIN {FS=OFS="\t"} $44 == -1 {$44 = "."} $44 != -1 && $44 != "." {$44 = "centromeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp17.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp17.txt -b $PERICENTROMERES | \
cut -f1-44,46 | \
awk 'BEGIN {FS=OFS="\t"} $45 == -1 {$45 = "."} $45 != -1 && $45 != "." {$45 = "pericentromeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp18.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp18.txt -b $TELOMERES | \
cut -f1-45,47 | \
awk 'BEGIN {FS=OFS="\t"} $46 == -1 {$46 = "."} $46 != -1 && $46 != "." {$46 = "telomeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp19.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp19.txt -b $STR | \
cut -f1-46,48 | \
awk 'BEGIN {FS=OFS="\t"} $47 == -1 {$47 = "."} $47 != -1 && $47 != "." {$47 = "STR"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp20.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp20.txt -b $VNTR | \
cut -f1-47,49 | \
awk 'BEGIN {FS=OFS="\t"} $48 == -1 {$48 = "."} $48 != -1 && $48 != "." {$48 = "VNTR"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp21.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp21.txt -b $SEGDUPS | \
cut -f1-48,50 | \
awk 'BEGIN {FS=OFS="\t"} $49 == -1 {$49 = "."} $49 != -1 && $49 != "." {$49 = "segdups"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp22.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp22.txt -b $REPEAT_MASKER | \
cut -f1-49,51 | \
awk 'BEGIN {FS=OFS="\t"} $50 == -1 {$50 = "."} $50 != -1 && $50 != "." {$50 = "repeat"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp23.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp23.txt -b $GAPS | \
cut -f1-50,52 | \
awk 'BEGIN {FS=OFS="\t"} $51 == -1 {$51 = "."} $51 != -1 && $51 != "." {$51 = "gap"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.txt -b $DEFRABB_HICONF  -f 1.0 | \
cut -f1-51,53 | \
awk 'BEGIN {FS=OFS="\t"} $52 == -1 {$52 = "."} $52 != -1 && $52 != "." {$52 = "hiconf"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp25.txt

## $1-$10 - SV info (Chr, Start_Pos, End_Pos, SV_Length, SV_Type, 1KGP_support, Genotype, Alt_reads, Ref_reads, Total_reads)
## $11-$34 - allele/pop freq ($34 = AlleleFreqAll)
## $35-$40 - HWE (GT_homWT, GThet, GT_homVar, p, q, HWE)
## $41-$52 - Annotations (gene, OMIM, exonic, cent, peri, tel, str, vntr, segdup, repeat, gap, hiconf)

#This prints the columns in the correct order
awk 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $34, $41, $42, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp25.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.txt

#This creates an output file with the correct column headers
echo -e "Chr\tStart_Pos\tEnd_Pos\tSV_Length\tSV_Type\t1KGP_support\tGenotype\tAlt_reads\tRef_reads\tTotal_reads\tAllele_Freq_ALL\tGenes\tOMIM\tExonic\tCentromeric\tPericentromeric\tTelomeric\tSTR\tVNTR\tSegdup\tRepeat\tGap\tHiConf\tPop_Count_AFR\tPop_Freq_AFR\tPop_Count_AMR\tPop_Freq_AMR\tPop_Count_EAS\tPop_Freq_EAS\tPop_Count_EUR\tPop_Freq_EUR\tPop_Count_SAS\tPop_Freq_SAS\tPop_Count_ALL\tPop_Freq_ALL\tAllele_Count_AFR\tAllele_Freq_AFR\tAllele_Count_AMR\tAllele_Freq_AMR\tAllele_Count_EAS\tAllele_Freq_EAS\tAllele_Count_EUR\tAllele_Freq_EUR\tAllele_Count_SAS\tAllele_Freq_SAS\tAllele_Count_ALL\tAllele_Freq_ALL\tGT_homWT\tGT_het\tGT_homVAR\tHWE-p\tHWE-q\tHWE" > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS.txt

#This adds the data to the output file
cat "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.txt >> "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS.txt

#This removes "?s" from OMIM annotations (hard to sort) - optional
sed -i 's/?//g' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS.txt

##This removes intermediate "temp" files 
find "$UNSOLVED_DIR" -type f -name "*temp*" -exec rm {} \;
find "$UNSOLVED_DIR" -type f -name "*jasmine_query_only*" -exec rm {} \;
rm -r output/

##This facilitates and ends the loop for the input VCF list
done < "$QUERY_VCFS"