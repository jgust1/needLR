#!/bin/bash

##Run from inside needLR_local

##Load environment with 
  #truvari v4.2.2 
  #bedtools v2.31.1 
  #bcftools v1.19 

##Default Truvari parameters
truvari_k="common"
truvari_s=50
truvari_S=10000000
truvari_B=50
truvari_r=2000
truvari_p=0
truvari_P=0.2
truvari_O=0.2

while getopts g:t:l:u:n:s:b:r:m:p: flag
do
    case "${flag}" in
        g) REF_GENOME=${OPTARG};;
        #t) TRUVARI_MERGED_CTRL_SET=${OPTARG};;
        #l) QUERY_VCFS=${OPTARG};;
        #u) MULTISAMPLE_VCF=${OPTARG};;
        #n) CTRL_SAMPLE_NAMES=${OPTARG};;
        #s) TOTAL_SAMPLES=${OPTARG};;
        b) PROBAND_VCF=${OPTARG};;
        #r) PARENTAL_VCF=${OPTARG};;
        m) MATERNAL_VCF=${OPTARG};;
        p) PATERNAL_VCF=${OPTARG};;
    esac
done

PROBAND_PREFIX=$(basename "$PROBAND_VCF" .vcf.gz)
MATERNAL_PREFIX=$(basename "$MATERNAL_VCF" .vcf.gz)
PATERNAL_PREFIX=$(basename "$PATERNAL_VCF" .vcf.gz)

QUERY_FILE_NAME=$(basename "$PROBAND_VCF")_needLR_3.4_trio
QUERY_SAMPLE_ID=$(basename "$QUERY_FILE_PATH" .vcf.gz)

UNSOLVED_DIR="needLR_output/$QUERY_FILE_NAME"
mkdir $UNSOLVED_DIR

#This is the Truvari merge of 450 1KGP samples
CTRL_MERGED="backend_files/UWONT_450_sniffles_2.5.2_T12_merge_sorted.vcf.gz"

##This is the list of the 1KGP sample names in the order they will be merged in Jasmine
CTRL_SAMPLE_NAMES="backend_files/ID_order_450.txt"

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

TOTAL_SAMPLES=450
AFR=147
AMR=61
EAS=81
EUR=68
SAS=93


#This preprocesses the proband VCF to be >=50bp, on full chromosomes, and FILTER=PASS
for vcf in $PROBAND_VCF; do
preprocessed_vcf=$(basename "$vcf" .gz)
bcftools view -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o "$UNSOLVED_DIR"/preprocessed_temp_"$preprocessed_vcf" "$vcf"
done

#This preprocesses the maternal VCF to be >=50bp, on full chromosomes, and FILTER=PASS
for vcf in $MATERNAL_VCF; do
preprocessed_vcf=$(basename "$vcf" .gz)
bcftools view -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o "$UNSOLVED_DIR"/preprocessed_temp_"$preprocessed_vcf" "$vcf"
done

#This preprocesses the paternal VCF to be >=50bp, on full chromosomes, and FILTER=PASS
for vcf in $PATERNAL_VCF; do
preprocessed_vcf=$(basename "$vcf" .gz)
bcftools view -i '(INFO/SVTYPE="BND") || (INFO/SVTYPE="INS" || INFO/SVTYPE="DEL" || INFO/SVTYPE="DUP" || INFO/SVTYPE="INV") && (INFO/SVLEN > 49 || INFO/SVLEN < -49)' -f PASS -r chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM -o "$UNSOLVED_DIR"/preprocessed_temp_"$preprocessed_vcf" "$vcf"
done

bgzip "$UNSOLVED_DIR"/preprocessed*

find "$UNSOLVED_DIR" -type f -name "preprocessed*gz" -exec tabix {} \;

#This creates the input list of VCFs for bcftools merge
echo "$UNSOLVED_DIR"/preprocessed_temp_$PROBAND_PREFIX*.vcf.gz > "$UNSOLVED_DIR"/list_temp.txt
echo "$UNSOLVED_DIR"/preprocessed_temp_$MATERNAL_PREFIX*.vcf.gz >> "$UNSOLVED_DIR"/list_temp.txt
echo "$UNSOLVED_DIR"/preprocessed_temp_$PATERNAL_PREFIX*.vcf.gz >> "$UNSOLVED_DIR"/list_temp.txt
echo "$CTRL_MERGED" >> "$UNSOLVED_DIR"/list_temp.txt

#bcftools merge
bcftools merge -m none --force-samples -l "$UNSOLVED_DIR"/list_temp.txt -Oz -o "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_with_UWONT450_temp.vcf.gz

tabix "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_with_UWONT450_temp.vcf.gz

#Truvari collapse
truvari collapse -i "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_bcftools_merged_with_UWONT450_temp.vcf.gz -o "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_temp.vcf -c "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_collapsed_temp.vcf -f $REF_GENOME -k $truvari_k -r $truvari_r -p $truvari_p -B $truvari_B -P $truvari_P -O $truvari_O --chain -s $truvari_s -S $truvari_S --passonly

##This parses the Truvari VCF output for relevant information
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/SVLEN\t%INFO/SVTYPE[\t%GT\t%DV\t%DR]\n' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_temp.vcf > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_temp.txt

#Only keep rows where column 7 is not "./."
#Column 7 is the genotype for the query sample
awk 'BEGIN {OFS="\t"} $7 != "./."' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_temp.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp.txt

##This removes BNDs from the analysis (too large to efficiently annotate)
cat "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_truvari_merged_query_only_temp.txt | grep -v BND > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp1.txt

##This parses for SVs <=10Mb
awk '($5 > -10000000 && $5 < 10000000) { print }' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp1.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp2.txt

##Rearrange the columns so that you keep the GT, DV, DR, (and total read depth) for any SVs in the query, and then order the GT outputs for the 1KGP samples

#$1 - Query Chr
#$2 - Query Pos
#$3 - Query REF allele
#$4 - Query ALT alele
#$5 - Query SV len
#$6 - Query SV type
#$7 - Query SV genotype
#$8 - Query SV variant reads
#$9 - Query SV reference reads
#$10 - Maternal genotype
#$11 - Maternal SV variant reads
#$12 - Maternal SV reference reads
#$13 - Paternal genotype
#$14 - Paternal SV variant reads
#$15 - Paternal SV reference reads

##Start column is the column after the parental sample info
START_COL=16
##End column is (total samples x 3) + (start column -2) 
END_COL=$(( (TOTAL_SAMPLES * 3) + 14 ))
STEP=3

awk -v start="$START_COL" -v end="$END_COL" -v step="$STEP" '
BEGIN {
  OFS="\t"
}
{
  printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
         $1, $2, $3, $4, $5, $6, $7, $8, $9, ($8 + $9), $10, $13)

  # Now print the variable columns in a loop
  for (i = start; i <= end; i += step) {
    printf("%s\t", $i)
  }
  printf("\n")
}
' "$UNSOLVED_DIR/$QUERY_FILE_NAME"_temp2.txt \
> "$UNSOLVED_DIR/$QUERY_FILE_NAME"_temp3.txt

##This creates a SUPP_VEC (like in Jasmine)

#$1 - Query Chr
#$2 - Query Pos
#$3 - Query REF allele
#$4 - Query ALT alele
#$5 - Query SV len
#$6 - Query SV type
#$7 - Query SV genotype
#$8 - Query SV variant reads
#$9 - Query SV reference reads
#$10 - Total reads
#$11 - Maternal genotype
#$12 - Paternal genotype
#$13-? - CTRl genotypes

#TEMPA is the total samples plus the last column before the first column of ctrl genotypes
TEMPA=$(( TOTAL_SAMPLES + 12 ))

awk -v tempa="$TEMPA" 'BEGIN {OFS="\t"} {
    result = "";
    for (i=13; i<=tempa; i++) {
        if ($i ~ /1/) {
            result = result "1";
        } else {
            result = result "0";
        }
    }
    print $0, result
}' "$UNSOLVED_DIR/$QUERY_FILE_NAME"_temp3.txt \
> "$UNSOLVED_DIR/$QUERY_FILE_NAME"_temp4.txt

##Rearrange columns so SUPP_VEC is $8
##START_COL = first column of control genotypes
START_COL=13
##END_COL = total samples + number of columns before ctrl genotypes
END_COL=$(( TOTAL_SAMPLES + 12 ))
SUPP_VECC=$(( TOTAL_SAMPLES + 13 ))
STEP=1

awk -v start="$START_COL" -v end="$END_COL" -v step="$STEP" -v suppvecc="$SUPP_VECC" 'BEGIN {OFS="\t"}
{printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
         $1, $2, "temp", $3, $4, $5, $6, $(suppvecc), $7, $8, $9, $10, $11, $12)

  # Print the range of columns in a loop
  for (i = start; i <= end; i += step) {
    printf("%s\t", $i)
  }
  printf("\n")
}
' "$UNSOLVED_DIR/$QUERY_FILE_NAME"_temp4.txt \
> "$UNSOLVED_DIR/$QUERY_FILE_NAME"_temp5.txt


##This calculates the population frequencies, allele frequencies, and allele types (for HWE)

TEMPAMR=$((AFR + 1))
TEMPEAS=$((TEMPAMR + AMR))
TEMPEUR=$((TEMPEAS + EAS))
TEMPSAS=$((TEMPEUR + EUR))

#14 is the number of columns before the control genotype array
TEMPAFRC=$((AFR + 14))

TEMPAMRB=$((TEMPAFRC + 1))
TEMPAMRC=$((TEMPAFRC + AMR))

TEMPEASB=$((TEMPAMRC + 1))
TEMPEASC=$((TEMPAMRC + EAS))

TEMPEURB=$((TEMPEASC + 1))
TEMPEURC=$((TEMPEASC + EUR))

TEMPSASB=$((TEMPEURC + 1))
TEMPSASC=$((TEMPEURC + SAS))

TEMPALL=$((TOTAL_SAMPLES + 14))

awk -v afr="$AFR" -v amr="$AMR" -v eas="$EAS" -v eur="$EUR" -v sas="$SAS" -v total_samples="$TOTAL_SAMPLES" -v tempamr="$TEMPAMR" -v tempeas="$TEMPEAS" -v tempeur="$TEMPEUR" -v tempsas="$TEMPSAS" -v tempafrc="$TEMPAFRC" -v tempamrb="$TEMPAMRB" -v tempamrc="$TEMPAMRC" -v tempeasb="$TEMPEASB" -v tempeasc="$TEMPEASC" -v tempeurb="$TEMPEURB" -v tempeurc="$TEMPEURC" -v tempsasb="$TEMPSASB" -v tempsasc="$TEMPSASC" -v tempall="$TEMPALL" 'BEGIN {OFS="\t"} {
   substring_AFR = substr($8, 1, afr)
   substring_AMR = substr($8, tempamr, amr)
   substring_EAS = substr($8, tempeas, eas)
   substring_EUR = substr($8, tempeur, eur)
   substring_SAS = substr($8, tempsas, sas)
   substring_ALL = substr($8, 1, total_samples)
   
   ## Count the number of 1s in the substring
   Pop_AFR = gsub("1", "", substring_AFR)
   Pop_AMR = gsub("1", "", substring_AMR)
   Pop_EAS = gsub("1", "", substring_EAS)
   Pop_EUR = gsub("1", "", substring_EUR)
   Pop_SAS = gsub("1", "", substring_SAS)
   Pop_ALL = gsub("1", "", substring_ALL)

  ## 15 is the first column of the control genotypes
   Allele_AFR = 0;
   Allele_AMR = 0;
   Allele_EAS = 0;
   Allele_EUR = 0;
   Allele_SAS = 0;
   Allele_ALL = 0;
   GT_homWT = 0;
   GT_het = 0;
   GT_homVAR = 0;
   for (i = 15; i <= tempafrc; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_AFR++;
       }
     }
   }
   for (i = tempamrb; i <= tempamrc; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_AMR++;
       }
     }
   }
   for (i = tempeasb; i <= tempeasc; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_EAS++;
       }
     }
   }
   for (i = tempeurb; i <= tempeurc; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_EUR++;
       }
     }
   }
   for (i = tempsasb; i <= tempsasc; i++) {
       n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_SAS++;
       }
     }
   }
   for (i = 15; i <= tempall; i++) {
      n = split($i, arr, "");
       for (j = 1; j <=n; j++) {
          if (arr[j] == "1") {
           Allele_ALL++;
       }
     }
   }
   for (i = 15; i <= tempall; i++) {
     if ($i == "./.") {
           GT_homWT++;
       }
     }
   for (i = 15; i <= tempall; i++) {
     if ($i == "0/1" || $i == "1/0") {
           GT_het++;
       }
     }
   for (i = 15; i <= tempall; i++) {
     if ($i == "1/1") {
          GT_homVAR++;
       }
     }

   print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, Pop_AFR, Pop_AMR, Pop_EAS, Pop_EUR, Pop_SAS, Pop_ALL, Allele_AFR, Allele_AMR, Allele_EAS, Allele_EUR, Allele_SAS, Allele_ALL, GT_homWT, GT_het, GT_homVAR
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp5.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp6.txt

##This matches the SUPP_VEC string to the sample list (and lists which 1KGP samples have the SV)
mapfile -t names < "$CTRL_SAMPLE_NAMES"

awk -v names="$(printf "%s\n" "${names[@]}")" '
BEGIN {
    split(names, nameArray, "\n");
    OFS="\t";
}
{
    # Get the binary string from column 8
    binary_string = $8;
    
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
    
    # Replace column 8 with the names string
    $8 = names_string;
    
    # Print the modified row
    print $0;
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp6.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp7.txt


#$27 - GT_homWT 
#$28 - GT_het 
#$29 - GT_homVAR

# This calculates p for HWE - WT alle freq ((2*homWT) + het)/(2* (homWT + het + homVar))
# Checks for demoninator of zero (if SV is counted but genotype is ./. for a 1KGP sample)

awk '{
  if ($27 + $28 + $29 != 0) {
    result = ((2 * $27) + $28) / (2 * ($27 + $28 + $29));
  } else {
    result = 0;
  }
  print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, result
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp7.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp8.txt

##This calculates q for HWE - (1-p)
awk '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, (1-$30)}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp8.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp9.txt

##This calculates np^2, 2npq, and nq^2 for HWE
TOTAL_HAPS=$((TOTAL_SAMPLES * 2))

awk -v total_haps="$TOTAL_HAPS" -v total_samples="$TOTAL_SAMPLES" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, (total_samples*$30*$30), (total_haps*$30*$31), (total_samples*$31*$31)}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp9.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp10.txt

#This calculates chi-squared value for HWE. (Can't divide by zero, so if expected is zero, chi-sq=NULL) 
awk '{
   if ($32 > 0 && $33 > 0 && $34 > 0) {
       result = ((($27 - $32) * ($27 - $32)) / $32) + ((($28 - $33) * ($28 - $33)) / $33) + ((($29 - $34) * ($29 - $34)) / $34);
       print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, result;
   } else {
       print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, "NULL";
   }
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp10.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp11.txt

## $14 - Pop_AFR 
## $15 - Pop_AMR 
## $16 - Pop_EAS
## $17 - Pop_EUR
## $18 - Pop_SAS
## $19 - Pop_ALL
## $20 - Allele_AFR
## $21 - Allele_AMR
## $22 - Allele_EAS
## $23 - Allele_EUR
## $24 - Allele_SAS
## $25 - Allele_ALL
## $26 - GT_homWT
## $27 - T_het
## $28 - GT_homVAR

AFRHAPS=$((AFR * 2))
AMRHAPS=$((AMR * 2))
EASHAPS=$((EAS * 2))
EURHAPS=$((EUR * 2))
SASHAPS=$((SAS * 2))

#This lists all population counts/frequencies, allele counts/frequencies, and HWE info
awk -v afr="$AFR" -v amr="$AMR" -v eas="$EAS" -v eur="$EUR" -v sas="$SAS" -v total_samples="$TOTAL_SAMPLES" -v afrhaps="$AFRHAPS" -v amrhaps="$AMRHAPS" -v eashaps="$EASHAPS" -v eurhaps="$EURHAPS" -v sashaps="$SASHAPS" -v total_haps="$TOTAL_HAPS" '{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, ($15/afr), $16, ($16/amr), $17, ($17/eas), $18, ($18/eur), $19, ($19/sas), $20, ($20/total_samples), $21, ($21/afrhaps), $22, ($22/amrhaps), $23, ($23/eashaps), $24, ($24/eurhaps), $25, ($25/sashaps), $26, ($26/total_haps), $27, $28, $29, $30, $31, $35}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp11.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp12.txt

#This determines if the SV is in HWE
awk '{
   if ($44 == "NULL") {
       $44 = "HWE_NULL";
   } else if ($44 > 3.84) {
       $44 = "HWE_FALSE";
   } else {
       $44 = "HWE_TRUE";
   }
   print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp12.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp13.txt

##This adds the absolute value of the length of the SV to the start position (DELs are listed as negative lengths in the vcf, but the range is from the start position plus that number of nucleotides)
awk 'BEGIN {OFS="\t"} { if ($7 == "INS") $3 = $2 + 1; else $3 = $2 + ($6 < 0 ? -1*$6 : $6); print }' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp13.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp14.txt

#This annotates each variant with gene info (+/-1kb to try to capture promoters, etc)
bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp14.txt -b $GENES > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp15.txt

#This sorts the file by gene name
sort -k47,47 "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp15.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp16.txt

#This joins the gene name with any OMIM annotations
join -t $'\t' -a 1 -1 47 -2 4 -o '1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,1.16,1.17,1.18,1.19,1.20,1.21,1.22,1.23,1.24,1.25,1.26,1.27,1.28,1.29,1.30,1.31,1.32,1.33,1.34,1.35,1.36,1.37,1.38,1.39,1.40,1.41,1.42,1.43,1.44,1.48,2.6' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp16.txt "$OMIM_GENES" > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp17.txt

#This condolidates rows that are identical except for the gene and lists the genes together in the last column (ie. a large deletion that spans 50 genes will be consolidated to one row)
awk 'BEGIN{FS=OFS="\t"} {key=$1 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7 FS $8 FS $9 FS $10 FS $11 FS $12 FS $13 FS $14 FS $15 FS $16 FS $17 FS $18 FS $19 FS $20 FS $21 FS $22 FS $23 FS $24 FS $25 FS $26 FS $27 FS $28 FS $29 FS $30 FS $31 FS $32 FS $33 FS $34 FS $35 FS $36 FS $37 FS $38 FS $39 FS $40 FS $41 FS $42 FS $43 FS $44} {a[key]=(a[key]!="")?(a[key] ", " $45):(a[key] $45)} {b[key]=(b[key]!="")?(b[key] ", " $46):(b[key] $46)} END {for (k in a) print k, a[k], b[k]}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp17.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp18.txt

#This replaces all empoty fields with a . (important for next command)
awk 'BEGIN {FS = OFS = "\t"} {for (i = 1; i <= NF; i++) if ($i == "") $i = "."} 1' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp18.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp19.txt

#These commands annotates each SV with info (for each bed file), removes unnecessary columns, labels annotated regions as such, and sorts the file (necessary for bedtools intersect)

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp19.txt -b $CODING_REGIONS | \
cut -f1-46,48 | \
awk 'BEGIN {FS=OFS="\t"} $47 == -1 {$47 = "."} $47 != -1 && $47 != "." {$47 = "exonic"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp20.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp20.txt -b $CENTROMERES | \
cut -f1-47,49 | \
awk 'BEGIN {FS=OFS="\t"} $48 == -1 {$48 = "."} $48 != -1 && $48 != "." {$48 = "centromeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp21.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp21.txt -b $PERICENTROMERES | \
cut -f1-48,50 | \
awk 'BEGIN {FS=OFS="\t"} $49 == -1 {$49 = "."} $49 != -1 && $49 != "." {$49 = "pericentromeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp22.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp22.txt -b $TELOMERES | \
cut -f1-49,51 | \
awk 'BEGIN {FS=OFS="\t"} $50 == -1 {$50 = "."} $50 != -1 && $50 != "." {$50 = "telomeric"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp23.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp23.txt -b $STR | \
cut -f1-50,52 | \
awk 'BEGIN {FS=OFS="\t"} $51 == -1 {$51 = "."} $51 != -1 && $51 != "." {$51 = "STR"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp24.txt -b $VNTR | \
cut -f1-51,53 | \
awk 'BEGIN {FS=OFS="\t"} $52 == -1 {$52 = "."} $52 != -1 && $52 != "." {$52 = "VNTR"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp25.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp25.txt -b $SEGDUPS | \
cut -f1-52,54 | \
awk 'BEGIN {FS=OFS="\t"} $53 == -1 {$53 = "."} $53 != -1 && $53 != "." {$53 = "segdups"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp26.txt -b $REPEAT_MASKER | \
cut -f1-53,55 | \
awk 'BEGIN {FS=OFS="\t"} $54 == -1 {$54 = "."} $54 != -1 && $54 != "." {$54 = "repeat"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp27.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp27.txt -b $GAPS | \
cut -f1-54,56 | \
awk 'BEGIN {FS=OFS="\t"} $55 == -1 {$55 = "."} $55 != -1 && $55 != "." {$55 = "Gap"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp28.txt

bedtools intersect -wa -wb -loj -a "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp28.txt -b $DEFRABB_HICONF  -f 1.0 | \
cut -f1-55,57 | \
awk 'BEGIN {FS=OFS="\t"} $56 == -1 {$56 = "."} $56 != -1 && $56 != "." {$56 = "hiconf"} 1' | \
sort | uniq > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp29.txt

## $1-$14 - SV info (Chr, Start_Pos, End_Pos, REF, ALT, SV_Length, SV_Type, 1KGP_support, Genotype, Alt_reads, Ref_reads, Total_reads, maternal gt, paternal)
## $15-$38 - allele/pop freq ($38 = AlleleFreqAll)
## $39-$44 - HWE (GT_homWT, GThet, GT_homVar, p, q, HWE)
## $45-$56 - Annotations (gene, OMIM, exonic, cent, peri, tel, str, vntr, segdup, repeat, CR/gaps, hiconf)

#This prints the columns in the correct order
awk -v query_sample_id="$QUERY_SAMPLE_ID" 'BEGIN {FS=OFS="\t"} {print $1, $2, $3, $4, $5, $6, $7, query_sample_id, $8, $9, $10, $11, $12, $13, $14, $38, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54, $55, $56, $15, $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $44}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp29.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp30.txt

# $1-$15 - SV info (Chr, Start_Pos, End_Pos, REF, ALT, SV_Length, SV_Type, query_id, 1KGP_support, Genotype, Alt_reads, Ref_reads, Total_reads, maternal gt, paternal)
# $16 - AlleleFreqAll
# $17-$28 - Annotations (gene, OMIM, exonic, cent, peri, tel, str, vntr, segdup, repeat, CR/gaps, hiconf)
# $29-$52 - allele/pop freq
# $53-$56 - HWE (GT_homWT, GThet, GT_homVar, HWE)


awk 'BEGIN {FS=OFS="\t"} {
    if ($14 == "./." && $15 == "./.") {
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "de novo", $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54, $55, $56
    } else {
        print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $13, $14, $15, "inherited", $16, $17, $18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $30, $31, $32, $33, $34, $35, $36, $37, $38, $39, $40, $41, $42, $43, $44, $45, $46, $47, $48, $49, $50, $51, $52, $53, $54, $55, $56
    }
}' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_temp30.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp1.txt

awk -F'\t' 'BEGIN { OFS="\t" }
{
  if (length($4) > 30000) {
    $4 = "[too long, see vcf]"
  }
  print
}
' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp1.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp2.txt

awk -F'\t' 'BEGIN { OFS="\t" }
{
  if (length($5) > 30000) {
    $5 = "[too long, see vcf]"
  }
  print
}
' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp2.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp3.txt

awk -F'\t' 'BEGIN { OFS="\t" }
{
  if (length($18) > 30000) {
    $11 = "[too long, see vcf]"
  }
  print
}
' "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp3.txt > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS_temp4.txt

#This creates an output file with the correct column headers
echo -e "Chr\tStart_Pos\tEnd_Pos\tREF\tALT\tSV_Length\tSV_Type\tQuery_ID\t1KGP_support\tGenotype\tAlt_reads\tRef_reads\tTotal_reads\tMaternal_genotype\tPaternal_genotype\tInheritance\tAllele_Freq_ALL\tGenes\tOMIM\tExonic\tCentromeric\tPericentromeric\tTelomeric\tSTR\tVNTR\tSegdup\tRepeat\tGap\tHiConf\tPop_Count_AFR\tPop_Freq_AFR\tPop_Count_AMR\tPop_Freq_AMR\tPop_Count_EAS\tPop_Freq_EAS\tPop_Count_EUR\tPop_Freq_EUR\tPop_Count_SAS\tPop_Freq_SAS\tPop_Count_ALL\tPop_Freq_ALL\tAllele_Count_AFR\tAllele_Freq_AFR\tAllele_Count_AMR\tAllele_Freq_AMR\tAllele_Count_EAS\tAllele_Freq_EAS\tAllele_Count_EUR\tAllele_Freq_EUR\tAllele_Count_SAS\tAllele_Freq_SAS\tAllele_Count_ALL\tAllele_Freq_ALL\tGT_homWT\tGT_het\tGT_homVAR\tHWE" > "$UNSOLVED_DIR"/"$QUERY_FILE_NAME"_RESULTS.txt

#This removes "?s" from OMIM annotations (hard to sort)
sed -i 's/?//g' "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp4.txt"

#This removes empty lines from the bottom of the txt file
sed -i -e :a -e '/^[[:space:]]*$/{$d;N;ba}' "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp4.txt"

#This adds the data to the output file
cat "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp4.txt" >> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.txt"

#This parses the RESULTS for SVs unique to the query (not in the control sample set)
awk -F'\t' 'BEGIN { OFS="\t" } NR == 1 || $17 == 0 { print }' "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.txt" > "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_unique.txt"

#This parses the RESULTS for SVs with an allele frequency <=0.01 in the control set
awk -F'\t' 'BEGIN { OFS="\t" } NR == 1 || $17 <=0.01 { print }' "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.txt" > "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_0.01.txt"


#This makes a vcf of the output
cat <<EOF > "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp.vcf"
##fileformat=VCFv4.2
##contig=<ID=chr1,length=248956422>
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
##contig=<ID=chr2,length=242193529>
##contig=<ID=chr20,length=64444167>
##contig=<ID=chr21,length=46709983>
##contig=<ID=chr22,length=50818468>
##contig=<ID=chr3,length=198295559>
##contig=<ID=chr4,length=190214555>
##contig=<ID=chr5,length=181538259>
##contig=<ID=chr6,length=170805979>
##contig=<ID=chr7,length=159345973>
##contig=<ID=chr8,length=145138636>
##contig=<ID=chr9,length=138394717>
##contig=<ID=chrM,length=16569>
##contig=<ID=chrX,length=156040895>
##contig=<ID=chrY,length=57227415>
##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the SV">
##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV length">
##INFO=<ID=SVTYPE,Number=1,Type=String,Description="SV type">
##INFO=<ID=Query_ID,Number=1,Type=String,Description="Query sample ID">
##INFO=<ID=Sample_support,Number=1,Type=String,Description="Samples in the 1000 Genomes Project samples (n=450) that share the SV with the query sample">
##INFO=<ID=Genotype,Number=1,Type=String,Description="Genotype of the SV in the query sample">
##INFO=<ID=Alt_reads,Number=1,Type=Integer,Description="Number of reads in the query sample supporting the SV">
##INFO=<ID=Ref_reads,Number=1,Type=Integer,Description="Number of reference reads in the query sample at the SV locus">
##INFO=<ID=Total_reads,Number=1,Type=Integer,Description="Total number of reads at the SV locus (in the query sample)">
##INFO=<ID=Maternal_genotype,Number=1,Type=Integer,Description="Parental genotype">
##INFO=<ID=Paternal_genotype,Number=1,Type=Integer,Description="Parental genotype">
##INFO=<ID=Inheritance,Number=1,Type=Integer,Description="Is the SV inherited from the parent?">
##INFO=<ID=Allele_Freq_ALL,Number=1,Type=Float,Description="Allele frequency of the SV in all 1KGP populations">
##INFO=<ID=Genes,Number=.,Type=String,Description="Genes that the SV overlaps (canonical hg38 coordinates, gencode)">
##INFO=<ID=OMIM,Number=.,Type=String,Description="OMIM phenotypes associated with any gene the SV overlaps (OMIM 8/2023)">
##INFO=<ID=Exonic,Number=1,Type=String,Description="If the SV overlaps exonic regions (canonical hg38 coordinates, gencode)">
##INFO=<ID=Centromeric,Number=1,Type=String,Description="If the SV overlaps centromeric regions (UCSC hg38)">
##INFO=<ID=Pericentromeric,Number=1,Type=String,Description="If the SV intersects with a pericentromeric region (+/-5Mb on either side of UCSC-defined centromere)">
##INFO=<ID=Telomeric,Number=1,Type=String,Description="If the SV intersects with a telomere (5Mb of either end of a chromosome)">
##INFO=<ID=STR,Number=1,Type=String,Description="If the SV intersects with a Short Tandem Repeat region (vamos original motifs, n=148)">
##INFO=<ID=VNTR,Number=1,Type=String,Description="If the SV intersects with a Variable Number Tandem Repeat region (vamos original motifs, n=148)">
##INFO=<ID=Segdup,Number=1,Type=String,Description="If the SV intersects with a segmental duplication (Genome in a Bottle v3.3)">
##INFO=<ID=Repeat,Number=1,Type=String,Description="If the SV intersects with a repeat region (UCSC hg38 repeat masker)">
##INFO=<ID=Gap,Number=1,Type=String,Description="If the SV intersects with a an hg38 gap region (UCSC hg38 mapping and sequencing: gap)">
##INFO=<ID=HiConf,Number=1,Type=String,Description="If the SV is fully contained within a high confidence region (Genome in a Bottle T2TQ100-V1.0_stvar)">
##INFO=<ID=Pop_Count_AFR,Number=1,Type=String,Description="How many 1KGP AFR ancestry samples have the SV (n=138)">
##INFO=<ID=Pop_Freq_AFR,Number=1,Type=String,Description="Frequency (%) of 1KGP AFR ancestry samples with SV (n=138)">
##INFO=<ID=Pop_Count_AMR,Number=1,Type=String,Description="How many 1KGP AMR ancestry samples have the SV (n=53)">
##INFO=<ID=Pop_Freq_AMR,Number=1,Type=String,Description="Frequency (%) of 1KGP AMR ancestry samples with SV (n=53)">
##INFO=<ID=Pop_Count_EAS,Number=1,Type=String,Description="How many 1KGP EAS ancestry samples have the SV (n=66)">
##INFO=<ID=Pop_Freq_EAS,Number=1,Type=String,Description="Frequency (%) of 1KGP EAS ancestry samples with SV (n=66)">
##INFO=<ID=Pop_Count_EUR,Number=1,Type=String,Description="How many 1KGP EUR ancestry samples have the SV (n=63)">
##INFO=<ID=Pop_Freq_EUR,Number=1,Type=String,Description="Frequency (%) of 1KGP EUR ancestry samples with SV (n=63)">
##INFO=<ID=Pop_Count_SAS,Number=1,Type=String,Description="How many 1KGP SAS ancestry samples have the SV (n=80)">
##INFO=<ID=Pop_Freq_SAS,Number=1,Type=String,Description="Frequency (%) of 1KGP SAS ancestry samples with SV (n=80)">
##INFO=<ID=Pop_Count_ALL,Number=1,Type=Integer,Description="Number of 1KGP individuals that carry the SV">
##INFO=<ID=Pop_Freq_ALL,Number=1,Type=Float,Description="Percent of 1KGP individuals that carry the SV">
##INFO=<ID=Allele_Count_AFR,Number=1,Type=Integer,Description="Frequency (%) of 1KGP AFR ancestry alleles with SV">
##INFO=<ID=Allele_Freq_AFR,Number=1,Type=Integer,Description="How many 1KGP AFR ancestry alleles have the SV">
##INFO=<ID=Allele_Count_AMR,Number=1,Type=Integer,Description="Frequency (%) of 1KGP AMR ancestry alleles with SV">
##INFO=<ID=Allele_Freq_AMR,Number=1,Type=Integer,Description="How many 1KGP AMR ancestry alleles have the SV">
##INFO=<ID=Allele_Count_EAS,Number=1,Type=Integer,Description="Frequency (%) of 1KGP EAS ancestry alleles with SV">
##INFO=<ID=Allele_Freq_EAS,Number=1,Type=Integer,Description="How many 1KGP EAS ancestry alleles have the SV">
##INFO=<ID=Allele_Count_EUR,Number=1,Type=Integer,Description="Frequency (%) of 1KGP EUR ancestry alleles with SV">
##INFO=<ID=Allele_Freq_EUR,Number=1,Type=Integer,Description="How many 1KGP EUR ancestry alleles have the SV">
##INFO=<ID=Allele_Count_SAS,Number=1,Type=Integer,Description="Frequency (%) of 1KGP SAS ancestry alleles with SV">
##INFO=<ID=Allele_Freq_SAS,Number=1,Type=Integer,Description="How many 1KGP SAS ancestry alleles have the SV">
##INFO=<ID=Allele_Count_ALL,Number=1,Type=Integer,Description="Number of alleles in the 1KGP cohort with the SV">
##INFO=<ID=Allele_Freq_ALL2,Number=1,Type=Float,Description="(2) Allele frequency of the SV in all 1KGP populations">
##INFO=<ID=GT_homWT,Number=1,Type=Integer,Description="Count of homozygous reference 1KGP genotypes">
##INFO=<ID=GT_het,Number=1,Type=Integer,Description="Count of heterozygous 1KGP genotypes for the SV">
##INFO=<ID=GT_homVAR,Number=1,Type=Integer,Description="Count of homozygous 1KGP genotypes for the SV">
##INFO=<ID=HWE,Number=1,Type=String,Description="Hardy-Weinberg Equilibrium result">
EOF

echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO' >> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp.vcf"

##Convert updated RESULTS.txt to VCF entries
awk -F"\t" '
BEGIN {
  OFS = "\t"
}
NR>1 {
  chrom = $1
  pos   = $2
  ref   = $4
  alt   = $5

  info = "END=" $3 \
         ";SVLEN=" $6 \
         ";SVTYPE=" $7 \
         ";Query_ID=" $8 \
         ";Sample_support=" $9 \
         ";Genotype=" $10 \
         ";Alt_reads=" $11 \
         ";Ref_reads=" $12 \
         ";Total_reads=" $13 \
         ";Maternal_genotype=" $14 \
         ";Paternal_genotype=" $15 \
         ";Inheritance=" $16 \
         ";Allele_Freq_ALL=" $17 \
         ";Genes=" $18 \
         ";OMIM=" $19 \
         ";Exonic=" $20 \
         ";Centromeric=" $21 \
         ";Pericentromeric=" $22 \
         ";Telomeric=" $23 \
         ";STR=" $24 \
         ";VNTR=" $25 \
         ";Segdup=" $26 \
         ";Repeat=" $27 \
         ";Gap=" $28 \
         ";HiConf=" $29 \
         ";Pop_Count_AFR=" $30 \
         ";Pop_Freq_AFR=" $31 \
         ";Pop_Count_AMR=" $32 \
         ";Pop_Freq_AMR=" $33 \
         ";Pop_Count_EAS=" $34 \
         ";Pop_Freq_EAS=" $35 \
         ";Pop_Count_EUR=" $36 \
         ";Pop_Freq_EUR=" $37 \
         ";Pop_Count_SAS=" $38 \
         ";Pop_Freq_SAS=" $39 \
         ";Pop_Count_ALL=" $40 \
         ";Pop_Freq_ALL=" $41 \
         ";Allele_Count_AFR=" $42 \
         ";Allele_Freq_AFR=" $43 \
         ";Allele_Count_AMR=" $44 \
         ";Allele_Freq_AMR=" $45 \
         ";Allele_Count_EAS=" $46 \
         ";Allele_Freq_EAS=" $47 \
         ";Allele_Count_EUR=" $48 \
         ";Allele_Freq_EUR=" $49 \
         ";Allele_Count_SAS=" $50 \
         ";Allele_Freq_SAS=" $51 \
         ";Allele_Count_ALL=" $52 \
         ";Allele_Freq_ALL2=" $53 \
         ";GT_homWT=" $54 \
         ";GT_het=" $55 \
         ";GT_homVAR=" $56 \
         ";HWE=" $57

  print chrom, pos, ".", ref, alt, ".", ".", info
}' "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.txt" \
>> "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp.vcf"

bcftools sort "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS_temp.vcf" > "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.vcf"

bgzip "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.vcf"
tabix "${UNSOLVED_DIR}/${QUERY_FILE_NAME}_RESULTS.vcf.gz"


##This removes intermediate "temp" files 
find "$UNSOLVED_DIR" -type f -name "*temp*" -exec rm {} \;
