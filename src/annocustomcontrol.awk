#! /usr/bin/env awk

#expects ncontrols, ncohort, step (for cohort, not controls), controlbreaks

BEGIN{
    chr=1;
    pos=2;
    tempcol=3;
    ref=4;
    alt=5;
    svlen=6;
    svtype=7;
    supportvec=8;
    svgt=9;
    svvarreads=10;
    svrefreads=11
    totalreads=12
    OFS="\t"
    if(ncohort>1){
      qsupportvec=8
      supportvec=9
      cohort_start=10
    }
    control_haps=ncontrols*2
    cohort_haps=ncohort*2
    cohort_end=supportvec+(ncohort*step)
    control_start=cohort_end+1
    control_end=control_start+ncontrols # these are for standard half open intervals

    nsubpops=split(controlbreaks, control_break_array, ",") # for 1kg this should be the value of AFR,AMR,EAS,EUR,SAS

}
{
  printroot=sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s",$chr, $pos, $tempcol, $ref, $alt, $svlen, $svtype)
  printsuffix=""
  # if there is a cohort, calculate cohort-level frequency statistics
  if( ncohort > 1){
    printroot=sprintf("%s\t%s\t%s", printroot, $qsupportvec, $supportvec)
    Pop_COHORT=0
    Allele_COHORT=0
    GT_coh_homWT=0
    GT_coh_het=0
    GT_coh_homVAR=0

    for (i=cohort_start; i<=cohort_end;i+=step){
      if( $i ~ /1/){
        Pop_COHORT++;
      }
      n=split($i,allele_array,"");
      for(j=1;j<=n;j++){
        if(allele_array[j] == "1"){
          Allele_COHORT++; }
        if ($i == "./." || $i == ".|." || $i == "0/0" || $i == "0|0"){
          GT_coh_homWT++; }
        if ($i == "0/1" || $i == "0|1" || $i == "1/0" || $i == "1|0"){
          GT_coh_het++; }
        if ($i == "1/1" || $i == "1|1") {
          GT_coh_homVAR++; }
      }
      for(k=0;k<step;k++){
        printroot=sprintf("%s\t%s",printroot,$(i+k))
      }
    }
    printsuffix=sprintf("%s\t%s\t%s\t%s\t%s\t%s",printsuffix,Pop_COHORT,Allele_COHORT,GT_coh_homWT,GT_coh_het,GT_coh_homVAR)
  } else {
    printroot=sprintf("%s\t%s\t%s\t%s\t%s\t%s", printroot, $supportvec, $svgt, $svvarreads, $svrefreads, $totalreads)
  }

  # always calculate the contrrol frequency statistics
  # if subpopulation indices were provided then calculate frequency stats for subpops
  if(nsubpops>1){
    subpop_pop_string=""
    subpop_allele_string=""
    istart=1
    jstart=control_start
    for(i in control_break_array){
      substring_control=substr($supportvec,istart,control_break_array[i])
      istart=istart+control_break_array[i]
      
      #count the number of 1s in the substring
      subpop_vec=gsub("1","",substring_control)
      if(subpop_vec==""){subpop_vec=0}
      subpop_pop_string=subpop_pop_string"\t"subpop_vec

      #allele counts
      jindex=jstart+control_break_array[i]
      allele_subpop=0
      for(j=jstart; j <= jindex; j++){
        n=split($j, allelearray,"");
        for(k=1;k<=n;k++){
          if(allelearray[k] == "1"){
            allele_subpop++;
          }
        }
      }
      jstart=jindex+1
      subpop_allele_string=subpop_allele_string"\t"allele_subpop
    }
    #subpop_pop_string will start with a tab
    #subpop_allele_string will start with a tab
    #contains POP_AFR, etc.
  }

  # compute control support ignoring sub populations
  substring_ALL = substr($supportvec, 1, ncontrols)
  Pop_ALL = gsub("1", "", substring_ALL)
  Allele_ALL = 0;

  GT_homWT = 0;
  GT_het = 0;
  GT_homVAR = 0;

  for (i = control_start; i <= NF; i++) {
    n = split($i, arr, "");
      for (j = 1; j <=n; j++) {
        if (arr[j] == "1") {
          Allele_ALL++;
      }
    }
    if ($i == "./." || $i == ".|." || $i == "0/0" || $i == "0|0") {
      GT_homWT++;
    }
    if ($i == "0/1" || $i == "0|1" || $i == "1/0" || $i == "1|0") {
      GT_het++;
    }
    if ($i == "1/1" || $i == "1|1") {
      GT_homVAR++;
    }
    printroot=sprintf("%s\t%s",printroot,$i)
  }

  if(nsubpops>1){
    printsuffix=sprintf("%s%s\t%s%s\t%s\t%s\t%s\t%s", printsuffix, subpop_pop_string, Pop_ALL, subpop_allele_string, Allele_ALL, GT_homWT, GT_het, GT_homVAR)
  }  else {
    printsuffix=sprintf("%s\t%s\t%s\t%s\t%s\t%s", printsuffix, Pop_ALL, Allele_ALL, GT_homWT, GT_het, GT_homVAR)
  }

  printf("%s%s\n", printroot, printsuffix) #expecting print stuffix to start with a tab
}
#the output format should now be the first 7 columns unaltered, the support vectors, the genotypes and reads dvdr for all cohort samples,
# cohort pop, cohort allele, gt hom, het, homvar for cohort
# then it should be all of the population gt columns, followed by subpopulation counts, pop_all, subpopulation allele counts, all allele,
# gt_homwt, gt_het,gt_homvar
