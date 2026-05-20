#! /usr/bin/env awk

#expects ncontrols, ncohort, step (for cohort, not controls), controlcounts

BEGIN{
    nsubpops=split(controlcounts, pop_totals,  ",");
    chr=1;
    pos=2;
    tempcol=3;
    ref=4;
    alt=5;
    svlen=6;
    svtype=7;
    suppvecc=8;
    svgt=9;
    svvarreads=10;
    svrefreads=11
    totalreads=12
    sample_col_end=12
    control_start=sample_col_end+1
    control_end=control_start+ncontrols-1
    subpop_start=control_end+1
    if(ncohort>1){
      qsupportvec=8
      suppvecc=9
      cohort_start=10
      cohort_end=(cohort_start+(ncohort*step))-1
      control_start=cohort_end+1
      control_end=control_start+ncontrols-1
      cohort_pop=control_end+1
      cohort_allele=cohort_pop+1
      gt_coh_homwt=cohort_pop+2
      gt_coh_het=cohort_pop+3
      gt_coh_homvar=cohort_pop+4
      
      subpop_start=gt_coh_homvar+1
    }
    subpop_end=subpop_start+nsubpops-1
    pop_all=subpop_end+1
    subpop_allele_start=pop_all+1
    subpop_allele_end=subpop_allele_start+nsubpops-1
    allele_all=subpop_allele_end+1
    if(nsubpops<2){
        pop_all=control_end+1
        if(ncohort>1){
            pop_all=gt_coh_homvar+1
        }
        allele_all=pop_all+1
    }

    OFS="\t"
}
NR==1{
    gt_homwt=NF-8
    gt_het=NF-7
    gt_homvar=NF-6
    hwe_p=NF-5
    hwe_q=NF-4
    np2=NF-3
    npq2=NF-2
    nq2=NF-1
    chi2=NF
}
{
    if(ncohort>1){
        hap=ncohort*2
        popfreq=$cohort_pop/ncohort
        hapfreq=$cohort_allele/hap
        $cohort_pop=$cohort_pop"\t"popfreq
        $cohort_allele=$cohort_allele"\t"hapfreq
    }
    if(nsubpops>1){
        for(i in pop_totals){
            popindex=subpop_start+i-1;
            popfreq=$popindex/pop_totals[i];
            $popindex=$popindex"\t"popfreq;

            hap=pop_totals[i]*2;
            hapindex=subpop_allele_start+i-1;
            hapfreq=$hapindex/hap;
            $hapindex=$hapindex"\t"hapfreq
        }
    }
    pop_all_freq=$pop_all/ncontrols
    allele_all_freq=$allele_all/(ncontrols*2)
    $pop_all=$pop_all"\t"pop_all_freq
    $allele_all=$allele_all"\t"allele_all_freq
   print $0
}
