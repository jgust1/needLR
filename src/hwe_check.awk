#! /usr/bin/env awk

#expects no arguments

BEGIN{
    chr=1;
    pos=2;
    tempcol=3;
    ref=4;
    alt=5;
    svlen=6;
    svtype=7;
    
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
    if ($np2 == "NULL") {
       $np2 = "HWE_NULL";
    } else if ($np2> 3.84) {
       $np2 = "HWE_FALSE";
    } else {
       $np2 = "HWE_TRUE";
    }
    print $0
}
