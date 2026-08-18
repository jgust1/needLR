#! /usr/bin/env awk

#expects totsamples

BEGIN{
    chr=1;
    pos=2;
    needlrid=4;
    tempcol=3;
    ref=5;
    alt=6;
    svlen=7;
    svtype=8;
    OFS="\t"
    
    tothaps=totsamples*2
    OFS="\t"
}
NR==1{
    gt_homwt=NF-2
    gt_het=NF-1
    gt_homvar=NF
}
{
    if ($gt_homwt + $gt_het + $gt_homvar != 0) {
        hwe_p = ((2 * $gt_homwt) + $gt_het) / (2 * ($gt_homwt + $gt_het + $gt_homvar));
    } else {
        hwe_p = 0; #28
    }
    hwe_q=1-hwe_p; #29
    np2=totsamples*hwe_p*hwe_p; #30
    npq2=totsamples*2*hwe_p*hwe_q; #31
    nq2=totsamples*hwe_q*hwe_q; #32
    if(np2<1 || npq2 < 1 || nq2 < 1){
        chi2="NULL"
    }else{
        chi2=((($gt_homwt - np2)^2) / np2) + \
            ((($gt_het - npq2)^2) / npq2) + \
            ((($gt_homvar - nq2)^2) / nq2);
    }
    print $0,hwe_p,hwe_q,np2,npq2,nq2,chi2

}
