#! /usr/bin/env awk

#expects step, ncontrols, ncohort

BEGIN{
    chr=1;
    pos=2;
    needlrid=3;
    ref=4;
    alt=5;
    svlen=6;
    svtype=7;
    svgt=8;
    svvarreads=9;
    svrefreads=10
    totalreads=11 #at end
    qend=10
    if(ncohort>1){
        qstart=8
        qend=qstart+(ncohort*3)-1
    }
    controlstart=qend+1
    controlend=controlstart+(ncontrols*3)-1
    genotypes=12
    OFS="\t"
}
{
    if(ncohort==1){
        tr=$svvarreads+$svrefreads
        $svrefreads=$svrefreads"\t"tr
    }else{
        for (i=qstart;i<=(qend-2);i+=step){
            tr=$(i+1)+$(i+2)
            $(i+2)=$(i+2)"\t"tr
        }
    }
    for(i=controlstart;i<=(controlend-2);i+=step){
        $(i+1)=""
        $(i+2)=""
    }
    $2=$2"\t."
    print $0
}
