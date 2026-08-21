#! /usr/bin/env awk

#expects no arguments

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
}
{
    if ($svtype == "INS"){ 
        $tempcol = $pos + 1;
    }else{
        $tempcol = $pos + ($svlen < 0 ? -1*$svlen : $svlen)
    }
    print $0
}
