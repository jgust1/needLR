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
{
    if ($svtype == "INS"){ 
        $tempcol = $pos + 1;
    }else{
        $tempcol = $pos + ($svlen < 0 ? -1*$svlen : $svlen)
    }
    print $0
}
