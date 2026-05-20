#! /usr/bin/env awk

#expects start, end, step

BEGIN{
    OFS="\t"
}
{
    result = "";
    for (i=start; i<=end; i+=step) {
        if ($i ~ /1/) {
            result = result "1";
        } else {
            result = result "0";
        }
    }
    printf("%s%s\n",$0,result) #there's no tab separator here because rearrange.awk create a file with a dangling tabs
}
