#! /usr/bin/env awk

#expects indexField

BEGIN{
    FS=OFS="\t"
}
{
    # Construct composite key from columns 1–indexField
    key = $1
    for (i = 2; i <= indexField; i++) key = key FS $i

    # Track unique values per key for every column after the index field
    for ( i=indexField+1;i<=NF;i++){
        multiFields[i][key][$i] = 1
    }
}
END{
    for (k in multiFields[indexField+1]) { #this iterates over the key names, which are the full rows
        split(k, f, FS)
        for (i=1; i<=indexField; i++) printf "%s%s", f[i], (i<indexField ? OFS : "") # prints original line nicely
        for ( j in multiFields ){ #iterate over the field indices in multiField (numbers)
            fv="";
            for(v in multiFields[j][k]) fv = fv ? fv";"v : v;
            printf("%s%s", OFS,fv);
        }
        printf("\n");
    }   
}