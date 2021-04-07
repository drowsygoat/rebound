BEGIN {
    print "S(L)AM converter started @ " strftime() > "/dev/stderr"
    OFS="\t"
    FS="\t"
}

# --- FUNCTIONS --- #
function tag_finder(tag,    i){
    if (substr($0,1,1)!="@"){
        for (i=1; i<=NF+1; i++){
            if (i==NF+1) {
                print "Can't find " tag ". Problem at record " $0 > "/dev/stderr"
                # exit 1
            }else if (substr($i,1,5)!=tag){
                continue
            }else{
                return $i
            }
        }
    }
}

function MD_extender(md,  out,a,b) {
    out=""
    # print md, "md in fucntion"
    match(md,/MD:Z:/)
    a=substr(md,RSTART+RLENGTH)

    while(match(a,/\^[A-Z]+/)) {
        a=substr(a,1,RSTART-1) "h" substr(a,RSTART+RLENGTH)
    }
    n = patsplit(a, b, /[0-9]+/, seps)
    for (i=1;i<=n;i++) {
        for(j=1;j<=b[i];j++) {
            out=out "."
        }
        if (seps[i] ~ /[ACTGN]/) {
            out=out seps[i]
        }else if (seps[i] ~ /h/) {
                out=out
        }
    }
    return out ####
}
function MM_extractor(read, ec,     a,b,j,i,k,out) {
    # returns only the matched+mismatched portion (out) of the original read (read) as an array cr
    out=""
    split(read, b, "")
    split(ec, a, "")
    j=1
    for (i in a) {
        # print a[i], i
        if (a[i] ~ /[I]/) {
            j++
            continue
        }else if (a[i] ~ /[S]/) {
            j++
            continue
        }else if (a[i] ~ /[DN]/) {
            continue
        }else if (a[i] ~ /[M]/) {
            out = out b[j]
            j++
            continue
        }else{
            print "Unknown character(s) in CIGAR string of read" $1 > "/dev/stderr"
        }
    }
    return out
}
function cigar_extender(cygaro,   out, b,i,j){
    # expands CIGAR string into a longer string (out) with a single letter per alignment position; works only for CIGARs with [SIDM] symbols
    out = ""
    n = patsplit(cygaro, b, /[^[:digit:]]/, seps)
    # this is looping through the CIGAR "segments" and gets their lengths (seps)
    for (i=0; i<=n-1; i++) {
        # this iteratively adds letter codes to each element of the CIGAR string being expanded
       for (j=1; j<=seps[i]; j++) {
           out=out b[i+1]
       }
   }
   return out
}

function MP_producer(MD_input, ec, read,      a, b){
    # outputs array (read_mis) with mismatched bases in the read, and their respective indices; if there are not gaps in the alignment, and if nothing was clipped from the read, these indices are identical to the indices of the mismatched reference bases; indexing starts from the start of the alignment;
    # arrays ref_mis and read_mis contain mismatched bases and their positions (with respect to start of the alignment); an offeset array is populated that contains offset values for consecutive pairs; in case of reads that have no I, D or S, the values in read_mis and ref_mis are identical
    sum_md=0; delete(ref_mis); delete(read_mis)
    match(MD_input,/MD:Z:/)
    a=substr(MD_input,RSTART+RLENGTH)
    while(match(a,/\^[A-Z]+/)) {
        a=substr(a,1,RSTART-1) "h" RLENGTH-1 "^" substr(a,RSTART+RLENGTH)
    }
    n=patsplit(a, b, /[[:alpha:]^]/, seps)
    for(i=1; i<=n; i++){
        if (seps[i] == ""){
            seps[i] = 0
        }
        if (b[i]=="^"){
            sum_md+=seps[i-1]
        }else if (b[i]=="h"){
            sum_md+=seps[i-1]
        }else{
            sum_md+=seps[i-1]+1
        }
        if (b[i]!="^" && b[i]!="h"){
            ref_mis[sum_md]=b[i]
        }
    }
    for (i in ref_mis){
        j=0
        k=0
        for (m=1; m<=length(ec); m++){
            if(substr(ec,m,1) == "M"){
                k++
                if (k==i*1){
                    break
                }
            }else if (substr(ec,m,1) == "S" || substr(ec,m,1) == "I") {
                j++
            }else if (substr(ec,m,1) == "D" || substr(ec,m,1) == "N") {
                j--
                k++ # only deletions need to be acounted for to match the positions in cigars
            }
                offset[i] = j
        }

        read_mis[i+j] = substr(read,i+j,1)

        mp = "MP:Z:"
        sep = ","
        cnt=1 # counter
        for (i in ref_mis) {
            j = offset[i]
            # if (ref_mis[i] == "A" && read_mis[i+j] == "A"){
            #     mp = mp "0:" i+j ":" i
            # }
            if (ref_mis[i] == "A" && read_mis[i+j] == "C"){
                mp = mp "1:" i+j ":" i
            }
            else if (ref_mis[i] == "A" && read_mis[i+j] == "G"){
                mp = mp "2:" i+j ":" i
            }
            else if (ref_mis[i] == "A" && read_mis[i+j] == "T"){
                mp = mp "3:" i+j ":" i
            }
            else if (ref_mis[i] == "A" && read_mis[i+j] == "N"){
                mp = mp "4:" i+j ":" i
            }
            else if (ref_mis[i] == "C" && read_mis[i+j] == "A"){
                mp = mp "5:" i+j ":" i
            }
            # else if (ref_mis[i] == "C" && read_mis[i+j] == "C"){
            #     mp = mp "6:" i+j ":" i
            # }
            else if (ref_mis[i] == "C" && read_mis[i+j] == "G"){
                mp = mp "7:" i+j ":" i
            }
            else if (ref_mis[i] == "C" && read_mis[i+j] == "T"){
                mp = mp "8:" i+j ":" i
            }
            else if (ref_mis[i] == "C" && read_mis[i+j] == "N"){
                mp = mp "9:" i+j ":" i
            }
            else if (ref_mis[i] == "G" && read_mis[i+j] == "A"){
                mp = mp "10:" i+j ":" i
            }
            else if (ref_mis[i] == "G" && read_mis[i+j] == "C"){
                mp = mp "11:" i+j ":" i
            }
            # else if (ref_mis[i] == "G" && read_mis[i+j] == "G"){
            #     mp = mp "12:" i+j ":" i
            # }
            else if (ref_mis[i] == "G" && read_mis[i+j] == "T"){
                mp = mp "13:" i+j ":" i
            }
            else if (ref_mis[i] == "G" && read_mis[i+j] == "N"){
                mp = mp "14:" i+j ":" i
            }
            else if (ref_mis[i] == "T" && read_mis[i+j] == "A"){
                mp = mp "15:" i+j ":" i
            }
            else if (ref_mis[i] == "T" && read_mis[i+j] == "C"){
                mp = mp "16:" i+j ":" i
            }
            else if (ref_mis[i] == "T" && read_mis[i+j] == "G"){
                mp = mp "17:" i+j ":" i
            }
            # else if (ref_mis[i] == "T" && read_mis[i+j] == "T"){
            #     mp = mp "18:" i+j ":" i
            # }
            else if (ref_mis[i] == "T" && read_mis[i+j] == "N"){
                mp = mp "19:" i+j ":" i
            }
            else if (ref_mis[i] == "N" && read_mis[i+j] == "A"){
                mp = mp "20:" i+j ":" i
            }
            else if (ref_mis[i] == "N" && read_mis[i+j] == "C"){
                mp = mp "21:" i+j ":" i
            }
            else if (ref_mis[i] == "N" && read_mis[i+j] == "G"){
                mp = mp "22:" i+j ":" i
            }
            else if (ref_mis[i] == "N" && read_mis[i+j] == "T"){
                mp = mp "23:" i+j ":" i
            }
            # else if (ref_mis[i] == "N" && read_mis[i+j] == "N"){
            #     mp = mp "24:" i+j ":" i
            # }
            if (cnt < length(ref_mis)){
                mp = mp ","
                cnt++
            }
        }

        if (length(ref_mis) == 0) {
            print "Something went wrong with generating MP tag" > "/dev/stderr"
            exit 1
        }
    }
}

function TCRA_producer(MM_extracted, md_extended,    i){
    delete count; delete mmms
    split(MM_extracted, mm_split, "")
    split(md_extended, md_split, "")
    for (i=1; i<=length(mm_split); i++) {
        if (md_split[i] == ".") {
            mmms[i] = mm_split[i] mm_split[i]
        }else{
            mmms[i] = md_split[i] mm_split[i]
        }
    }
    for (i in mmms) {
        if       (mmms[i] == "AA"){
            count["AA"]++
        }else if (mmms[i] == "AC"){
            count["AC"]++
        }else if (mmms[i] == "AG"){
            count["AG"]++
        }else if (mmms[i] == "AT"){
            count["AT"]++
        }else if (mmms[i] == "AN"){
            count["AN"]++
        }else if (mmms[i] == "CA"){
            count["CA"]++
        }else if (mmms[i] == "CC"){
            count["CC"]++
        }else if (mmms[i] == "CG"){
            count["CG"]++
        }else if (mmms[i] == "CT"){
            count["CT"]++
        }else if (mmms[i] == "CN"){
            count["CN"]++
        }else if (mmms[i] == "GA"){
            count["GA"]++
        }else if (mmms[i] == "GC"){
            count["GC"]++
        }else if (mmms[i] == "GG"){
            count["GG"]++
        }else if (mmms[i] == "GT"){
            count["GT"]++
        }else if (mmms[i] == "GN"){
            count["GN"]++
        }else if (mmms[i] == "TA"){
            count["TA"]++
        }else if (mmms[i] == "TC"){
            count["TC"]++
        }else if (mmms[i] == "TG"){
            count["TG"]++
        }else if (mmms[i] == "TT"){
            count["TT"]++
        }else if (mmms[i] == "TN"){
            count["TN"]++
        }else if (mmms[i] == "NA"){
            count["NA"]++
        }else if (mmms[i] == "NC"){
            count["NC"]++
        }else if (mmms[i] == "NG"){
            count["NG"]++
        }else if (mmms[i] == "NT"){
            count["NT"]++
        }else if (mmms[i] == "NN"){
            count["NN"]++
        }
    }
    if ($2 ~ /^129|137|161|163|65|73|97|99|0$/){
        tc = "TC:i:" count["TC"]+0
    }else if ($2 ~ /^113|145|147|153|177|81|83|89|16$/) {
        tc= "TC:i:" count["AG"]+0
    }
    sp = ","
    ra = "RA:Z:" count["AA"]+0 sp count["AC"]+0 sp count["AG"]+0 sp count["AT"]+0 sp count["AN"]+0 sp count["CA"]+0 sp count["CC"]+0 sp count["CG"]+0 sp count["CT"]+0 sp count["CN"]+0 sp count["GA"]+0 sp count["GC"]+0 sp count["GG"]+0 sp count["GT"]+0 sp count["GN"]+0 sp count["TA"]+0 sp count["TC"]+0 sp count["TG"]+0 sp count["TT"]+0 sp count["TN"]+0 sp count["NA"]+0 sp count["NC"]+0 sp count["NG"]+0 sp count["NT"]+0 sp count["NN"]+0 sp
}

function TCRA_producer_simple(read, ec,      a, count) {
    # use "simple" function also if reads were entirely clipped (cases of 100% overlap) == only S in CIGARS; this reads may actually be tossed completely, but they can be used for DE analysis when one wants to use exactly the same pool of reads for both, DE and SLAM
    # no mathes and mismatches == full clip
    if ($6 ~ /^[^M]+$/){
        tc_simple = "TC:i:0"
        ra_simple = "RA:Z:0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,"
    }else{
            # "clean" CIGAR case
        if ($6 ~ /^[^SDI]+$/) {
            a = $10
            # "typical" CIGAR cases with insertions, deletions and partial clipping
        }else{
            a = MM_extractor($10, ec)
        }
        delete(count)
        split(a, a_split, "")
        for (i in a_split){
            if       (a_split[i] == "A"){
                count["A"]++
            }else if (a_split[i] == "C"){
                count["C"]++
            }else if (a_split[i] == "G"){
                count["G"]++
            }else if (a_split[i] == "T"){
                count["T"]++
            }else if (a_split[i] == "N"){
                count["N"]++
            }
        }
        tc_simple = "TC:i:0"
        sp = ","
        ra_simple = "RA:Z:" count["A"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["C"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["G"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["T"]+0 sp 0 sp 0 sp 0 sp 0 sp 0 sp count["N"]+0 sp
    }
}

# ---- main script ---- #

 {
    # get MD tag field if present
    if ($0 ~ /^.*MD:Z:.*$/){
        md = tag_finder("MD:Z:")
    }

    # change YI to XI for slamdunk compatibility; keep in mind the values will be expressed as percentagesa fractions instead; only if BBmap YI tags are present
    if ($0 ~ /^.*YI:f:.*$/){
        gsub("YI:f:", "XI:f:", $0)
    }

    # skip header and modify @RG by adding DS. SM can be supplied as argument. To uses files normally produced by nf-core pipeline. Alternatively, STAR + Picard MarkDuplicates
    if ($1 ~ /^@/){
        if ($1 ~ /^@RG/){
            print $1, $2, "DS:{'sequenced':0,'mapped':0,'filtered':0,'mqfiltered':0,'idfiltered':0,'nmfiltered':0,'multimapper':0,'dedup':0,'snps':0,'annotation':'','annotationmd5':''}", $3
        }else{
               print $0
        }
    # skip unmapped
    }else if ($2 ~ /^101|133|165|69|77|141$/){
    # use "simple" function that skips making MP if no mismatches are present (MP would otherwise be empty); second "if" condition is duplicated in the TCRA_producer_simple() function, maybe removed at some point...
    }else if (gensub(/MD:Z:/, "", "g" , md) ~ /^[0-9]*(\^[ACTGN]+[0-9]+)*$/){
        ec=cigar_extender($6)
        TCRA_producer_simple($10, ec)
        print $0, tc_simple, ra_simple
    }else{
        ec=cigar_extender($6)
        MM_extracted=MM_extractor($10, ec)
        md_ext=MD_extender(md)
        MP_producer(md, ec, $10)
        TCRA_producer(MM_extracted, md_ext)
        print $0, tc, ra, mp
        }
    }
END{
    print "S(L)AM converter ended @ " strftime() > "/dev/stderr"
}
