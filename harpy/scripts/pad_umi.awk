#! /usr/bin/awk -f

BEGIN {
    FS = "\t"

    # Precompute pad sequences for p=0..7
    pad[0]="TTTTTTT"; pad[1]="CCCCCC"; pad[2]="GGGGG"; pad[3]="AAAA"
    pad[4]="TTT";     pad[5]="CC";     pad[6]="GG";     pad[7]=""

    # Precompute quality pads (7-p) × 'I'
    qpad[0]="IIIIIII"; qpad[1]="IIIIII"; qpad[2]="IIIII"; qpad[3]="IIII"
    qpad[4]="III";     qpad[5]="II";     qpad[6]="I";     qpad[7]=""

    t = 0
}

# --- First file: sample_info.txt ---
NR==FNR {
    # Extract first token without index()
    rid = $1
    sub(/ .*/, "", rid)  # Remove everything after first space

    # Compute pad length p
    p = ($2 == -1 || $3 < 51 || $3 > 58) ? 7 : $3 - 51

    expID[++t] = "@" rid
    plen[t] = p
    next
}

# --- Second file: FASTQ (zcat) ---
{
    x = FNR % 4

    # Header line (x==1)
    if (x == 1) {
        t++
        p = plen[t]

        # Extract header token
        hdr = $0
        sub(/ .*/, "", hdr)

        if (hdr != expID[t]) {
            print "Read name mismatch at read " t "! Expected " expID[t] " but found R1=" hdr ". Exit!" > "/dev/stderr"
            exit 1
        }
        print
        next
    }

    # Sequence line (x==2)
    if (x == 2) {
        print (p == 6) ? pad[p] substr($0, 2) : pad[p] $0
        next
    }

    # '+' line (x==3)
    if (x == 3) {
        print
        next
    }

    # Quality line (x==0)
    print qpad[p] $0
}