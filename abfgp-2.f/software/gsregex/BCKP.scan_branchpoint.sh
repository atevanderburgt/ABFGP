# Find (degenerate) branch point sequence(s) in (multi)fasta sequences
# Uses the ScanForMatches algorith
# Creates GFF on STDOUT for the sequence(s) on STDIN
# according to Kupfer e.a. 2004; Introns and splicing elements in five diverse fungi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PATH_SFM_PAT="$DIR/sfmpat_branchpoint"

# full path to scan_for_matches; it must be provided as first argument to this script!
PATH_SFM_EXE=$1;

if [ "$PATH_SFM_EXE" == "" ];
then
    echo "usage: cat SINGLE_FASTA | $0 <PATH_TO_SCAN_FOR_MATCHES>"
    exit 2;
fi


while read LINE; do echo ${LINE}; done | $PATH_SFM_EXE $PATH_SFM_PAT | cat -E | sed 's/ \$/#/' | tr -d ">\n[]" | tr "$:,#" "   \n" | sed 's/[CTUY][TU][AGR]A[CTUY]$/1/i' | sed 's/.[TU].A.$/-1/i' | awk -F' ' '{ printf "%s\tBPregex\tBPregex\t%s\t%s\t%s\t+\t.\n", $1, $2, $3, $4 }'

