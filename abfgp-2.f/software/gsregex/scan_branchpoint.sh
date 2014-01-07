#!/usr/bin/env bash
# using getopts

help_prog()
{
  printf "Usage: cat SINGLE_FASTA | %s -e PATH_TO_SCAN_FOR_MATCHES [ -h ] \n" $0 >&2
  echo "Options: required"
  echo "   -e   full path to scan_for_matches executable"
  exit 1;
}

exe_sfm="";
while getopts 'e:h' OPTION
do
    case $OPTION in
    e)  exe_sfm="$OPTARG"
        ;;
    h)  help_prog
        exit 2
        ;;
    ?)  help_prog
        exit 2
        ;;
    esac
done

if [ "$exe_sfm" == "" ];
then
    help_prog;
    exit 2;
fi


# Find (degenerate) branch point sequence(s) in (multi)fasta sequences
# Uses the ScanForMatches algorith
# Creates GFF on STDOUT for the sequence(s) on STDIN
# according to Kupfer e.a. 2004; Introns and splicing elements in five diverse fungi
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PATH_SFM_PAT="$DIR/sfmpat_branchpoint"

while read LINE; do echo ${LINE}; done | $exe_sfm $PATH_SFM_PAT | cat -E | sed 's/ \$/#/' | tr -d ">\n[]" | tr "$:,#" "   \n" | sed 's/[CTUY][TU][AGR]A[CTUY]$/1/i' | sed 's/.[TU].A.$/-1/i' | awk -F' ' '{ printf "%s\tBPregex\tBPregex\t%s\t%s\t%s\t+\t.\n", $1, $2, $3, $4 }'

