# is argument a !single fasta! !dna! file? 
if [ "$1" ]
then
  twolines=$(grep ">" -m 2 -A 1 $1 | wc -l | sed 's/^[ \t]*//');
  if [ "$twolines" != 2 ]; then echo 0; exit 2; fi;
  # check for DNA sequence symbols
  nondnachars=$(cat $1 | sed 1d | sed 's/[nATGC \t\n]*//ig' | wc -w | sed 's/^[ \t]*//');
  if [ "$nondnachars" -eq "0" ]; then echo 1; else echo 0; fi;
else
  echo 0;
fi
