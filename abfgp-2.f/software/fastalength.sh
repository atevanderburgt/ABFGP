# output !single fasta! sequence length
if [ "$1" ]
then
  cat $1 | sed 1d | tr -d " \t\n\r" | wc -c | sed 's/^[ \t]*//'
else
  echo 0
fi
