while read LINE; do
   echo ${LINE} | awk '{ if (substr($0,1,1)==">") { print $0"@" } else { print $0"#" } }';
done | tr -d "\n" | tr -d "#" | sed "s/>/\n>/g" | sed "s/@/\n/g" | sed '/^$/d';
printf "\n";
exit 0;

#while read LINE; do
#   echo ${LINE};
#done | awk '{ if (substr($0,1,1)==">") { print $0"@" } else { print $0"#" } }' | tr -d "\n" | tr -d "#" | sed "s/>/\n>/g" | sed "s/@/\n/g" | sed '/^$/d';
#exit 0;
