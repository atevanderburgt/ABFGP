# read single !line! of gff and output feature length
read -a gff;
fstart=${gff[3]}
fstop=${gff[4]}
expr $fstop - $fstart + 1
