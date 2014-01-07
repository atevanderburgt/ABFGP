./configure
make
cd src
date
./cbalignp -i test2.fasta -R -P -S > results.out
kwrite results.out
date
cd ..
