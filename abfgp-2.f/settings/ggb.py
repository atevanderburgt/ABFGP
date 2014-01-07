################################################################################
# specify accessibility to GGB database for direct visualisation of results 
################################################################################

DB_CONNECTION = 'dbi:mysql'
# Specify database name, host,user & pwd for direct storing of results in GGB database
# The funtion gff2db() in gff/__init__.py will format a command line like this:
# perl /YOUR/PATH/TO/ABFGP/abfgp-2.0/gff/load_gff.pl --dsn dbi:mysql:<DB_NAME>:<DB_HOST> --user <DB_USER> --pass <DB_PASS> <GFF_FNAME> --fasta <FASTA_FNAME>;
DB_NAME       = 'XXXXXXXXXX'    # database name
DB_HOST       = 'XXXXXXXXXX'    # database host
DB_USER       = 'XXXXXXXXXX'    # database user
DB_PASS       = 'XXXXXXXXXX'    # database pass

