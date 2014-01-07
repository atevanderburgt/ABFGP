import miscList
import sys

def nextfasta( fin ):
    header = ""
    sequence = ""

    while 1:
        otell = fin.tell()

        line = fin.readline()
        if not line:
            break
        qline = line.strip()
        if len( qline ) > 0:
            if qline[0] == ">":
                if header == "":
                    header = qline
                else:
                    fin.seek( otell )
                    break
            elif header != "":
                sequence+=qline
    return header, sequence


def load_sequences( source ):
    headers = []
    sequences = []

    print "loading fasta:", source

    fin = open( source )
    while 1:
        header, sequence = nextfasta( fin )
        if header == "":
            break

        headers.append( header )
        sequences.append( sequence )

    fin.close()

    print "\t", len( headers )

    return headers, sequences


def prep_launch_fasta( headers, sequences, dest ):

    print "building:", dest

    fout = open( dest, "w+t" )

    i = 0
    while i < len( headers ):
        j = i + 1
        while j < len( headers ):
            fout.write( headers[i] + "\n" )
            fout.write( sequences[i] + "\n\n" )
            fout.write( headers[j] + "\n" )
            fout.write( sequences[j] + "\n\n" )
            j+=1
        i+=1
    fout.close()

def prep_report_file( headers, sequences, dest ):
    uheaders = headers[:]
    uheaders.sort()

    sizes = [0] * len( uheaders )

    i = 0
    while i < len( headers ):
        F = miscList.FindInList_Binary( headers[i],
                                        uheaders )
        if F > -1:
            sizes[F] = len( sequences[i] )
        i+=1

    fout = open( dest, "w+t" )
    i = 0
    while i < len( uheaders ):
        P = [ uheaders[i] ]
        P.append( str( sizes[i] ) )
        fout.write( "\t".join( P ) + "\n" )
        i+=1
    fout.close()

def handle_all( source, dest_fasta, dest_report ):

    headers, sequences = load_sequences( source )

    print "buildin fasta:", dest_fasta

    prep_launch_fasta( headers,
                       sequences,
                       dest_fasta )

    print "building report file"
    prep_report_file( headers,
                      sequences,
                      dest_report )
    

if __name__ == "__main__":
    if len( sys.argv ) != 4:
        print "source, launch_fasta, report file"
        sys.exit()
    else:
        handle_all( sys.argv[1],
                    sys.argv[2],
                    sys.argv[3] )
        
