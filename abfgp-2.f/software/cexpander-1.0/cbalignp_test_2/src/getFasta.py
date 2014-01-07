import MySQLdb

h1 = [ "localhost", "root" ]
h2 = [ "db1.ab.wurnet.nl", "edouard" ]

U = h2

dbA = MySQLdb.connect( host = U[0], user = U[1], db =  "ASgene_Athal" )
cursorA = dbA.cursor()

dbO = MySQLdb.connect( host = U[0], user = U[1], db = "ASgene_Osativa" )
cursorO = dbO.cursor()

dbE = MySQLdb.connect( host = U[0], user = U[1], db = "edouard" )
cursorE = dbE.cursor()

def getAnchors( Anumber, dbCursor ):
    cm = """
        SELECT Anchor1, Anchor2 FROM gT_Alignments WHERE
            Anumber = %s """ % ( Anumber )
    zzz = dbCursor.execute( cm )
    if zzz == 0:
        return "", ""
    else:
        return dbCursor.fetchall()[0]

def getSequence( Anchor, dbCursor ):
    cm = """
        SELECT Sequence FROM predicted_proteome WHERE
            name = '%s' """ % ( Anchor )
    zzz = dbCursor.execute( cm )
    if zzz == 0:
        return ""
    else:
        return dbCursor.fetchall()[0][0]


def Hfasta( Anumber, dest, dbCursorO, dbSpec1, dbSpec2 ):
    anchors = getAnchors( Anumber, dbCursorO )

    S1 = getSequence( anchors[0], dbSpec1 )
    S2 = getSequence( anchors[1], dbSpec2 )

    fout = open( dest, "w+t" )
    fout.write( ">ara\n" )
    fout.write( S1 + "\n\n" )
    fout.write( ">osa\n" )
    fout.write( S2 + "\n\n" )

    fout.close()

Hfasta( 1513, "TF.fasta", cursorE, cursorA, cursorO )
