### miscList.py just contains some misc- list-related fuctions ( basics )

def loadList( Fname, keepEmty ):
    """ just load the lines into a list """
    fin = open( Fname )
    lst = [] #empty list
    while 1:
        line = fin.readline()
        if not line:
            break
        qline = line.strip()
        if keepEmty == True:
            lst.append( qline )
        else:
            if len(qline) > 0:
                lst.append(qline)
    fin.close()
    return lst


def saveList( Fname, lst, aCh ):
    """ simply save the list """
    fout = open( Fname, "w+t" )
    for i in lst:
        fout.write( i + aCh )
    fout.close()

def uniqueList( lst ):
    """ sort and keep unique ones """
    lst.sort()
    rt = []
    if len( lst ) > 0:
        rt.append( lst[0] )

    if len( lst ) > 1:
        for i in range( 1, len( lst ) ):
            if lst[i] != lst[i - 1]:
                rt.append( lst[i] )
    return rt

def addListToList( lstTo, lstFrom ):
    """ just add a list to a list """
    for h in lstFrom:
        lstTo.append( h )
    return lstTo

def FindInList_Linear( S, lst ):
    """ Linear search ( only use when binairy search is not applicable """
    fnd = -1
    J = 0
    while J < len( lst ):
        if S == lst[J]:
            fnd = J
            break
        else:
            J = J + 1
    return fnd

def FindInList_Binary( S, lst ):
    """ make sure lst is SORTED ( much faster than linear search ) """
    if len( lst ) < 1:
        return -1 #update 26-1-2010
    fnd = -1
    if S < lst[0] or S > lst[ len(lst) - 1 ]:
        return -1
    l = 0
    r = len( lst ) - 1
    while l <= r:
        mid = ( l + r ) / 2
        if S > lst[mid]:
            l = mid + 1
        elif S < lst[mid]:
            r = mid - 1
        else:
            fnd = mid
            break
    return fnd

def FindInList_Binary2( S, lst ):
    """ make sure lst is SORTED ( much faster than linear search ) """
    if len( lst ) < 1:
        return -1 #update 26-1-2010
    fnd = -1
    if S < lst[0]: #or S > lst[ len(lst) - 1 ]:
        return -1, 0
    if S > lst[ -1 ]:
        return -1, len( lst )
    
    l = 0
    r = len( lst ) - 1
    while l <= r:
        mid = ( l + r ) / 2
        if S > lst[mid]:
            l = mid + 1
        elif S < lst[mid]:
            r = mid - 1
        else:
            fnd = mid
            break
    return fnd, l

def FindInList_SinL( S, lst ): # 20 - 9 - 2006 
    """ search linear allways shortest in longest """
    fnd = -1
    for k in range( len( lst ) ):
        if len(S) > len( lst[k] ):
            if S.find( lst[k] ) > -1:
                fnd = k
                break
        else:
            if lst[k].find( S ) > -1:
                fnd = k
                break
    return fnd



### interlist comparison

def similarListItems( list1, list2 ):
    """ returns list of common items in list1 and list2 """
    sims = []
    for i in list1:
        if FindInList_Linear( i, list2 ) > -1:
            sims.append( i )
    return sims

def highestItems( list1 ):
    s = list1[0]
    for k in list1:
        if k > s:
            s = k
    ll = []

    for i in range( len( list1 ) ):
        if list1[i] == s:
            ll.append( i )
    return ll

def listMaxSpecific( nrList, itList ):
    g = nrList[ itList[0] ]

    for k in range( len( itList ) ):
        if nrList[ itList[k] ] > g:
            g = nrList[ itList[0] ]

    return g

def listMinSpecific( nrList, itList ):
    g = nrList[ itList[0] ]

    for k in range( len( itList ) ):
        if nrList[ itList[k] ] < g:
            g = nrList[itList[0] ]

    return g

def printList( lst ):
    for item in lst:
        print item





        
    
