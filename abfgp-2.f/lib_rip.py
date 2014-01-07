"""
"""

def RIPindex(seq):
    """ """
    cnt_ca=0
    cnt_ta=0
    cnt_tg=0
    seq = seq.upper()
    offset=0
    while seq.find("CA",offset) >= 0:
        cnt_ca+=1
        offset = seq.find("CA",offset)+1
    offset=0
    while seq.find("TA",offset) >= 0:
        cnt_ta+=1
        offset = seq.find("TA",offset)+1
    offset=0
    while seq.find("TG",offset) >= 0:
        cnt_tg+=1
        offset = seq.find("TG",offset)+1
    # return ratio; max([x,1]) to prevent ZeroDivisionError
    return float(cnt_ta) / max([ cnt_ca + cnt_tg, 1])

# end of function RIPindex

