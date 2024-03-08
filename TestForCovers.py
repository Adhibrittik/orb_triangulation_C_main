def CuspSeqs(destseq):
    cusp_seqs = [] # a sequence of cusp sequences
    lorb = int(len(destseq)/4)
    if lorb != int(lorb):
        print('destseq error')
        return
    tetslist = []
    for i in range(lorb):
        tetslist.append(lorb - 1 - i)
    current_tet = tetslist.pop()
    this_cusp = [current_tet]
    if tetslist == []: cusp_seqs.append(this_cusp)
    while tetslist != []:
        for i in range(3):
            addtocusp = destseq[4*current_tet+i+1]
            if this_cusp.count(addtocusp) == 0:
                this_cusp.append(addtocusp)
                tetslist.remove(addtocusp)
        if this_cusp[len(this_cusp)-1] == current_tet:
            cusp_seqs.append(this_cusp)
            current_tet = tetslist.pop()
            this_cusp = [current_tet]
        else:
            for i in range(len(this_cusp)):
                if this_cusp[i] == current_tet:
                    current_tet = this_cusp[i+1]
                    break
        if tetslist == []: cusp_seqs.append(this_cusp)
    return cusp_seqs

def Covers(biggie,smalls):
    covering_dest_seqs = [] # a sequence of sequences, each encoding a cover
    lbig = int(len(biggie)/4) # the number of tets in the covering orbifold
    lsmall = int(len(smalls)/4) # the number of tets in the orbifold to be covered

    for i in range(lsmall):
        seq_i = [i] # encodes a (potential) cover sending tet 0 of biggie to
            # tet i of smalls. If the for loop below completes successfully,
            # seq_i will be a list of length lbig whose jth entry records the
            # destination of tet j of biggie
        goodflag = 1 # a kludge for breaking out of an inner then an outer loop
        counter = 0
        for j in range(lbig):
            for k in range(4):
                if biggie[4*j+k] == len(seq_i): # this is the case when the kth
                    # face of the jth tet abuts a tetrahedron whose destination
                    # has not yet been specified.
                    seq_i.append(smalls[4*seq_i[j]+k]) # in this case, the
                        # abutted tet's destination is the tet of smalls abutting
                        # the corresponding face of the jth tet's destination
                else:
                    if seq_i[biggie[4*j+k]] != smalls[4*seq_i[j]+k]: # if the
                        # kth face of the jth tet abuts a tet whose destination
                        # has been specified, and the specified destination does
                        # not match the tet of smalls abutting the kth face of
                        # the jth tet's destination, cut off the construction
                        goodflag = 0
                        break
                counter = counter + 1
            if goodflag == 0: break
        if counter == len(biggie): covering_dest_seqs.append(seq_i)

    return covering_dest_seqs


def CuspCovers(covering_dest_seq,cusp_seqs_up,cusp_seqs_down):
    cusp_pairs = []
    for up_seq in cusp_seqs_up:
        i = covering_dest_seq[up_seq[0]]
        for down_seq in cusp_seqs_down:
            if down_seq.count(i) != 0:
                cusp_pairs.append([up_seq,down_seq])
                break
    return cusp_pairs
