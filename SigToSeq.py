def List_minus(A,B):
    '''returns a sublist of A exdluding elements from list B'''
    K=[]
    for x in A:
        if x not in B:
            K.append(x)
    return K

def Adjacent_bary_simplex(bary_simplex):
    '''returns the list of barycentric simplices that respectively faces
    v,f and e of bary_simplex are glued to'''
    i=bary_simplex[0]
    j=bary_simplex[1]
    k=bary_simplex[2]
    l=bary_simplex[3]
    m=list({0,1,2,3}-{j,k,l})[0]
    return [(i,j,l,k), (i,m,k,l),(i,j,k,m)]

def Bary_simplices(i):
    '''returns a list of all the 24 barycentric simplices of
    the ith ideal tetrahedron'''
    L=[]
    for j in {0,1,2,3}:
        for k in {0,1,2,3}-{j}:
            for l in {0,1,2,3}-{j,k}:
                L.append((i,j,k,l))
    return L

def Gluing_image(t,i,j,k):
    '''returns the image of k under the face gluing
    map that sends the jth face of the ith
    ideal tetrahedron to the face it is paired with
    in the Regina triangulation t'''
    g=t.simplex(i).adjacentGluing(j).inverse()
    return g.pre(k)

def Glued_tet_index(t,i,j):
    '''returns the index of the ideal tetrahedron
    glued to the jth face of ith ideal tetrahedron
    in the Regina triangulation t'''
    return t.simplex(i).adjacentSimplex(j).index()           

def Des_seq_cusp_info(t):
    '''returns an orbifold destination sequence of the Regina tetrahedral triangulation t'''

    def Paired_bary_simplex(i,j,k,l):
        '''returns the barycentric simplex that is glued to the barycentric simplex
        (i,j,k,l) along the j face of (i,j,k,l)'''
        i1=Glued_tet_index(t,i,j)
        j1=Gluing_image(t,i,j,j)
        k1=Gluing_image(t,i,j,k)
        l1=Gluing_image(t,i,j,l)
        return (i1,j1,k1,l1)

    def Glued_bary_simplex_pair(bary_smplx):
        '''Given barycentric simplex bary_smplx
        this returns the set containing bary_smplx and the barycentric simplex
        that bary_smplx is glued to'''

        i=bary_smplx[0]
        j=bary_smplx[1]
        k=bary_smplx[2]
        l=bary_smplx[3]
        Paired_bary_simplex(i,j,k,l)
        return {(i,j,k,l),Paired_bary_simplex(i,j,k,l)}


    '''lists all Glued_bary_simplex_pair possible:it's a list of sets consisting of two 4-tuples'''
    Glued_pair_list=[]
    for i in range(t.size()):
        for x in Bary_simplices(i):
            if Glued_bary_simplex_pair(x) not in Glued_pair_list:
                Glued_pair_list.append(Glued_bary_simplex_pair(x))


    def Adjacent_pairs(simplex_pair):
        '''returns the list of four lists each representing the Glued_bary_simplex_pair
        that respectively v,f_0,f_1 and e faces of simplex_pair (a list of two 4-tuples) are glued to'''
        
        v1=Adjacent_bary_simplex(simplex_pair[1])[0]
        f1=Adjacent_bary_simplex(simplex_pair[1])[1]
        e1=Adjacent_bary_simplex(simplex_pair[1])[2]
        f0=Adjacent_bary_simplex(simplex_pair[0])[1]
        v_adjacent=[v1,list(Glued_bary_simplex_pair(v1)-{v1})[0]]
        e_adjacent=[e1,list(Glued_bary_simplex_pair(e1)-{e1})[0]]
        f0_adjacent=[f1,list(Glued_bary_simplex_pair(f1)-{f1})[0]]
        f1_adjacent=[list(Glued_bary_simplex_pair(f0)-{f0})[0],f0]
        L=[v_adjacent,f0_adjacent,f1_adjacent,e_adjacent]
        return L
    
    Current_simplex_pair=list(Glued_pair_list[0])
    '''Current_simplex_pair is the list representing Glued_bary_simplex_pair currently traversing'''

    Untraversed_pairs=List_minus(Glued_pair_list,[Glued_pair_list[0]])
    '''Untraversed_pairs is the list of Glued_bary_simplex_pair which are not traversed yet'''

    Current_adjacents=Adjacent_pairs(Current_simplex_pair)
    '''Current_adjacents is the list of four lists
    each representing the Glued_bary_simplex_pair glued to respectively
    v,f0,f1 and e faces of Current_simplex_pair'''
    
    Cumul_pairs=[Current_simplex_pair]
    '''Cumul_pairs is the cumulative list of all the Glued_bary_simplex_pair
    traversed as well as seen till the Glued_bary_simplex_pair
    e face of the Current_simplex_pair is glued to'''
    
    cusp=[(0,(Current_simplex_pair[0][0], Current_simplex_pair[0][2]), (Current_simplex_pair[1][0], Current_simplex_pair[1][2]))]
    '''If the Current_simplex_pair is [(i,j,k,l),(i',j',k',l')] then cusp=[(0,(i,k),(i',k')]'''
    
    Nontraversed_seen_pairs=[]
    '''Nontraversed_seen_pairs is the list of non-traversed but seen
    Glued_bary_simplex_pair till the Glued_bary_simplex_pair
    e face of the Current_simplex_pair is glued to'''

    Des_seq_block=[]
    '''Des_seq_block is part of the destination sequence for the Current_simplex_pair'''


    '''Setting up the lists Cumul_pairs, Nontraversed_seen_pairs, Des_seq_block and cusp in Step 1'''
    for x in Current_adjacents:
        if x !=Current_simplex_pair:
            Cumul_pairs.append(x)
            Nontraversed_seen_pairs.append(x)
            '''Glued_bary_simplex_pair glued to the faces of Current_simplex_pair and differet from Current_simplex_pair
            are added to the lists Cumul_pairs and Nontraversed_seen_pairs'''
        Des_seq_block.append(Cumul_pairs.index(x))
        '''Des_seq_block is the list of indices of the Glued_bary_simplex_pair that are glued to v,f_0,f_1 and e faces of
        the Current_simplex_pair'''
        cusp.append((Cumul_pairs.index(x), (x[0][0], x[0][2]), (x[1][0], x[1][2])))
        '''for the m-th seen Glued_bary_simplex_pair in the Cumul_pairs
        [(i,j,k,l),(i',j',k',l')], (m, (i,k),(i',k')) is added to cusp'''

    '''cumul_des_seq_blocks is the list of cumulative Des_seq_block
    of all the Glued_bary_simplex_pair traversed so far'''
    cumul_des_seq_blocks=[Des_seq_block]
                          
    '''Updating Current_simplex_pair, Untraversed_pairs,Current_adjacents,
    Nontraversed_seen_pairs, cusp, Des_seq_block and cumul_des_seq_blocks
    from Step 2 onwards'''
    for label in range(1,12*t.size()):
        Des_seq_block=[]
        if len(Nontraversed_seen_pairs)>0:
            Current_simplex_pair=Nontraversed_seen_pairs[0]
        else:
            Current_simplex_pair=list(Untraversed_pairs[0])
        Untraversed_pairs.remove(set(Current_simplex_pair))
        Current_adjacents=Adjacent_pairs(Current_simplex_pair)
        if len(Nontraversed_seen_pairs)>0:
            Nontraversed_seen_pairs.pop(0)
        for x in Current_adjacents:
            if x not in Cumul_pairs:
                Cumul_pairs.append(x)
                cusp.append((Cumul_pairs.index(x), (x[0][0], x[0][2]), (x[1][0], x[1][2])))
                Nontraversed_seen_pairs.append(x)
            Des_seq_block.append(Cumul_pairs.index(x))
        cumul_des_seq_blocks.append(Des_seq_block)

    K=[]
    for i in range(len(cumul_des_seq_blocks)):
        for j in range(4):
            K.append(cumul_des_seq_blocks[i][j])
    '''K is the destination sequence and it is obtained from cumul_des_seq_blocks by removing
    the sqaure brackets around every element (which are Des_seq_block) of cumul_des_seq_blocks'''
                
    return (K, cusp)

def Des_seq(t):
    '''returns the destination sequence of the orbifold triangulation
    of the Regina tetrahedral triangulation t'''
    return Des_seq_cusp_info(t)[0]

def Cusp_info(t):
    '''returns the list of all (m, (i,k),(i',k')) where [(i,j,k,l),(i',j',k',l')]
    is the m-th tetrahedron/Glued_bary_simplex_pair in the orbifold triangulation
    of the Regina tetrahedral manifold trianguation t'''
    return Des_seq_cusp_info(t)[1]


def Cusp_class(t, x):
    '''Cusp_class lists all tetrahedron-ideal vertex pair (i,k) which corresponds
    to the same cusp as the tetrahedron-ideal vertex pair x of the Regina tetrahedral triangulation t'''
    
    equiv_class=[x]
    image_known=[]
    for pair in equiv_class:
        if pair not in image_known:
            for l in {0,1,2,3}-{pair[1]}:
                if (Glued_tet_index(t,pair[0],l),Gluing_image(t,pair[0], l, pair[1])) not in equiv_class:
                    equiv_class.append((Glued_tet_index(t,pair[0],l),Gluing_image(t,pair[0], l, pair[1])))
                    image_known.append(pair)
    return equiv_class



def Mfd_cusp_classes(t):
    '''Mfd_cusp_classes returns the list of all distinct Cusp_class(t,x)
    where x is a tetrahedron-ideal vertex pair
    in the Regina tetrahedral triangulation t'''

    All_vert=[]
    for i in range(t.size()):
        for k in range(4):
            All_vert.append((i,k))
        
    cusps=[]
    pairs_indexed=[]
        
    for x in All_vert:
        if x not in pairs_indexed:
            cusps.append(Cusp_class(t,x))
            pairs_indexed+=Cusp_class(t,x)

    return cusps


def Which_cusp(t, tet_vertex):
    '''Which_cusp reurns the index of the cusp/Cusp_class containing tetrahedron-ideal vertex pair
    tet_vertex in Mfd_cusp_classes(t)'''

    r=-1
    for x in Mfd_cusp_classes(t):
        if tet_vertex in x:
            r=Mfd_cusp_classes(t).index(x)
            break
    return r
            

def Mfd_cusp_index(t,orb_tet_num):
    '''Mfd_cusp_index returns the index of the cusp in Mfd_cusp_classes(t) associated to
    the Glued_bary_simplex_pair in the orbifold triangulation of the
    Regina tetrahedral triangulation t indexed by orb_tet_num''' 
    r=-2
    for y in Cusp_info(t):
        if y[0]==orb_tet_num:
            r=Which_cusp(t, y[1])
            break
    return r      




