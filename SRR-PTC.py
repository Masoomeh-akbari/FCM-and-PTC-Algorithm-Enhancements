# Last version SRR-PTC
#********************************

#  Prob_TC
"""*********************************************************************************"""
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
import copy
import pdb
import csv
class Walk:

    """Represents a signed directed walk in a bipolar weighted digraph"""

    def __init__(self, vertexMultiSet, edgeSet, length):
        self.vertexMultiSet = vertexMultiSet
        self.edgeSet = edgeSet
        self.length = length

    def getVertexMultiSet(self):
        """Returns the vertex sequence of the walk. Vertices appear in the
        order they were added into the walk (reverse order)."""
        return self.vertexMultiSet

    def getEdgeSet(self):
        """Returns the set of relevant arcs (arcs of weight not equal to 1) of the walk."""
        return self.edgeSet

    def getLength(self):
        """Returns the length (number of arcs traversed) of the walk."""
        return len(self.vertexMultiSet) - 1

    def getMultiplicity(self, a, b):
        """Returns the multiplicity of vertex a in the vertex multiset of the walk minus one copy of the vertex b."""
        count = 0
        d = copy.deepcopy(self.vertexMultiSet)
        d.remove(b)
        for i in range(0, len(d)):
            if (d[i] == a):
                count = count + 1
        return count

def arc_index_prob(n,A,B):
    """From input matrix M, construct the matrices of positive (A) and negative (B) weights
    when the input consists of ch (=1 or 2) matrices. Return the arc -> index map 
    and arc -> probability map."""

    arcToIndexMap = {}
    arcProbabilityList = {}
    index = 0
   
   
    for s in range(0, n):
        for t in range(0, n):
            if A[s][t] != 0:
                arc = s, t, 1
                arcToIndexMap[arc] = index
                arcProbabilityList[index] = A[s][t]
                index = index + 1
            if B[s][t] != 0:
                arc = s, t, -1
                arcToIndexMap[arc] = index
                arcProbabilityList[index] = B[s][t]
                index = index + 1

    return arcToIndexMap, arcProbabilityList

def writePTransitiveClosureToFile(outfilename, n, A, B):
    """ Write output matrices A and B to a file."""
    
    
    fo = open(outfilename, "wb")
    fo.write(bytes("Positive Matrix: \n", 'UTF-8'))
    for i in range(n):
        for j in range(n):
            fo.write(bytes(str(round(A[i][j], 6)), 'UTF-8'))
            if j != n - 1:
                fo.write(bytes(",", 'UTF-8'))
        fo.write(bytes("\n", 'UTF-8'))
    fo.write(bytes("Negative Matrix: \n", 'UTF-8'))
    for i in range(n):
        for j in range(n):
            fo.write(bytes(str(round(B[i][j], 6)), 'UTF-8'))
            if j != n - 1:
                fo.write(bytes(",", 'UTF-8'))
        fo.write(bytes("\n", 'UTF-8'))
    fo.close()

def isSubset(a, b): 
    """Return True if event b is a subevent of event a, and False otherwise.
    An event e is a list with at least 2 components: e[0] is the characterictic vector of
    the set of included arcs, and e[1] is the characteristic vector of the set of excluded arcs."""
    
    if (a[0] & b[0] == a[0]) and (a[1] & b[1] == a[1]):
        return True
    else:
        return False


def minWalkGenHelper(n, A1, WRepSign1, WRepSign2, l, s, t, i, sign):
    """Generate new (s,t)-walks of length l for WRepSign2 from (i,t)-walks of length l-1
    in WRepSign1 and an (s,i,sign)-arc

    n : number of vertices
    A1 : matrix of weights of sign 'sign'
    wRepSign1 : matrix of sets of representatives of minimal (s,t)-walks of sign 'Sign1'
    wRepSign2 : matrix of sets of representatives of minimal (s,t)-walks of sign 'Sign2'
    l : length of the walk generated
    s : source vertex
    t : target vertex
    i : second vertex of the new walks
    sign : sign of the first arc of the new walks; Sign2 = sign * Sign1"""

    m = len(WRepSign1[i][t])
    for k in range(m):
        walk = WRepSign1[i][t][k]
        if(walk.getLength() == l - 1 and walk.getMultiplicity(s, t) <= 1):
            vertexSet = copy.deepcopy(walk.getVertexMultiSet())
            vertexSet.append(s)
            edgeSet = copy.deepcopy(walk.getEdgeSet())
            """ this is a relevant arc"""
            if(A1[s][i] < 1):
                edge = s, i, sign
                edgeSet.add(edge)
            newalk = Walk(vertexSet, edgeSet, l)

            """ include newalk in WRepSign2 if it does not contain an existing walk"""
            """ if newalk is to be included in WRepSign2, delete all existing walks that contain newalk"""

            goon = True
            q = len(WRepSign2[s][t])
            j = 0
            while j < q and goon: 
                if WRepSign2[s][t][j].getEdgeSet().issubset(newalk.getEdgeSet()):
                    goon = False
                else:
                    if newalk.getEdgeSet().issubset(WRepSign2[s][t][j].getEdgeSet()):
                        WRepSign2[s][t].remove(WRepSign2[s][t][j])
                        q = q - 1
                        j = j - 1
                    j = j+1
            if goon: 
                WRepSign2[s][t].append(newalk)


def genRepOfMinWalksHelper(n, A, B, wRepPos, wRepNeg, l):
    """Generate representatives of minimal directed walks of length l

       n : number of vertices
       A : matrix of positive weights
       B : matrix of negative weights
       l : length of the walks
       wRepPos : matrix of sets of representatives of minimal (s,t,+1)-walks
       wRepNeg : matrix of sets of representatives of minimal (s,t,-1)-walks"""

    
    for s in range(0, n):
        for t in range(0, n):
            for i in range(0, n):
                """ first arc is positive """
                if A[s][i] != 0:
                    """ creating positive walks """
                    minWalkGenHelper(n, A, wRepPos, wRepPos, l, s, t, i, 1)
                    """ creating negative walks """
                    minWalkGenHelper(n, A, wRepNeg, wRepNeg, l, s, t, i, 1)

                """ first arc is negative """
                if B[s][i] != 0:
                    """ creating positive walks """
                    minWalkGenHelper(n, B, wRepNeg, wRepPos, l, s, t, i, -1)
                    """ creating negative walks """
                    minWalkGenHelper(n, B, wRepPos, wRepNeg, l, s, t, i, -1)


def genRepOfMinWalks(n, A, B):
    """Generate a set of representatives of all minimal directed (s,t,sigma)-walks,
    for all (s,t,sigma)

    n : number of vertices  
    A : matrix of positive weights   
    B : matrix of negative weights"""

    """ initialization """

    wRepPos = [[[] for i in range(n)] for i in range(n)]
    wRepNeg = [[[] for i in range(n)] for i in range(n)]
    
    for i in range(0, n):
        for j in range(0, n):

            if A[i][j] != 0:
                edgeSet = set()
                if A[i][j] < 1:
                    edge = i, j, 1
                    edgeSet.add(edge)
                walk = Walk([j, i], edgeSet, 1)
                """ vertices are added at the back, in reverse order"""
                wRepPos[i][j].append(walk)

            if B[i][j] != 0: 
                edgeSet = set()
                if B[i][j] < 1:
                    edge = i, j, -1
                    edgeSet.add(edge)
                walk = Walk([j, i], edgeSet, 1)
                wRepNeg[i][j].append(walk)

    for l in range(2, 2*n+1):
        genRepOfMinWalksHelper(n, A, B, wRepPos, wRepNeg, l)

    return wRepPos, wRepNeg


def boolClosure(n, c, d, arcToIndexMap, arcProbabilityList, epsilon):
    """Compute probabilistic transitive closure of a bipolar weighted digraph using the Boolean Algebra method.
       n = number of vertices
       c = matrix of representatives of minimal positive walks
       d = matrix of representatives of minimal negative walks
       epsilon = max error ( = 0 for the exact method) """ 

    pPos = [[0 for i in range(n)] for i in range(n)]
    pNeg = [[0 for i in range(n)] for i in range(n)]

    for s in range(n):
        for t in range(n):
            if epsilon == 0:
                pPos[s][t] = TCweight(c[s][t], arcToIndexMap, arcProbabilityList)
                pNeg[s][t] = TCweight(d[s][t], arcToIndexMap, arcProbabilityList)
            else:
                pPos[s][t] = approxTCweight(c[s][t], arcToIndexMap, arcProbabilityList, 0, epsilon)
                pNeg[s][t] = approxTCweight(d[s][t], arcToIndexMap, arcProbabilityList, 0, epsilon)

    return pPos, pNeg
    """ matrices of positive and negative weights, resp.,  of the transitive closure """


def TCweight(walkList, arcToIndexMap, arcProbabilityList):
    """Compute the probability of existence of at least one walk in the list WalkList"""

    eventsQueue = []
    """ create events corresponding to the walks in the list, and store them in a queue """
    """ event = [included_edges,excluded_edges,number_of_included_edges] """
    for i in range(0, len(walkList)):
        edgeset = list(walkList[i].getEdgeSet())
        incl = 0
        excl = 0
        """ construct the char. vector 'incl' of 'edgeset' """
        for j in range(0, len(edgeset)):
            incl = incl + pow(2, arcToIndexMap[edgeset[j]])
        event = []
        event.append(incl)
        event.append(excl)
        event.append(len(edgeset))
        eventsQueue.append(event)

    p = 0
    eventsQueue = sorted(eventsQueue, key=lambda e: e[2])
    eventsQueue.reverse()
    """ delete the last component of each event (number of included edges) """
    for i in range(len(eventsQueue)):
        eventsQueue[i].pop()
    """ isolate an event e with the smallest number of terms e[2] """
    while len(eventsQueue) != 0:
        e = eventsQueue.pop()
        p = p + computeEventProbability(e, arcProbabilityList)
        newEventsQueue = []
        for j in range(len(eventsQueue)):
            """ merge (take intersection of the negation of) the isolated event e with
            the j-th event 'e1' from the queue """
            e1 = eventsQueue[j]
            includedEdges = e[0]
            excludedEdges = e[1]
            i = 0
            while includedEdges != 0:
                """ find position i of the first included edge """
                while includedEdges % 2 == 0:
                    includedEdges = includedEdges >> 1
                    i = i + 1
                if pow(2, i) & e1[0] == 0:
                    """ included edge i is not an included edge of event e1 """
                    e2 = copy.deepcopy(e1)
                    if pow(2, i) & e1[1] == 0:
                        """ included edge i is not an excluded edge of event e1 """
                        e2[1] = e2[1] + pow(2,i)
                    newEventsQueue.append(e2)
                includedEdges = includedEdges >> 1
                i = i + 1
            i = 0
            while excludedEdges != 0:
                """ find position i of the first excluded edge """
                while excludedEdges % 2 == 0:
                    excludedEdges = excludedEdges >> 1
                    i = i + 1
                if pow(2, i) & e1[1] == 0:
                    """ excluded edge i is not an excluded edge of event e1 """
                    e2 = copy.deepcopy(e1)
                    if pow(2, i) & e1[0] == 0:
                        """ excluded edge i is not an included edge of event e1 """
                        e2[0] = e2[0] + pow(2,i) 
                    newEventsQueue.append(e2)
                excludedEdges = excludedEdges >> 1
                i = i + 1
            
        i=0        
        while i < len(newEventsQueue):
            j = i+1
            while j < len(newEventsQueue):
                if isSubset(newEventsQueue[i], newEventsQueue[j]):
                    newEventsQueue.remove(newEventsQueue[j])
                else:
                    if isSubset(newEventsQueue[j], newEventsQueue[i]):
                       newEventsQueue.remove(newEventsQueue[i])
                       i=i-1
                       j=len(newEventsQueue)
                    else:
                        j=j+1
            i=i+1
        eventsQueue = copy.deepcopy(newEventsQueue)
        eventsQueue.reverse()
    return p


def approxTCweight(walkList, arcToIndexMap, arcProbabilityList, track, epsilon):
    """Compute the probability of existence of at least one walk in the list to
    precision epsilon"""

    eventsQueue = []
    """ current upper bound for error """
    r = 0
    """ create events and add them to queue """
    for i in range(len(walkList)):
        edgeset = list(walkList[i].getEdgeSet())
        incl = 0
        excl = 0
        for j in range(len(edgeset)):
            incl = incl + pow(2, arcToIndexMap[edgeset[j]])
        event = []
        event.append(incl)
        event.append(excl)
        prob = computeEventProbability(event, arcProbabilityList)
        event.append(prob)
        eventsQueue.append(event)
        r = r + prob

    p = 0
    if track == 1:
        print("Upper bound for error: r=", r)
        print("Length of events queue: ", len(eventsQueue))
        
    while r > epsilon:
        """ isolate the event with largest probability """
        eventsQueue = sorted(eventsQueue, key=lambda e: e[2])
        e = eventsQueue.pop()
        p = p + e[2]
        newEventsQueue = []

        for j in range(len(eventsQueue)):
            """ merge (take intersection of the negation of) the isolated event e with
              the j-th event 'e1' from the queue """
            e1 = eventsQueue[j]
            includedEdges = e[0]
            excludedEdges = e[1]
            i = 0
            while includedEdges != 0:
                """ find position i of the first included edge """
                while includedEdges % 2 == 0:
                    includedEdges = includedEdges >> 1
                    i = i + 1
                e2 = copy.deepcopy(e1)
                if pow(2, i) & e1[0] == 0:
                    if pow(2, i) & e1[1] == 0: 
                        e2[1] = e1[1] + pow(2,i)
                        e2[2] = e1[2] * (1 - arcProbabilityList[i])
                    newEventsQueue.append(e2)
                includedEdges = includedEdges >> 1
                i = i + 1
            i = 0
            while excludedEdges != 0:
                """ find position i of the first excluded edge """
                while excludedEdges % 2 == 0:
                    excludedEdges = excludedEdges >> 1
                    i = i + 1
                e2 = copy.deepcopy(e1)
                if pow(2, i) & e1[1] == 0: 
                   if pow(2, i) & e1[0] == 0:
                       e2[0] = e1[0] + pow(2,i)
                       e2[2] = e1[2] * arcProbabilityList[i]
                   newEventsQueue.append(e2)
                excludedEdges = excludedEdges >> 1
                i = i + 1
                
        """ remove redundant events """
        i=0        
        while i < len(newEventsQueue):
            j = i+1
            while j < len(newEventsQueue):
                if isSubset(newEventsQueue[i], newEventsQueue[j]):
                    newEventsQueue.remove(newEventsQueue[j])
                else:
                    if isSubset(newEventsQueue[j], newEventsQueue[i]):
                       newEventsQueue.remove(newEventsQueue[i])
                       i=i-1
                       j=len(newEventsQueue)
                    else:
                        j=j+1
            i=i+1
        eventsQueue = copy.deepcopy(newEventsQueue)

        r = 0
        x = len(eventsQueue)
        for i in range(0, x):
            r = r + eventsQueue[i][2]
        if track == 1:
            print("Upper bound for error: r=", r)
            print("Length of events queue: ", len(eventsQueue))
        
    return p


def computeEventProbability(event, arcProbabilityList):
    """Compute the probability of an event."""
    
    probability = 1
    includedEdges = event[0]
    i = 0
    while includedEdges != 0:
        while includedEdges % 2 == 0:
            includedEdges = includedEdges >> 1
            i = i + 1
        probability = probability * arcProbabilityList[i]
        includedEdges = includedEdges >> 1
        i = i + 1
    i = 0
    excludedEdges = event[1]
    while excludedEdges != 0:
        while excludedEdges % 2 == 0:
            excludedEdges = excludedEdges >> 1
            i = i + 1
        probability = probability * (1 - arcProbabilityList[i])
        excludedEdges = excludedEdges >> 1
        i = i + 1
    return probability
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
"""*********************************************************************************"""
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@



def add_vertex(D,u):
    # add vertex u, but no new arcs, to digraph D
    n=len(D)+1
    N = [[0 for j in range(n)] for i in range(n)]

    for i in range(n):
        for j in range(n):
            if i < u and j < u:
                N[i][j] = D[i][j]

            elif i < u and j > u:
                N[i][j] = D[i][j-1]

            elif i > u and j < u:
                N[i][j] = D[i-1][j]

            elif i > u and j > u:
                N[i][j] = D[i-1][j-1]

    return N

#*************************************************************************************** 

def NegWalk_Existence(A,B):

    # Warshall's Algorithm for the existence of a negative directed (s,t,sigma)-walk,
    # for all (s,t,sigma)
    # A: matrix of positive weights
    # B: matrix of negative weights
    # WP: 0-1 matrix of positive walks (existence only)
    # WN: 0-1 matrix of negative walks (existence only)


    n=len(A)
    WP=[[0 for i in range(n)] for i in range(n)]
    WN=[[0 for i in range(n)] for i in range(n)]

    for s in range(n):
        for t in range(n):
            if A[s][t]!=0:
                WP[s][t]=1

            if B[s][t]!=0:
                WN[s][t]=1

    for k in range(n):
        for i in range(2):
             for s in range(n):
                for t in range(n):
                    if (WP[s][k]==1 and WP[k][t]==1) or (WN[s][k]==1 and WN[k][t]==1):
                         WP[s][t]=1

                    if (WP[s][k]==1 and WN[k][t]==1) or (WN[s][k]==1 and WP[k][t]==1):
                         WN[s][t]=1                        

    return WN


#****************************************************************************************


def Ind(v,D):
    # compute indegree of vertex v in the digraph with weighted adjacency matrix D
    ind=0
    for i in range(len(D)):
        if D[i][v]!=0:
            ind=ind+1
    return(ind)

#*******************************************************************************    
def Outd(v,D):
    # compute outdegree of vertex v in the digraph with weighted adjacency matrix D
    outd=0
    for i in range(len(D)):
        if D[v][i]!=0:
            outd=outd+1
    return(outd)

#****************************************************************************

def in_neighbour(u,A,B):

    # Find the last in-neighbour v1 of vertex u and the weight of arc (v1,u)
    # in the digraph with weighted adjacency matrix D

    v1=None

    for i in range(len(A)):
        
        if A[i][u]!=0:    
            v1=i
            return v1,A[v1][u]
            break
        elif B[i][u]!=0:    
            v1=i
            return v1,-B[v1][u]
            break
    
#*****************************************************************************

def out_neighbour(u,A,B):

    # Find the last out-neighbour v2 of vertex u and the weight of arc (u,v2)
    # in the digraph with weighted adjacency matrix D

 
    v2=None

    for i in range(len(A)):
        if A[u][i]!=0:
            v2=i
            return v2,A[u][v2]
            break
        elif B[u][i]!=0:    
            v2=i
            return v2,-B[u][v2]
            break



#***************************************************************************************
def getInputMatrix(ch, M):

    # Convert the input matrix M into the weighted adjacency matrices of positive (A)
    # and negative weight (B)

    if ch == 1:
        n=len(M)
        A = [[0 for i in range(n)] for i in range(n)]
        B = [[0 for i in range(n)] for i in range(n)]
        for s in range(n):
            for t in range(n):
                if M[s][t] > 0:
                    A[s][t] = M[s][t]

                elif M[s][t] < 0:
                    B[s][t] = -M[s][t]

    else:
        n=len(M)//2
        A = [[0 for i in range(n)] for i in range(n)]
        B = [[0 for i in range(n)] for i in range(n)]
        for s in range(n):
            for t in range(n):
                if M[s][t] != 0:
                    A[s][t] = M[s][t]

                if M[s+n][t] != 0:    
                    B[s][t] = M[s+n][t]

    return A,B



#*********************************************************************************

def Reduced_Digraph(A,B,u):

    # delete vertex u from the digraph with weighted adjacency matrices A (positive)
    # and B (negative)

    A1 = []
    for i, Ni in enumerate(A):
        if i != u:
            A1.append(Ni[:u]+Ni[u+1:])

    B1 = []
    for i, Ni in enumerate(B):
        if i != u:
            B1.append(Ni[:u]+Ni[u+1:])

    return A1, B1



#***************************************************************************************

def vertex_type(A,B,u):

    # Determine the type of vertex u in the digraph with weighted adjacency matrices
    # A (positive) and B (negative)
    # t=1: vertex of indegree 1 and outdegree 0
    # t=2: vertex of indegree 0 and outdegree 1
    # t=3: vertex of indegree 1 and outdegree 1 not in a negative cycle
    # t=0: none of the above

    Id1=Ind(u,A)
    Od1=Outd(u,A)
    Id2=Ind(u,B)
    Od2=Outd(u,B)
    t=0
    Id=Id1+Id2
    Od=Od1+Od2

    if Id==1:    
        if Od==0:        
            t=1
            
        elif Od==1:
            WN=NegWalk_Existence(A,B)
            
            if WN[u][u]==0:
                t=3
            
    elif Id==0 and Od==1:
        t=2       

    return t  



#************************************************************************************************


def Reduction(A,B):
    # Deleting all vertices of types 1, 2, and 3 (giving priority to the first two)
    # in a digraph with weighted adjacency matrices A (positive) and B (negative)
    # Return the reduced digraph and the stack of deleted vertices
    # (each vertex with additional info required to recover the transitive closure)

    Stack=[]
    Reduction_completed=False

    while len(A)>2 and Reduction_completed==False :
       
        save_vertex_type3=None
        type_of_u=0
        
        for u in range(len(A)):

            type_of_u=vertex_type(A,B,u)
                
            if type_of_u==1:
                vertex,prob=in_neighbour(u,A,B)
                Stack.append([u,1,[vertex,prob]])
                A,B = Reduced_Digraph(A,B,u)
                break
            elif type_of_u==2:
                vertex,prob=out_neighbour(u,A,B)
                Stack.append([u,2,[vertex,prob]])
                A,B = Reduced_Digraph(A,B,u)
                break

            elif type_of_u==3: #save the vertex to delete only if there are no type 1 or 2 vertices
                save_vertex_type3=u


        if (type_of_u in {0,3}) and save_vertex_type3!=None:

            u=save_vertex_type3
            type_of_u=3

            v1,prob1=in_neighbour(u,A,B)
            v2,prob2=out_neighbour(u,A,B)
            Stack.append([u,3,[v1,prob1],[v2,prob2]])

            if prob1>0 and prob2>0:                
               
                A[v1][v2]=A[v1][v2]+A[v1][u]*A[u][v2]-A[v1][v2]*A[v1][u]*A[u][v2]

            elif prob1>0 and prob2<0:
              
                B[v1][v2]=B[v1][v2]+A[v1][u]*B[u][v2]-B[v1][v2]*A[v1][u]*B[u][v2]


            elif prob1<0 and prob2>0:
              
                B[v1][v2]=B[v1][v2]+B[v1][u]*A[u][v2]-B[v1][v2]*B[v1][u]*A[u][v2]


            elif prob1<0 and prob2<0:
           
                A[v1][v2]=A[v1][v2]+B[v1][u]*B[u][v2]-A[v1][v2]*B[v1][u]*B[u][v2]
                

            A,B = Reduced_Digraph(A,B,u)
            

        elif type_of_u==0 and save_vertex_type3==None: 

            Reduction_completed=True

    return A,B,Stack



            
#*********************************************************************************
def SRR_PTC(A,B,epsilon):  
    
    A,B,Stack=Reduction(A,B) 
    # A and B are now the weighted adjacency matrices of the reduced digraph
    
    # calling Prob_TC
    n=len(A)
    arcToIndexMap, arcProbabilityList = arc_index_prob(n,A,B)
    c, d = genRepOfMinWalks(n, A, B)
    A, B=boolClosure(n, c, d, arcToIndexMap, arcProbabilityList,epsilon)
    
    # recovery algorithm
    while Stack!=[]:
        stored=Stack.pop()
        u=stored[0] #stored vertex
        A=add_vertex(A,u)
        B=add_vertex(B,u)
        if stored[1]==1: # type of stored vertex u
       
            v=stored[2][0] # in-neighbour
            p=stored[2][1] # weight of arc (v,u)
            if p>0:
                A[v][u]=p
                B[v][u]=B[v][v]*p
                for s in range(len(A)):
                    if s!=v and s!=u:
                        A[s][u]=A[s][v]*p
                        B[s][u]=B[s][v]*p
            elif p<0:
                p=abs(p)
                B[v][u]=p
                A[v][u]=B[v][v]*p
                for s in range(len(A)):
                    if s!=v and s!=u:
                        A[s][u]=B[s][v]*p
                        B[s][u]=A[s][v]*p
            
        elif stored[1]==2: # type of stored vertex u
           
            v=stored[2][0] # out-neighbour
            p=stored[2][1] # weight of arc (v,u)
            if p>0:
                A[u][v]=p
                B[u][v]=B[v][v]*p
                for t in range(len(A)):
                    if t!=v and t!=u:
                        A[u][t]=p*A[v][t]
                        B[u][t]=p*B[v][t]
            if p<0:
                p=abs(p)
                B[u][v]=p
                A[u][v]=B[v][v]*p
                for t in range(len(A)):
                    if t!=v and t!=u:
                        A[u][t]=p*B[v][t]
                        B[u][t]=p*A[v][t]
            
        elif stored[1]==3: # type of stored vertex u

            v1=stored[2][0] # in-neighbour of u
            p1=stored[2][1] # weight of arc (v1,u)
            v2=stored[3][0] # out-neighbour of u
            p2=stored[3][1] # weight of arc (u,v2)
           
                
            if p1>0:
                A[v1][u]=p1
                B[v1][u]=B[v1][v1]*p1
                for s in range(len(A)):
                    if s!=v1 and s!=u:
                        A[s][u]=A[s][v1]*p1
                        B[s][u]=B[s][v1]*p1
            elif p1<0:
                p1=abs(p1)
                B[v1][u]=p1
                A[v1][u]=B[v1][v1]*p1
                for s in range(len(A)):
                    if s!=v1 and s!=u:
                        A[s][u]=B[s][v1]*p1
                        B[s][u]=A[s][v1]*p1
            
            if p2>0:
                A[u][v2]=p2
                B[u][v2]=B[v2][v2]*p2
                for t in range(len(A)):
                    if t!=v2 and t!=u:
                        A[u][t]=p2*A[v2][t]
                        B[u][t]=p2*B[v2][t]
            if p2<0:
                p2=abs(p2)
                B[u][v2]=p2
                A[u][v2]=B[v2][v2]*p2
                for t in range(len(A)):
                    if t!=v2 and t!=u:
                        A[u][t]=p2*B[v2][t]
                        B[u][t]=p2*A[v2][t]
            if v1==v2:
                A[u][u]=abs(p1*p2)                
            else:
                m1=A[v2][v1]
                m2=B[v2][v1]

                if p1*p2>0:
                    A[u][u]=abs(p1*p2*m1)

                elif p1*p2<0:
                    A[u][u]=abs(p1*p2*m2)
            
    return A,B
#**************************************************************************
""" main procedure """

print('\nCompute the probabilistic transitive closure of a bipolar weighted digraph.')
print('-------------------------------------------------------------------------')
print('\nInput Options:\n')
print('1. Input comprises of a single matrix.')
print('2. Input comprises of two matrices.\n')
print('Enter your choice of input:',)
ch = int(input())
print('\nEnter the number of vertices:',)
n = int(input())
print('\nEnter the input file path:',)
filename = input()
#*****************************
option = 1
while option == 1:
    print('\nEnter the rounding difference for the input matrix. It should be a decimal between 0 (inclusive) and 0.5 (exclusive). 0 results in no rounding.',)
    delta = float(input())
    
    """ create input matrices of arc weights """

    input1 = []
    with open(filename, 'r') as csvfile:
        inputReader = csv.reader(csvfile, delimiter=',')
        for row in inputReader:
            input1.append(row)

    for i in range(ch * n):
        for j in range(n):
            input1[i][j] = float(input1[i][j])
            if input1[i][j] < delta-1:
                input1[i][j] = -1
            elif input1[i][j] > -delta and input1[i][j] < delta:
                input1[i][j] = 0
            elif input1[i][j] > 1-delta:
                input1[i][j] = 1

    A,B=getInputMatrix(ch, input1)

    again = 1
    while again == 1:
        
        print('\nEnter the output file path for transitive closure:',)
        outfilename1 = input()
        print('\nAlgorithm Options:\n')
        print('1. Exact Boolean Algebra Approach.')
        print('2. Approximative Boolean Algebra Approach (recommended for all but the smallest input matrices).\n')
        print('Enter your choice of algorithm:',)
        ch2 = int(input())

        if ch2 == 2:
            print('\nEnter required precision:',)
            epsilon = float(input())    
        else:
            epsilon = 0

        """ compute transitive closure """

        e,f=SRR_PTC(A,B,epsilon)
        writePTransitiveClosureToFile(outfilename1,n, e, f)
        
        print('\nEnter 1 if you want to compute transitive closure again for the same rounding difference.',)
        again = int(input())

    print('\nEnter 1 if you want to compute transitive closure for a different rounding difference.',)
    option = int(input())

