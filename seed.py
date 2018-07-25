import sensor, node, fold

import random, time, os, subprocess, sys

RNA_PATH = '/Users/abonham/Downloads/RNAstructure'

sys.path.append(f'{RNA_PATH}/exe')
os.environ['DATAPATH']=f'{RNA_PATH}/data_tables'
import RNAstructure

class Seed:
    '''
    __init__() is the constructor for the class seed.

    Parameters:
        initData      <-- A List. The graph data obtained from the seed file (see the
                          comment for make_graph() to see how this data is structured).
        recNodeName   <-- A string. The name of the node which is to contain the
                          recognition sequence. As with all node names, it is a number.
                          eg. the string '1', '2', '3' etc.
        recSeq        <-- A string. The recognition sequence.
        bindingState  <-- An integer, 0 or 1, representing the binding state of the
                          recognition sequence. 0 is for DS, 1 is for SS.
        seedName      <-- A string. The name of the seed graph
        maxSensorSize <-- An integer. The maximum number of bases the user allows the
                          sensor to be.

    Returns:
        A Seed object

    '''

    def __init__(self, initData, recNodeName, recSeq, bindingState, seedName, maxSensorSize):
        self.name = seedName
        self.head = node.SSNode(None)
        self.nodes = {}
        self.recNodeName = recNodeName
        self.recSeq = recSeq
        self.bindingState = bindingState
        self.maxSensorSize = maxSensorSize
        self.make_graph(initData, self.head, self.nodes, recNodeName, recSeq)

    '''
        make_graph() uses the seed graph data to construct the graph. The data looks like
        this:
            [2 1 3 3 5,     <- This line represents the data for one node
            3 2 2,          
            5 2 4,
            4 5 7 7 0,
            7 4 4]

        Each node is represented by a number. Non zero evens corrispond to DSNodes, odds
        corrispond to SSNodes, and 0 corrispods to None. The first node is always
        represended by 1, which must corrispond to an SSNode. In cases where the first
        node is actually a DSNode, the SSNode '1' will simply have length 0 (ie no base
        pairs are in the node).

        The first number in the data for a particular node will be the nodes number. The
        next number is the number representing the node's progenitor. The following
        numbers are the other nodes which can be linked through this one. If the current
        node is an SSNode, there will only be one other number, the downstream DSNode
        of this SSNode. If the current node is a DSNode, then there will be three more
        numbers. The first is the link to the MidSSNode1, the next to MidSSNode2, the last
        is to the downstreamSSnode. (In all cases, the links are listed in the order in
        which they would be encountered if the DNA strand were traced in the 5' to 3'
        direction.) Examples are given below.

        2 1 3 3 5

        The first number is a 2, so this is a DSNode, we're calling it '2'.
        '2' has a progenitor (also called the upstream SSNode, attatched to the 5' end of
        the "leading segment"). We're calling it '1'.
        '2' has a midSSNode1 (the SSNode that will be attatched to the 3' end of the
        "leading segment"). We're calling this SSNode, '3'.
        '2' also has a midSSNode2 (the SSNode that will be attached to 5' end of the
        "lagging segment"). This SSNode is also called '3', so it's the same as midSSNode1.
        '2' has a downstream SSNode (attatched to the 5' end of the "lagging segment"),
        we're calling this SSNode '5'.


        3 2 2

        The first number is a 3, which is odd, so this is a SSNode, and we're calling it
        '3.'
        '3' has progenitor '2', so '3' is attached at it's 5' end to the DSNode called '2.'
        '3' has downstream DSNode '2', so '3' is attatched at it's 3' end to the DSNode
        called '2.' (We can see that '3' must be a loop.)

        There are a few more rules for the seed graph data:
            1) The first line of data must represent the first real node of the graph. (ie
               it must be the node in which the first base pair from the 5' end will be.)
            2) If the graph starts with a DSNode, its progenitor must be lited as '1'.
            3) If the graph starts with an SSNode, the first node must be named '1', and
               its progenitor must be '0.'


        Parameters:
            data         <-- A list. The data described above.
            current      <-- An object of the class Node. The node currently under
                             construction.
            nodes        <-- A dictionary of Node objects. This dictionary will eventually
                             contain all the nodes used in this seed graph.
            recNodeName  <-- A string. The name (which is also a key in the nodes
                             dictionary) of the node which will contain the recognition
                             sequence.
            recSeq       <-- A string. The recognition sequence.

        Returns:
            Nothing
    '''
    def make_graph(self, data, current, nodes, recNodeName, recSeq):
        firstNodeData = data[0].split()
        if int(firstNodeData[0]) % 2 == 0:
            current.set_length(0)
            nodes[firstNodeData[1]] = current
            nodes[firstNodeData[0]] = node.DSNode(current)
            current.set_downstreamDSNode(nodes[firstNodeData[0]])
            self.build_the_rest(data, nodes)
        else:
            nodes[firstNodeData[0]] = current
            prev = current
            current = node.DSNode(prev)
            prev.set_downstreamDSNode(current)
            nodes[firstNodeData[2]] = current
            self.build_the_rest(data[1:], nodes)

    '''
    build_the_rest() is an auxiliary function of make_graph.
    '''
    def build_the_rest(self, data, nodes):
        if not data:
            #The list "data" is empty, we are out of data
            return
        #"data" has at least one element
        currentLine = data[0].split()
        current = nodes[currentLine[0]]
        links = []
        for i, v in enumerate(currentLine[2:]):
            if v in nodes:
                links.append(nodes[v])
                if i == 2: #have to set prog. the prog this node was made with is actually a child
                    nodes[v].set_progenitor(current)
            elif v == '0':
                links.append(None)
                nodes[v] = None
            elif int(v) % 2 == 0: #New DSNode
                nodes[v] = node.DSNode(current)
                links.append(nodes[v])
            else:#new SSNode 
                nodes[v] = node.SSNode(current)
                links.append(nodes[v])
        current.set_links(links)
        self.build_the_rest(data[1:], nodes)

    '''
        build_sensor() first builds a 'Sensor' sequence using the 'Seed' of 'self'
        and an object of the 'random.Random' class. The sensor sequence is
        constructed. In some cases the seed graph may require the use of more bases
        than permitted by the user. In that case, this function returns None. Otherwise
        the sensors sequence is fed to RNAstructure's FoldSingleStrand. A 'Sensor' object
        is then constructed from the resulting data. This object is returned.

        Parameters:
            core     <-- An integer, the number of the processor using this function.
            version  <-- An integer, the number of the sensor being built on this seed.
                         (ie. if this is the 1047th sensor built by processor 3 from seed
                          graph number 4: core = 3 and version = 1047).
            tempdir  <-- A string, the name of the temporary directory where unafold
                         files should be stored
        Returns:
            a Sensor object
    '''
    def build_sensor(self, core, version, tempdir, baseSeq):
        rand = random.Random()
        self.generate_node_sizes(rand)
        self.populate_nodes(rand)
        seq = ''.join(self.get_sequence())
        seq = seq.upper()
        #some graphs may result in sequences of larger length than maxSensorSize set by user
        if len(seq) > self.maxSensorSize:
            return None
            
        (leadingRecDat, laggingRecDat) = self.nodes[self.recNodeName].get_rec_seq_data()
        
        # Create an equivalent RNSAstructure object and then create a virtual ct file as list

        rna_obj = RNAstructure.RNA.fromString(seq, backbone='dna')
        rna_obj.FoldSingleStrand(percent=15,window=0)
        
        ct_list = []
        for structNum in range(rna_obj.GetStructureNumber()):

            length = rna_obj.GetSequenceLength()
            free_energy = rna_obj.GetFreeEnergy(structNum+1) # 1 vs 0 indexing

            # Write the header line of the 'ct file'
            header = f'{length} dG = {free_energy} {structNum+1}'
            ct_list.append(header)

            # Loop through bases
            for base in range(rna_obj.GetSequenceLength()):
                base_num = base + 1
                base_type = rna_obj.GetNucleotide(base_num)
                base_back = base_num - 1
                if base_num == rna_obj.GetSequenceLength()+1:
                    base_forward = 0
                else:
                    base_forward = base_num + 1
                base_pair = rna_obj.GetPair(base_num,structurenumber=structNum+1)
                
                #Write the line of the 'ct file'
                line = f'{base_num} {base_type} {base_back} {base_forward} {base_pair} {base_num}'
                ct_list.append(line)
        
        sen = sensor.Sensor(ct_list, leadingRecDat, laggingRecDat, self.bindingState, self.name, baseSeq)

        return sen

    '''
        generate_node_sizes() semi-randomly determines the size of the sensor, based on
        this number, a size for each node which represets physical DNA is assigned. Each
        node is given a minimum length of three bases or three base-pairs, depending on
        node type. One of the nodes contains the recognition sequence, its minimum size
        is the length of that sequence. Based on these minimums and the sensor size, the
        number of 'used bases' is calculated and subtracted from the sensor size. This new
        number is the number of bases left which are then assigned, randomly, to nodes.
        This method of determinig node size, while slightly complex, avoids many issues
        of other methods which comprimize the impartiality of random node size selection
        because the sensor has a size limit.

        Parameters:
            rand   <-- a pointer to an object of the class 'random'

        Returns:
            Nothing
    '''
    def generate_node_sizes(self, rand):
        self.nodes[self.recNodeName].set_length(len(self.recSeq)) #min len of node with recSeq
        #print "rec seq node len is " + str(self.nodes[self.recNodeName].get_length())
        MAX_SIZE = self.maxSensorSize
        MIN_SIZE = 20
        size = rand.randint(MIN_SIZE, MAX_SIZE)

        MIN_NODE_SIZE = 3 #to allow for loop SSNodes? 

        realNodes = {}
        for n in self.nodes: #initializing "real" (ie rep. physical DNA) nodes to min size
            current = self.nodes[n]
            if current == None: #this is not a 'real' node
                continue
            length = current.get_length()
            
            if length == 0:# this is not a 'real' node
                continue
            
            if length == -1: #is empty
                if current.get_state() == fold.Fold.SEQ_STATE["DS"]: #is DS
                    realNodes[n] = (current, MIN_NODE_SIZE)
                    size -= MIN_NODE_SIZE*2 #DS node uses 2X the number of bps
                else: #is SS
                    realNodes[n] = (current, MIN_NODE_SIZE)
                    size -= MIN_NODE_SIZE
            else: #is not empty (ie. has recognition seq.)
                if current.get_state() == fold.Fold.SEQ_STATE["DS"]: #is DS
                    realNodes[n] = (current, length)
                    size -= 2*length #DS node uses 2X the number of bps
                else: # is SS
                    realNodes[n] = (current, length)
                    size -= length
                    
        keys = [n for n in realNodes] # a list of the 'key' names in the realNodes dict
        while size > 0: #increasing the size of random nodes until size limit is reached
            key = random.choice(keys)
            (current, length) = realNodes[key]
            if current.get_state() == fold.Fold.SEQ_STATE["DS"]: #is DS
                realNodes[key] = (current, length+1)
                size -=2 #DS node uses 2X the number of bps
            else: #is SS
                realNodes[key] = (current, length+1)
                size -=1

        for r in realNodes: #assigning the new sizes to the respective nodes
            (n, s) = realNodes[r]
            n.set_length(s)
            
        
    '''
        populate_nodes() populates the empty nodes with DNA bases (ie. A, C, T, or G)
        this method requires that all nodes, which are not None, have a length.

        Parameters:
            rand    <-- a pointer to an object of a class from the module 'random'
        Returns:
            Nothing
    '''
    def populate_nodes(self, rand):
        for n in self.nodes.values():
            if n == None:
                continue
            length = n.get_length()
            seq = []
            #if this is the node with the recognition sequence we treat it differently
            if n == self.nodes[self.recNodeName]: #is node with recognition sequence
                extra = n.get_length() - len(self.recSeq) #the length not required for the recSeq
                if extra != 0:
                    relLocRecSeq = rand.randint(1, extra) #the position of the recSeq in the node
                    n.set_relLocRecStart(relLocRecSeq)
                    #print "TADA: rel loc of rec seq is " + str(relLocRecSeq)
                    n.set_relLocRecEnd(extra - (relLocRecSeq-1)+1)
                    #print 'there are ' + str(extra - (relLocRecSeq -1)) + ' spots left.'
                    #print 'rec seq size is ' + str(len(self.recSeq))
                    #print 'node size is ' + str(n.get_length())
                    #print 'node is ' + str(n)
                    seq = self.generate_rand_DNA_string(relLocRecSeq-1, rand)
                    end = self.generate_rand_DNA_string(extra - (relLocRecSeq-1), rand)
                    seq.extend(self.recSeq)
                    seq.extend(end)
                else:
                    #print "TADA: relitive loc of rec seq is 1."
                    n.set_relLocRecStart(1)
                    n.set_relLocRecEnd(1)
                    seq = list(self.recSeq)
            else: #this node does not contain the recognition sequence
                seq = self.generate_rand_DNA_string(length, rand)
            n.set_seq(seq)

        
    '''
        generate_rand_DNA_string() generates a list of pseudo-randomly selected
        DNA bases (ie. A, C, T, or G) of a specified size.

        Parameters:
            size   <-- an integer, the size of the desired list
            rand   <-- a pointer to an object of the 'random' class

        Returns:
            a list of random DNA letters
    '''        
    def generate_rand_DNA_string(self, size, rand):
        if size ==0:
            return []
        return [ rand.choice(['A', 'T', 'C', 'G']) for i in range(0, size) ]

    '''
        get_sequence() returns the sequence represented by the populated nodes of
        the seed graph up to a given node. If the nodes are not populated, the function
        will return an empty sequence.

        Parameters:
            None
        Returns:
            A string. The sequence that was constructed using this seed graph. 

    '''
    def get_sequence(self):
        prev = self.head
        (sequence, current) = prev.get_seq_and_next_node(None, 0)

        while current != None:
            (seq, nxt) = current.get_seq_and_next_node(prev, len(sequence))
            sequence.extend(seq)
            prev       = current
            current    = nxt

        return sequence
        
