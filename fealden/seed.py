from __future__ import annotations

import random

from . import fold, node, sensor, structure


class Seed:

    """
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

    """

    def __init__(
        self,
        initData: list[str],
        recNodeName: str,
        recSeq: str,
        bindingState: int,
        seedName: str,
        maxSensorSize: int,
    ) -> None:
        """Initialize new Seed obj"""
        self.name = seedName
        self.head = node.SSNode(None)
        self.nodes: dict[str, node.DSNode | node.SSNode | None] = {}
        self.recNodeName = recNodeName
        self.recSeq = recSeq
        self.bindingState = bindingState
        self.maxSensorSize = maxSensorSize
        self.make_graph(initData, self.head, self.nodes, recNodeName, recSeq)

    def make_graph(
        self,
        data: list[str],
        current: node.SSNode | node.DSNode,
        nodes: dict[str, node.DSNode | node.SSNode | None],
        recNodeName: str,
        recSeq: str,
    ) -> None:
        """
        make_graph() uses the seed graph data to construct the graph. Data example:
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
        numbers. The 1st is the link to the MidSSNode1, the next to MidSSNode2, the last
        is to the downstreamSSnode. (In all cases, the links are listed in the order in
        which they would be encountered if the DNA strand were traced in the 5' to 3'
        direction.) Examples are given below.

        2 1 3 3 5

        The first number is a 2, so this is a DSNode, we're calling it '2'.
        '2' has a progenitor (also called the upstream SSNode, attached to the 5' end of
        the "leading segment"). We're calling it '1'.
        '2' has a midSSNode1 (the SSNode that will be attatched to the 3' end of the
        "leading segment"). We're calling this SSNode, '3'.
        '2' also has a midSSNode2 (the SSNode that will be attached to 5' end of the
        "lagging segment"). This SSNode is also called '3', so same as midSSNode1.
        '2' has a downstream SSNode (attatched to the 5' end of the "lagging segment"),
        we're calling this SSNode '5'.


        3 2 2

        The first number is a 3, which is odd, so this is a SSNode, and we're calling it
        '3.'
        '3' has progenitor '2', so '3' is attchd at it's 5' end to DSNode called '2.'
        '3' has downstream DSNode '2', so '3' is attatched at it's 3' end to the DSNode
        called '2.' (We can see that '3' must be a loop.)

        There are a few more rules for the seed graph data:
            1) The 1st line of data must represent the first real node of the graph. (ie
               it must be the node in which the 1st base pair from the 5' end will be.)
            2) If the graph starts with a DSNode, its progenitor must be lited as '1'.
            3) If the graph starts with an SSNode, the first node must be named '1', and
               its progenitor must be '0.'


        Parameters:
            data         <-- A list. The data described above.
            current      <-- An object of the class Node. The node currently under
                             construction.
            nodes        <-- A dictionary of Node objects. This dictionary will
                             contain all the nodes used in this seed graph.
            recNodeName  <-- A string. The name (which is also a key in the nodes
                             dictionary) of the node which will contain the recognition
                             sequence.
            recSeq       <-- A string. The recognition sequence.

        Returns:
            Nothing
        """

        firstNodeData = data[0].split()
        if int(firstNodeData[0]) % 2 == 0:
            current.set_length(0)
            nodes[firstNodeData[1]] = current
            nodes[firstNodeData[0]] = node.DSNode(current)
            assert not isinstance(current, node.DSNode)
            current.set_downstreamDSNode(nodes[firstNodeData[0]])
            self.build_the_rest(data, nodes)
        else:
            nodes[firstNodeData[0]] = current
            prev = current
            current = node.DSNode(prev)
            assert not isinstance(prev, node.DSNode)
            prev.set_downstreamDSNode(current)
            nodes[firstNodeData[2]] = current
            self.build_the_rest(data[1:], nodes)

    def build_the_rest(
        self, data: list[str], nodes: dict[str, node.SSNode | node.DSNode | None]
    ) -> None:
        """
        build_the_rest() is an auxiliary function of make_graph.
        """
        if not data:
            # The list "data" is empty, we are out of data
            return
        # "data" has at least one element
        currentLine = data[0].split()
        current = nodes[currentLine[0]]
        links: list[node.Node | node.DSNode | None] = []
        for i, v in enumerate(currentLine[2:]):
            if v in nodes:
                links.append(nodes[v])
                if i == 2:
                    # have to set prog.
                    # the prog this node was made with is actually a child
                    nodes[v].set_progenitor(current)  # type: ignore[union-attr]
            elif v == "0":
                links.append(None)
                nodes[v] = None
            elif int(v) % 2 == 0:  # New DSNode
                nodes[v] = node.DSNode(current)
                links.append(nodes[v])
            else:  # new SSNode
                nodes[v] = node.SSNode(current)
                links.append(nodes[v])
        assert current is not None
        current.set_links(links)
        self.build_the_rest(data[1:], nodes)

    def build_sensor(
        self, core: int, version: int, baseSeq: str
    ) -> sensor.Sensor | None:
        """
        build_sensor() first builds a 'Sensor' sequence using the 'Seed' of 'self'
        and an object of the 'random.Random' class. The sensor sequence is
        constructed. In some cases the seed graph may require the use of more bases
        than permitted by the user. In that case, this function returns None. Otherwise
        the sensors sequence is fed to RNAstructure's FoldSingleStrand. A 'Sensor' obj
        is then constructed from the resulting data. This object is returned.

        Parameters:
            core     <-- An integer, the number of the processor using this function.
            version  <-- An integer, the number of the sensor being built on this seed.
                         (ie. if this is the 1047th sensor built by proc 3 from seed
                          graph number 4: core = 3 and version = 1047).

        Returns:
            a Sensor object
        """
        self.generate_node_sizes()
        self.populate_nodes()
        seq = "".join(self.get_sequence())
        seq = seq.upper()
        # some graphs may result in sequences of larger length
        # than maxSensorSize set by user
        if len(seq) > self.maxSensorSize:
            return None

        (leadingRecDat, laggingRecDat) = self.nodes[
            self.recNodeName
        ].get_rec_seq_data()  # type: ignore[union-attr]

        # Create an RNSAstructure object
        RNA_obj = structure.RNAfolder(seq)
        sen_in = seq.lower(), RNA_obj.structure_dict
        sen = sensor.Sensor(
            sen_in, leadingRecDat, laggingRecDat, self.bindingState, self.name, baseSeq
        )

        return sen

    def generate_node_sizes(self) -> None:
        """
        generate_node_sizes() semi-randomly determines the size of the sensor, based on
        this number, a size for each node which represets physical DNA is assigned. Each
        node is given a minimum length of three bases or three base-pairs, depending on
        node type. One of the nodes contains the recognition sequence, its minimum size
        is the length of that sequence. Based on these minimums and the sensor size, the
        number of 'used bases' is calced and subtracted from the sensor size. This new
        number is the number of bases left which are then assigned, randomly, to nodes.
        This method of determinig node size, while slightly complex, avoids many issues
        of other methods which comprimize the impartiality of random node size selection
        because the sensor has a size limit.

        Parameters:
            rand   <-- a pointer to an object of the class 'random'

        Returns:
            Nothing
        """
        self.nodes[self.recNodeName].set_length(  # type: ignore[union-attr]
            len(self.recSeq)
        )  # min len of node with recSeq
        # print "rec seq node len is " + str(self.nodes[self.recNodeName].get_length())
        MAX_SIZE = self.maxSensorSize
        MIN_SIZE = 20
        size = random.randint(MIN_SIZE, MAX_SIZE)

        MIN_NODE_SIZE = 3  # to allow for loop SSNodes?

        realNodes: dict[str, tuple[node.Node, int]] = {}
        for n in self.nodes:
            # initializing "real" (ie rep. physical DNA) nodes to min size
            current = self.nodes[n]
            if current is None:  # this is not a 'real' node
                continue
            length = current.get_length()

            if length == 0:  # this is not a 'real' node
                continue

            if length == -1:  # is empty
                if current.get_state() == fold.Fold.SEQ_STATE["DS"]:  # is DS
                    realNodes[n] = (current, MIN_NODE_SIZE)
                    size -= MIN_NODE_SIZE * 2  # DS node uses 2X the number of bps
                else:  # is SS
                    realNodes[n] = (current, MIN_NODE_SIZE)
                    size -= MIN_NODE_SIZE
            else:  # is not empty (ie. has recognition seq.)
                if current.get_state() == fold.Fold.SEQ_STATE["DS"]:  # is DS
                    realNodes[n] = (current, length)
                    size -= 2 * length  # DS node uses 2X the number of bps
                else:  # is SS
                    realNodes[n] = (current, length)
                    size -= length

        keys = [n for n in realNodes]  # a list of the 'key' names in the realNodes dict
        while size > 0:
            # increasing the size of random nodes until size limit is reached
            key = random.choice(keys)
            (current, length) = realNodes[key]  # type: ignore[assignment]
            assert current is not None
            if current.get_state() == fold.Fold.SEQ_STATE["DS"]:  # is DS
                realNodes[key] = (current, length + 1)
                size -= 2  # DS node uses 2X the number of bps
            else:  # is SS
                realNodes[key] = (current, length + 1)
                size -= 1

        for r in realNodes:  # assigning the new sizes to the respective nodes
            (node, s) = realNodes[r]
            node.set_length(s)

    def populate_nodes(self) -> None:
        """
        populate_nodes() populates the empty nodes with DNA bases (ie. A, C, T, or G)
        this method requires that all nodes, which are not None, have a length.

        Parameters:
            rand    <-- a pointer to an object of a class from the module 'random'
        Returns:
            Nothing
        """
        for n in self.nodes.values():
            if n is None:
                continue
            length = n.get_length()
            seq = []
            # if this is the node with the recognition sequence we treat it differently
            if n == self.nodes[self.recNodeName]:  # is node with recognition sequence
                extra = n.get_length() - len(self.recSeq)
                # the length not required for the recSeq
                if extra != 0:
                    relLocRecSeq = random.randint(
                        1, extra
                    )  # the position of the recSeq in the node
                    n.set_relLocRecStart(relLocRecSeq)
                    # print "TADA: rel loc of rec seq is " + str(relLocRecSeq)
                    n.set_relLocRecEnd(extra - (relLocRecSeq - 1) + 1)
                    # print 'there are' + str(extra - (relLocRecSeq -1)) + 'spots left.'
                    # print 'rec seq size is ' + str(len(self.recSeq))
                    # print 'node size is ' + str(n.get_length())
                    # print 'node is ' + str(n)
                    seq = self.generate_rand_DNA_string(relLocRecSeq - 1)
                    end = self.generate_rand_DNA_string(extra - (relLocRecSeq - 1))
                    seq.extend(self.recSeq)
                    seq.extend(end)
                else:
                    # print "TADA: relitive loc of rec seq is 1."
                    n.set_relLocRecStart(1)
                    n.set_relLocRecEnd(1)
                    seq = list(self.recSeq)
            else:  # this node does not contain the recognition sequence
                seq = self.generate_rand_DNA_string(length)
            n.set_seq(seq)  # type: ignore[arg-type]

    def generate_rand_DNA_string(self, size: int) -> list[str]:
        """
        generate_rand_DNA_string() generates a list of pseudo-randomly selected
        DNA bases (ie. A, C, T, or G) of a specified size.

        Parameters:
            size   <-- an integer, the size of the desired list
            rand   <-- a pointer to an object of the 'random' class

        Returns:
            a list of random DNA letters
        """
        if size == 0:
            return []
        return [random.choice(["A", "T", "C", "G"]) for i in range(0, size)]

    def get_sequence(self) -> str:
        """
        get_sequence() returns the sequence represented by the populated nodes of
        the seed graph up to a given node. If the nodes are not populated, the function
        will return an empty sequence.

        Parameters:
            None
        Returns:
            A string. The sequence that was constructed using this seed graph.

        """
        prev = self.head
        (sequence, current) = prev.get_seq_and_next_node(None, 0)  # type: ignore

        while current is not None:
            (seq, nxt) = current.get_seq_and_next_node(
                prev, len(sequence)
            )  # type: ignore
            sequence.extend(seq)  # type: ignore[union-attr]
            prev = current  # type: ignore
            current = nxt

        return sequence
