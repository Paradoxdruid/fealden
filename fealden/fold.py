import math
import sys
from typing import ClassVar

from . import node


class Fold:
    """The constructor for Fold.

    Parameters:
        foldData    <- a list of fold data in the format of the second
                        value of the tuple which is output by simplify_input.
        deltaG      <- a float, the deltaG of the fold (we got this from unafold)
        recSeq      <- a string, the recognition sequence

    """

    # the states a sequence can be in
    SEQ_STATE: dict[str, int] = {"DS": 0, "SS": 1, "MIXED": 3}
    RT: ClassVar[float] = 8.3144598 * (1.0 / 4184.0) * 298.0

    def __init__(
        self, foldData: list[list[int]], deltaG: float, recSeq: dict[str, int]
    ) -> None:
        """Initialize new Fold obj."""
        self.head = node.SSNode(None)
        self.deltaG = deltaG
        self.conc = math.e ** (-self.deltaG / Fold.RT)
        self.foldData = foldData
        # for bp i, ptrList[i-1]= ptr to node that bp i belongs to
        self.ptrList: list[None | node.SSNode | node.DSNode] = [None] * len(
            self.foldData
        )
        self.construct_graph_SSNode(self.head, 0)
        self.recSeq = recSeq
        self.recSeqState: int = self.get_rec_seq_state()

    def construct_graph_SSNode(
        self, currentNode: node.SSNode, currentIndex: int
    ) -> None:
        """
        consturuct_graph_SSNode is the algorithem for constructing the
        fold graph given that the current node is an SSNode.

        Parameters:
            currentNode  <- the SSNode we're currently on.
            currentIndex <- the index in self.graphData which we're currently at.
        Returns:
            Nothing
        """
        currentNode.set_start(currentIndex + 1)
        isLastNode = True
        length = 0
        for i, v in enumerate(self.foldData[currentIndex::]):
            if v[1] == 0:
                self.ptrList[i + currentIndex] = currentNode
                length = i + 1  # keep track in case this is the last node
            else:
                # print i
                isLastNode = False
                currentNode.set_length(i)
                nextNode: node.DSNode | node.SSNode | None = None
                if self.ptrList[v[1] - 1] is None:
                    nextNode = node.DSNode(currentNode)
                    self.construct_graph_DSNode_strand1(nextNode, currentIndex + i)
                else:
                    nextNode = self.ptrList[v[1] - 1]
                    assert not isinstance(nextNode, node.SSNode)
                    self.construct_graph_DSNode_strand2(
                        nextNode, currentIndex + i, currentNode
                    )

                currentNode.set_downstreamDSNode(nextNode)
                break
        if isLastNode:
            currentNode.set_length(length)

    def construct_graph_DSNode_strand1(
        self, currentNode: node.DSNode, currentIndex: int
    ) -> None:
        """
        construct_graph_DSNode_strand1 is the algorithem for building a
        DSNode into the graph given that the base pairs in the DSNode that
        we're processing are all on the "first" strand of the DSNode. (ie
        the location number of these bps on the strand is less then that
        of their response bps).

        Parameters:
            currentNode     <- the DSNode we're currently working on
            currentIndex    <- the index in the self.foldData that we're currently on
        Returns:
            Nothing
        """
        currentNode.set_strand1Start(currentIndex + 1)
        prevPair = self.foldData[currentIndex][1] + 1
        for i, v in enumerate(self.foldData[currentIndex::]):
            if v[1] == prevPair - 1:
                self.ptrList[i + currentIndex] = currentNode
                self.ptrList[prevPair - 2] = currentNode
                prevPair -= 1
            else:
                currentNode.set_length(i)
                nextNode = node.SSNode(currentNode)
                self.construct_graph_SSNode(nextNode, currentIndex + i)
                currentNode.set_midSSNode1(nextNode)

                break

    def construct_graph_DSNode_strand2(
        self,
        currentNode: node.DSNode | None,
        currentIndex: int,
        prevNode: node.SSNode | node.DSNode | None,
    ) -> None:
        """
        construct_graph_DSNode_strand2 is the algorithem for building a
        DSNode into the graph given that the base pairs in the DSNode that
        we're processing are all on the "second" strand of the DSNode. (ie
        the location number of these bps on the strand is greater then that
        of their response bps).

        Parameters:
            currentNode     <- the DSNode we're currently working on
            currentIndex    <- the index in the self.foldData that we're currently on
            prevNode        <- the last SSNode we've worked on
        Returns:
            Nothing
        """
        assert currentNode is not None
        currentNode.set_strand2Start(currentIndex + 1)
        currentNode.set_midSSNode2(prevNode)
        prevPair = self.foldData[currentIndex][1] + 1
        for i, v in enumerate(self.foldData[currentIndex::]):
            if (v[1] == prevPair - 1) and (v[1] != 0):
                # self.ptrList[i+currentIndex] = currentNode
                # #already set when building strand 1
                # self.ptrList[prevPair-2] = currentNode #already set when
                # building strand 1
                prevPair -= 1
            else:
                nextNode = node.SSNode(currentNode)
                self.construct_graph_SSNode(nextNode, currentIndex + i)
                currentNode.set_downstreamSSNode(nextNode)
                # don't set length, it's already been set while building strand
                # 1
                break

    def get_distance(self, index1: int, index2: int) -> int:
        """
        get_distance (is a proper distance metric which) captures the approx. spacial
        distance between two base pairs

        Parameters:
            firstIndex  <- an integer, the index of the first bp (the smaller index)
            secondIndex <- an integer, the index of the second bp (the larger index)
        Returns:
            distance   <- an integer, the calculated distance
        """
        if index1 > index2:
            # TODO: Double-check this implementation
            # temp = index1
            index1 = index2
            index2 = index1

        node1 = self.ptrList[index1 - 1]
        node2 = self.ptrList[index2 - 1]

        if node1 == node2:
            assert node1 is not None
            return node1.get_index_distance(index1, index2)
        dist = sys.maxsize
        assert node1 is not None
        links = node1.get_links()
        for each in links:
            distToIndex2 = self.get_dist_to_index(index2, [node1], node1, each)

            assert node1 is not None
            tempDist = distToIndex2 + node1.get_index_to_link_dist(index1, each, 0)
            if tempDist < dist:
                dist = tempDist
        return dist

    def get_dist_to_index(
        self,
        index: int,
        traversed: list[node.Node],
        previous: node.Node,
        current: node.Node | node.DSNode | None,
    ) -> int:
        """
        get_dist_to_index() gets the distance from the start of the current
        node to the base pair at the index referenced.

        Parameters:
            index       <- the target, an integer
            traversed   <- a list of nodes we have previously traversed
                           (To avoid going in circles).
            previous    <- the node most recently traversed
            current     <- the node being traversed

        """
        if current is None:  # reached end and have not found index
            return sys.maxsize
        if current.contains(index):
            return current.get_index_to_link_dist(index, previous, 1)
        dist = sys.maxsize
        links = current.get_links()
        traversed.append(current)
        for each in links:
            if each in traversed:
                continue
            tempDist = current.get_distance(previous, each) + self.get_dist_to_index(
                index, [j for j in traversed], current, each
            )
            if tempDist < dist:
                dist = tempDist

        return dist

    def get_rec_seq_state(self) -> int:
        """
        get_rec_seq_state() gets the state (ie DS, SS, or Mixed) into which the rec
        sequence has folded.

        Parameters:
            none
        Returns:
            the state  <- an item from the Fold.SEQ_STATE list.
        """
        recSeqPtrList = self.ptrList[self.recSeq["start"] - 1 : self.recSeq["end"] - 1]
        assert recSeqPtrList != []
        startingNode = recSeqPtrList[0]

        # TODO Implement Mixed Req Seqs
        # recSeqMixed = False
        # for p in recSeqPtrList:
        #     if p != startingNode:
        #         recSeqMixed = True
        # if recSeqMixed:
        # print "Rec seq mixed"
        #    return Fold.SEQ_STATE["MIXED"]
        # print "Rec seq " + str(startingNode.get_state())
        assert startingNode is not None
        return startingNode.get_state()
