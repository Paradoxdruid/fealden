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
        self, fold_data: list[list[int]], deltaG: float, rec_seq: dict[str, int]
    ) -> None:
        """Initialize new Fold obj."""
        self.head = node.SSNode(None)
        self.deltaG = deltaG
        self.conc = math.e ** (-self.deltaG / Fold.RT)
        self.fold_data = fold_data
        # for bp i, ptrList[i-1]= ptr to node that bp i belongs to
        self.ptr_list: list[None | node.SSNode | node.DSNode] = [None] * len(
            self.fold_data
        )
        self.construct_graph_SSNode(self.head, 0)
        self.rec_seq = rec_seq
        self.rec_seq_state: int = self.get_rec_seq_state()

    def __repr__(self) -> str:
        return f"Fold: head={repr(self.head)}, deltaG={self.deltaG}, conc={self.conc},\
 foldData={self.fold_data}, ptrList={self.ptr_list}, recSeq={self.rec_seq},\
 recSeqState={self.rec_seq_state}"

    def construct_graph_SSNode(
        self, current_node: node.SSNode, current_index: int
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
        current_node.set_start(current_index + 1)
        is_last_node = True
        length = 0
        for i, v in enumerate(self.fold_data[current_index::]):
            if v[1] == 0:
                self.ptr_list[i + current_index] = current_node
                length = i + 1  # keep track in case this is the last node
            else:
                is_last_node = False
                current_node.set_length(i)
                next_node: node.DSNode | node.SSNode | None = None
                if self.ptr_list[v[1] - 1] is None:
                    next_node = node.DSNode(current_node)
                    self.construct_graph_DSNode_strand1(next_node, current_index + i)
                else:
                    next_node = self.ptr_list[v[1] - 1]
                    assert not isinstance(next_node, node.SSNode)
                    self.construct_graph_DSNode_strand2(
                        next_node, current_index + i, current_node
                    )

                current_node.set_downstream_DSNode(next_node)
                break
        if is_last_node:
            current_node.set_length(length)

    def construct_graph_DSNode_strand1(
        self, current_node: node.DSNode, current_index: int
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
        current_node.set_strand_1_start(current_index + 1)
        prev_pair = self.fold_data[current_index][1] + 1
        for i, v in enumerate(self.fold_data[current_index::]):
            if v[1] == prev_pair - 1:
                self.ptr_list[i + current_index] = current_node
                self.ptr_list[prev_pair - 2] = current_node
                prev_pair -= 1
            else:
                current_node.set_length(i)
                next_node = node.SSNode(current_node)
                self.construct_graph_SSNode(next_node, current_index + i)
                current_node.set_mid_SSNode1(next_node)

                break

    def construct_graph_DSNode_strand2(
        self,
        current_node: node.DSNode | None,
        current_index: int,
        prev_node: node.SSNode | node.DSNode | None,
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
        assert current_node is not None
        current_node.set_strand_2_start(current_index + 1)
        current_node.set_mid_SSNode2(prev_node)
        prev_pair = self.fold_data[current_index][1] + 1
        for i, v in enumerate(self.fold_data[current_index::]):
            if (v[1] == prev_pair - 1) and (v[1] != 0):
                # self.ptrList[i+currentIndex] = currentNode
                # #already set when building strand 1
                # self.ptrList[prevPair-2] = currentNode #already set when
                # building strand 1
                prev_pair -= 1
            else:
                next_node = node.SSNode(current_node)
                self.construct_graph_SSNode(next_node, current_index + i)
                current_node.set_downstream_SSNode(next_node)
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

        node1 = self.ptr_list[index1 - 1]
        node2 = self.ptr_list[index2 - 1]

        if node1 == node2:
            assert node1 is not None
            return node1.get_index_distance(index1, index2)
        dist = sys.maxsize
        assert node1 is not None
        links = node1.get_links()
        for each in links:
            dist_to_index_2 = self.get_dist_to_index(index2, [node1], node1, each)

            assert node1 is not None
            temp_dist = dist_to_index_2 + node1.get_index_to_link_dist(index1, each, 0)
            if temp_dist < dist:
                dist = temp_dist
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
            temp_dist = current.get_distance(previous, each) + self.get_dist_to_index(
                index, traversed, current, each
            )
            if temp_dist < dist:
                dist = temp_dist

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
        rec_seq_ptr_list = self.ptr_list[
            self.rec_seq["start"] - 1 : self.rec_seq["end"] - 1
        ]
        assert rec_seq_ptr_list != []
        starting_node = rec_seq_ptr_list[0]

        # TODO Implement Mixed Req Seqs
        # recSeqMixed = False
        # for p in recSeqPtrList:
        #     if p != startingNode:
        #         recSeqMixed = True
        # if recSeqMixed:
        #    return Fold.SEQ_STATE["MIXED"]

        assert starting_node is not None
        return starting_node.get_state()
