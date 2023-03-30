from __future__ import annotations

from typing import ClassVar


class Node:

    """
    Node is the class from which DSNode and SSNode inherit. It is a convenient place for
    getters and setters that are identical in both. It also contains a list of empty
    methods that are implemented in both classes.
    """

    def set_length(self, length: int) -> None:
        self.length = length

    def set_seq(self, seq: str) -> None:
        self.seq = seq
        self.length = len(seq)

    def set_relLocRecStart(self, loc: int) -> None:
        self.relLocRecStart = loc

    def set_relLocRecEnd(self, loc: int) -> None:
        self.relLocRecEnd = loc

    def get_length(self) -> int:
        return self.length

    def get_response(self, seq: str) -> list[str]:
        resp = {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "a": "t",
            "c": "g",
            "t": "a",
            "g": "c",
        }
        r = [resp[c] for c in seq]
        return r[::-1]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}: length={self.length},sequence={self.seq}"

    # The following methods are all empty. These methods are implemented in SSNode and
    # DSNode, so these empty methods are over-written. This is for reference, see the
    # actual implementation of these methods for their comments.

    def set_links(self, links: list[DSNode | Node | None]) -> None:
        ...

    def set_progenitor(self, prog: DSNode | SSNode | None) -> None:
        """To be implemented in DSNode and SSNode"""

    def get_state(self) -> int:  # type: ignore[empty-body]
        ...

    def get_rec_seq_data(self) -> list[dict[str, int]]:  # type: ignore[empty-body]
        ...

    def get_seq_and_next_node(
        self, prog: Node, prev_length: int
    ) -> tuple[list[str] | str, SSNode | Node | None] | None:
        ...

    def get_links(self) -> list[Node | None]:  # type: ignore[empty-body]
        ...

    def get_distance(  # type: ignore[empty-body]
        self, link1: SSNode | Node | None, link2: SSNode | Node | None
    ) -> int:
        ...

    def get_index_distance(  # type: ignore[empty-body]
        self, index1: int, index2: int
    ) -> int:
        ...

    def get_location_of(self, index: int) -> int:  # type: ignore[empty-body]
        ...

    def contains(self, index: int) -> bool:  # type: ignore[empty-body]
        ...

    def get_index_to_link_dist(  # type: ignore[empty-body]
        self, index: int, link: Node | DSNode | None, num: int
    ) -> int:
        ...


class DSNode(Node):

    """
    ----------------------------------------------------------------------
    DSNode holds all the information for a double stranded node in a
    graph which represents the folded structure of a strand of DNA.
    (The code in this class is self evident, and thus uncommented.)
    It stores: a pointer to it's upstream SSNode
               a pointer to the SSNode attatched from strand 1
               a pointer to the SSNode attatched from strand 2
               a pointer to the downstream SSNode
               the size of this node (ie. the number of bps in it)
               where strand 1 starts in the overall sequence
               where strand 2 starts in the overall sequence
    ---------------------------------------------------------------------
    """

    DIST_MULTIPLIER: ClassVar[int] = 2
    # The constant multiplier for the distance beteween bps on DS nodes

    def __init__(self, progenitor: Node | None, length: int = -1):
        """Initialize new DSNode."""
        self.upstreamSSNode = progenitor
        self.midSSNode1: SSNode | Node | None = None
        self.midSSNode2: SSNode | Node | None = None
        self.downstreamSSNode: SSNode | Node | None = None
        self.strand1Start = -1
        self.strand2Start = -1
        self.length = length
        self.seq = ""
        # set only if this node contains recSeq, is rel to begining of node
        self.relLocRecStart = -1
        self.relLocRecEnd = -1
        # set only if this node contains recSeq, is rel to end of node
        self.recSeqStart = -1  # set only if this node contains recSeq, is abs loc
        self.recRespStart = -1  # set only if this node contains recSeq, is abs loc

    def set_midSSNode1(self, ssnode: SSNode | DSNode | None) -> None:
        self.midSSNode1 = ssnode

    def set_midSSNode2(self, ssnode: SSNode | DSNode | None) -> None:
        self.midSSNode2 = ssnode

    def set_downstreamSSNode(self, ssnode: SSNode | None) -> None:
        self.downstreamSSNode = ssnode

    def set_strand1Start(self, start: int) -> None:
        self.strand1Start = start

    def set_strand2Start(self, start: int) -> None:
        self.strand2Start = start

    def set_links(self, links: list[DSNode | Node | None]) -> None:
        if len(links) != 3:
            print("Error in DSNode set links, wrong number of links given")
        self.midSSNode1 = links[0]
        self.midSSNode2 = links[1]
        self.downstreamSSNode = links[2]

    def set_progenitor(self, prog: DSNode | SSNode | None) -> None:
        self.upstreamSSNode = prog

    def get_state(self) -> int:  # replace with get_type()
        return 0  # DS

    def get_rec_seq_data(self) -> list[dict[str, int]]:
        recSeqSize = self.length - (self.relLocRecStart - 1 + self.relLocRecEnd - 1)
        return [
            {"start": self.recSeqStart, "end": self.recSeqStart + recSeqSize},
            {"start": self.recRespStart, "end": self.recRespStart + recSeqSize},
        ]

    def get_seq_and_next_node(
        self, prog: Node | None, prev_length: int
    ) -> tuple[list[str] | str, SSNode | Node | None] | None:
        """
        get_seq_and_next_node() returns a tuple consisting of:
            1) The sequence which follows the given progenitor node
            2) A pointer to the node who's sequence directly follows the sequence in 1)

        Parameter:
            prog   <-- a pointer to a node object (the progenitor of the desired seq)
            prev_length <-- length of the sequence up to this point

        Returns:
            (seq, next)     <-- a tuple with the properties described above
        """
        if prog == self.upstreamSSNode:
            self.strand1Start = prev_length + 1
            if self.relLocRecStart != -1:  # this node contains the recognition sequence
                self.recSeqStart = self.strand1Start + self.relLocRecStart - 1
            return (self.seq, self.midSSNode1)
        if prog == self.midSSNode2:
            self.strand2Start = prev_length + 1
            if self.relLocRecEnd != -1:  # this node contains the recognition sequence
                self.recRespStart = self.strand2Start + self.relLocRecEnd
            return (self.get_response(self.seq), self.downstreamSSNode)
        # error, the only two upstream nodes are accounted for
        print("Error in get_seq_and_next_node() of DSNode: False progenitor given.")
        return None

    def get_links(self) -> list[SSNode | Node | None]:
        """
        get_links() returns a list of all the links associated with this node.

        Parameters:
            None
        Returns:
            [l1, l2, l3, l4]    <-- a list of links. The first is the upstreamSSNode,
                                    then the downstream SSNode, then the first midstream
                                    SSNode (which would connect to the upstream
                                    SSNode if this DSNode were unzipped), and lastly the
                                    second midstream SSNode (which would connect to the
                                    downstream SSNode were this DSNode unzipped).
        """
        return [
            self.upstreamSSNode,
            self.downstreamSSNode,
            self.midSSNode1,
            self.midSSNode2,
        ]

    def get_distance(
        self, link1: SSNode | Node | None, link2: SSNode | Node | None
    ) -> int:
        """
        get_distance() gets the distance between two links. Because this is a DSNode,
        that distance can either be the length of the node, if the links are on opposing
        sides, or it can be 2*DSNode.DIST_MULTIPLIER, if they are on the same side.

        Parameters:
            link1   <-- a link to an SSNode
            link2   <-- a link to a different SSNode
        Returns:
            An integer, the distance according to the used metric.
        """
        if (
            (link1 is self.upstreamSSNode and link2 is self.downstreamSSNode)
            or (link2 is self.upstreamSSNode and link1 is self.downstreamSSNode)
            or (link1 is self.midSSNode1 and link2 is self.midSSNode2)
            or (link2 is self.midSSNode1 and link1 is self.midSSNode2)
        ):
            # both links are on the same side of current node, return distance as 2
            # to account for the two base pairs that link up the two nodes
            return 2 * DSNode.DIST_MULTIPLIER
        # links are on opp. sides of this node, so the distance between them is the
        # length of this node
        return self.length * DSNode.DIST_MULTIPLIER

    def get_index_distance(self, index1: int, index2: int) -> int:
        """
        get_index_distance() gets the dist between the indices of two bps on this node.
        Since this is a double stranded node, the index distance is index2-index1.

        Parameters:
            index1 <- an integer, the smaller index
            index2 <- an integer, the larger index
        """
        if not (self.contains(index1) and self.contains(index2)):
            print("Error in get_index_distance, DSNode")
            return -1

        return (
            abs(self.get_location_of(index2) - self.get_location_of(index1))
            * DSNode.DIST_MULTIPLIER
        )

    def get_location_of(self, index: int) -> int:
        """
        get_location_of() obtains the location of a bp with a particular index within
        the node. This would be the distance to the strand1Start
        or from the strand 2 end position, depending upon which
        strand the bp in question resides.

        Parameters:
            index   <- an integer, the index of the bp in question

        Returns:
            an integer, the location of the bp
        """
        if index > (self.strand1Start + self.length):
            return self.length - (index - self.strand2Start)
        return index - self.strand1Start + 1

    def contains(self, index: int) -> bool:
        """
        contains() determines if the given indexed bp is contained in this
        node.

        Parameters:
            index    <- an integer, the index of the bp in question

        Returns:
            a boolean
        """

        return (
            index < (self.strand1Start + self.length) and index >= self.strand1Start
        ) or (index < (self.strand2Start + self.length) and index >= self.strand2Start)

    def get_index_to_link_dist(
        self, index: int, link: Node | DSNode | None, num: int
    ) -> int:
        """
        get_index_to_link_dist() calculates the distance between an index in a node and
        a link at an end of the node. This method is only used by the distance funcs,
        and it has a slightly odd feature: In order to calculate the distance between
        index1 and index2, we must calc the distances from index1 to a link that will
        lead to index2, and from a link, ariving from index1, to index2. In the 1st case
        we do not wish to count index1 as part of the distance. In the second case we do
        count index2. The reason for this disparity can be seen if one tries to calc
        the distance between 1 and 10 as follows:
            1-2-3-4-5-6-Link-7-8-9-10
        We need to calculate on either side of the link seperatly. We can count 5 ints
        between 1 and the Link, and three between 10 and the link. This gives us a value
        of 8. But the dist between 1 and 10 (ie. 10-1) is obviously 9. This disparity
        is due to the fact that, in calcuating the euclidian distance between 1 and 10,
        the number 10 is included. For this reason, an extra parameter (num) has been
        added to this function. Num should be 0 if the function is calculating a dist
        between an idx and a link (ie. between 1 and Link in our example), and it should
        be 1 if calculating between a link and an index (between Link and 10 in our
        example). This is used to make the needed correction.

        Parameters:
            index   <-- an integer, the index in question
            link    <-- a pointer to an SSNode, the link in question
            num     <-- 1 or 0, as discussed above

        Returns:
            an integer, the distance.

        """
        if not self.contains(index):  # bad index
            print("Bad Index: DSNode get_index_to_link_dist from node.py")
            return -1
        loc = self.get_location_of(index)
        if link is self.upstreamSSNode or link is self.downstreamSSNode:
            return (loc - 1 + num) * DSNode.DIST_MULTIPLIER
        if link is self.midSSNode1 or link is self.midSSNode2:
            return (self.length - loc + num) * DSNode.DIST_MULTIPLIER
        # bad link
        print("Bad Link: DSNode get_index_to_link_dist node.py")
        return -2


class SSNode(Node):
    """
    -------------------------------------------------------------------
    SSNode holds all the information for a single stranded node in a
    graph which represents the folded structure of a strand of DNA.
    (The code in this class is self evident, and thus uncommented.)
    It stores: a pointer to it's upstream DSNode
               a pointer to the downstream DSNode
               the size of this node (ie. the number of bps in it)
               where this strand starts in the overall sequence
    --------------------------------------------------------------------
    """

    def __init__(self, progenitor: DSNode | Node | None, length: int = -1):
        self.upstreamDSNode = progenitor
        self.downstreamDSNode: Node | None = None
        self.start = -1
        self.length = length
        self.seq = ""
        # set only if this node contains recSeq, is rel to start of node
        self.relLocRecStart = -1
        self.relLocRecEnd = -1
        # set only if this node contains recSeq, is rel to end of node
        self.recSeqStart = -1  # set only if this node contains recSeq

    def set_downstreamDSNode(self, dsnode: DSNode | SSNode | None) -> None:
        self.downstreamDSNode = dsnode

    def set_links(self, links: list[DSNode | Node | None]) -> None:
        if len(links) != 1:
            print("Error in SSNode set_links, wrong number of links given.")
        self.downstreamDSNode = links[0]

    def set_start(self, start: int) -> None:
        self.start = start

    def set_progenitor(self, prog: DSNode | Node | None) -> None:
        self.upstreamDSNode = prog

    def get_state(self) -> int:
        return 1  # SS

    def get_rec_seq_data(self) -> list[dict[str, int]]:
        recSeqSize = self.length - (self.relLocRecStart - 1 + self.relLocRecEnd - 1)
        return [
            {"start": self.recSeqStart, "end": self.recSeqStart + recSeqSize},
            {"start": -1, "end": -1},
        ]  # response seq DNE for SSNodes

    def get_seq_and_next_node(
        self, prog: Node, prev_length: int
    ) -> tuple[str, Node | None] | None:
        """
        get_seq_and_next_node() returns a tuple consisting of:
            1) The sequence which follows the given progenitor node
            2) A pointer to the node who's sequence directly follows the sequence in 1)

        Parameter:
            prog   <-- a pointer to a node object (the progenitor of the desired seq)
            prev_length <-- length of seq to this point
        Returns:
            (seq, next)     <-- a tuple with the properties described above
        """
        if prog == self.upstreamDSNode:
            self.start = prev_length + 1
            if self.relLocRecStart != -1:  # this node contains the recognition sequence
                self.recSeqStart = self.start + self.relLocRecStart - 1
            return (self.seq, self.downstreamDSNode)
        # error, the only upstream node is accounted for
        print("Error in get_seq_and_next_node() of SSNode: False progenitor given.")
        return None

    def get_links(self) -> list[Node | DSNode | None]:
        """
        get_links() returns a list of all the links associatd with this node.

        Parameters:
            None
        Returns:
            [l1,l2] <-- a list wherein the first link is to the upstream DSNode
                        and the second link is to the downstream DSNode
        """
        return [self.upstreamDSNode, self.downstreamDSNode]

    def get_distance(self, link1: Node | None, link2: Node | None) -> int:
        """
        get_distance() gets the distance between 2 links of this node.

        Parameters:
            link1   <-- a DSNode link
            link2   <-- a different DSNode link
        """
        # this is an SSNode so distance beween links is always the length of
        # this node
        return self.length

    def get_index_distance(self, index1: int, index2: int) -> int:
        """
        get_index_distance() gets the distance between two bps on this node based upon
        their respective indices. Since this is an SSNode, that is calculated as
        (index2-index1)/2

        Parameters:
            index1 <-- the smaller index
            index2 <-- the larger index
        """
        if not (self.contains(index1) and self.contains(index2)):
            print("Error from get_index_distance() SSNode")
            return -1
        return self.get_location_of(index2) - self.get_location_of(index1)

    def get_location_of(self, index: int) -> int:
        """
        get the location of a bp with a particular index within
        the node. This would be the distance to the start of the
        strand.

        Parameters:
            index   <- an integer, the index of the bp in question

        Returns:
            an integer, the location of the bp
        """
        return index - self.start

    def contains(self, index: int) -> bool:
        """
        contains determines if the given indexed bp is contained in this
        node.

        Parameters:
            index   <- an integer, the index of the bp in question

        Returns:
            a boolean, yes if the index is contained in this node.
        """
        if index < (self.length + self.start) and index >= self.start:
            return True
        return False

    def get_index_to_link_dist(
        self, index: int, link: Node | DSNode | None, num: int
    ) -> int:
        """
        get_index_to_link_dist() calculates the distance between an index in a node and
        a link at an end of the node. This method is only used by the dist functions,
        and it has a slightly odd feature: In order to calculate the distance between
        index1 and index2, we must calculate the dists from index1 to a link that will
        lead to index2, and from a link, ariving from index1, to index2. In the 1st case
        we do not wish to count index1 as part of the distance. In the second case we do
        count index2. The reason for this disparity can be seen if one tries to calc
        the distance between 1 and 10 as follows:
            1-2-3-4-5-6-Link-7-8-9-10
        We need to calculate on either side of the link seperatly. We can count 5 ints
        between 1 and the Link, and three between 10 and the link. This gives us a value
        of 8. But the dist between 1 and 10 (ie. 10-1) is obviously 9. This disparity
        is due to the fact that, in calcuating the euclidian distance between 1 and 10,
        the number 10 is included. For this reason, an extra parameter (num) has been
        added to this function. Num should be 0 if the function is calculating a dist
        between an idx and a link (ie. between 1 and Link in our example), and it should
        be 1 if calculating between a link and an index (between Link and 10 in our
        example). This is used to make the needed correction.

        Parameters:
            index   <-- an integer, the index in question
            link    <-- a pointer to a DSNode, the link in question
            num     <-- 1 or 0, as discussed above

        Returns:
            an integer, the distance.

        """
        if not self.contains(index):  # bad index
            print("Bad Index: SSNode get_index_to_link_dist from node.py")
            return -1
        distToDownstream = self.start + self.length - index - 1 + num
        distToUpstream = index - self.start + num
        if link is self.upstreamDSNode and link is self.downstreamDSNode:  # loop node
            return (
                distToUpstream
                if distToUpstream < distToDownstream
                else distToDownstream
            )
        if link is self.upstreamDSNode:
            return distToUpstream
        if link is self.downstreamDSNode:
            return distToDownstream
        # bad link
        print("Bad Link: SSNode get_index_to_link_dist from node.py")
        return -2
