/*
 * Copyright 2013, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef STREE_H_
#define STREE_H_

#include <iostream>
#include <stdint.h>
#include "ds.h"
#include "sstring.h"

/**
 * Encapsulates a suffix tree, using Ukkonen's algorithm to build it.
 */
class SuffixTree {

public:
	
	typedef int64_t TOff;
	typedef int64_t TNodeId;
	typedef int TChar;
	
	static const TNodeId INVALID_NODE;
	
protected:
	
	/**
	 * Node of a suffix tree.  For a node, we store (1) offset and (2) length
	 * of the substring labeling the edge leading *to* the node from its
	 * parent.  We also store (3) a link to the parent (NULL in the case of the
	 * root), (4) a suffix link (NULL in the case of root or leaf), and (5)
	 * outgoing edges for each character of the DNA alphabet, including N.
	 */
	struct Node {
		
		Node() {
			off_ = len_ = std::numeric_limits<TOff>::max();
		}
		
		Node(TOff off, TOff len, TOff id) {
			init(off, len, id);
		}
		
		/**
		 * Initialize with a reference offset and length.
		 */
		void init(TOff off, TOff len, TOff id) {
			reset();
			off_ = off;
			len_ = len;
			id_ = id;
		}
		
		/**
		 * Reset node to uninitialized state, with no parent, suffix link or
		 * children.
		 */
		void reset() {
			parent_ = slink_ =
			out_[0] = out_[1] = out_[2] =
			out_[3] = out_[4] = out_[5] = INVALID_NODE;
			id_ = 0;
			// 0=$, 1=A, 2=C, 3=G, 4=T, 5=N
		}
		
		/**
		 * Return true iff this is a leaf node.
		 */
		bool isLeaf() const {
			return out_[0] == INVALID_NODE &&
			       out_[1] == INVALID_NODE &&
				   out_[2] == INVALID_NODE &&
				   out_[3] == INVALID_NODE &&
				   out_[4] == INVALID_NODE &&
				   out_[5] == INVALID_NODE;
		}
		
		/**
		 * Return the number of children the node has.  Should never be 1.
		 */
		size_t numChildren() const {
			size_t nc = 0;
			for(size_t i = 0; i <= 5; i++) {
				if(out_[i] != INVALID_NODE) {
					nc++;
				}
			}
			return nc;
		}
		
		/**
		 * Return true iff this is the root node.
		 */
		bool isRoot() const {
			return len_ == 0;
		}
		
		/**
		 * Return true iff this is an internal node.
		 */
		bool isInternal() const {
			return !isLeaf() && !isRoot();
		}
		
		/** 
		 * Return true iff this currently has a non-NULL suffix link.  During
		 * suffix tree construction, this might return false for an internal node.
		 */
		bool hasSuffixLink() const {
			return slink_ != INVALID_NODE;
		}
		
		TOff off_; // offset of substring labeling edge from parent
		TOff len_; // length of substring labeling edge from parent
		TOff id_;  // identifier; suffix offset for leaf nodes
		TNodeId parent_; // parent link (NULL for root)
		TNodeId slink_;  // suffix link (NULL for root, leaf)
		TNodeId out_[6]; // children on $, A, C, G, T, N
	};

	/**
	 * A location in the suffix tree.  Could be on a node or in the middle of
	 * an edge.
	 */
	struct Locus {

		Locus() { }
		
		Locus(
			TNodeId nodeAbove,
			TNodeId nodeBelow,
			int outChar,
			TOff charsBelow)
		{
			init(nodeAbove, nodeBelow, outChar, charsBelow);
		}
		
		/**
		 * Initialize with a reference offset and length.
		 */
		void init(
			TNodeId nodeAbove,
			TNodeId nodeBelow,
			int outChar,
			TOff charsBelow)
		{
			na_ = nodeAbove;
			nb_ = nodeBelow;
			oc_ = outChar;
			below_ = charsBelow;
		}
		
		/**
		 * Reset node to uninitialized state, with no parent, suffix link or
		 * children.
		 */
		void reset() {
			na_ = nb_ = INVALID_NODE;
			oc_ = -1;
			below_ = 0;
		}
		
		/**
		 * Check that the locus is internally consistent.
		 */
		bool repOk() const {
			return true;
		}
		
		TNodeId na_;
		TNodeId nb_;
		int oc_;
		TOff below_;
	};
	
public:
	
	SuffixTree() { reset(); }
	
	/**
	 * Initialize with string t.
	 */
	SuffixTree(const BTDnaString& t) {
		init(t.buf(), t.length());
	}
	
	/**
	 * Initialize with string t with given length.
	 */
	SuffixTree(const char *t, TOff tlen) {
		init(t, tlen);
	}

	/**
	 * Initialize with string t.
	 */
	void init(const BTDnaString& t) {
		init(t.buf(), t.length());
	}
	
	/**
	 * Initialize with string t with given length.
	 */
	void init(const char *t, TOff length) {
		reset();
		assert_gt(length, 0);
		t_ = t;
		tlen_ = length;
		rootid_ = nodes_.alloc();
		bottomid_ = nodes_.alloc();
		nodes_[rootid_].init(0, 0, 0);
		nodes_[bottomid_].init(0, tlen_ + 1, 0);
		nodes_[bottomid_].parent_ = rootid_;
		nodes_[rootid_].out_[t_[0]+1] = bottomid_;
		ukkonen();
		assert(repOk());
	}
	
	/**
	 * Reset to uninitialized state.
	 */
	void reset() {
		nodes_.clear();
		suftmp_.clear();
		t_ = NULL;
		tlen_ = 0;
		rootid_ = bottomid_ = INVALID_NODE;
	}
	
	/**
	 * Return true iff suffix tree is initialized.
	 */
	bool inited() const {
		return t_ != NULL;
	}
	
	/**
	 * Check for internal consistency.
	 */
	bool repOk() const;
	
	/**
	 * Find maximal exact matches (MEMs) between the indexed string t and the
	 * given pattern p.  Only report MEMs that are at least l characters long.
	 */
	void mems(
		const BTString& p,                     // query
		EList<Triple<TOff, TOff, TOff> >& mem, // len, t-off, q-off
		TOff l)                                // min MEM len
	{
		mems(p.buf(), p.length(), mem, l);
	}

	/**
	 * Find maximal exact matches (MEMs) between the indexed string t and the
	 * given pattern p.  Only report MEMs that are at least l characters long.
	 */
	void mems(
		const char *p,                         // query
		TOff plen,                             // query length
		EList<Triple<TOff, TOff, TOff> >& mem, // len, t-off, q-off
		TOff l);                               // min MEM len
	
protected:

	/**
	 * Create new node and return its id.
	 */
	TNodeId newNode(TOff off, TOff len, TOff id) {
		TOff allocid = nodes_.alloc();
		nodes_[allocid].init(off, len, id);
		return allocid;
	}

	/**
	 * Split an edge.
	 */
	TNodeId splitEdge(
		TNodeId nodeid, // node above edge to be split
		TChar outc,     //
		TOff depth)     // # chars below the parent node to make the split
	{
		assert_gt(depth, 0);
		TNodeId childid = nodes_[nodeid].out_[outc];
		assert_lt(depth, nodes_[childid].len_);
		TNodeId midid = newNode(nodes_[childid].off_, depth, 0);
		Node& node = nodes_[nodeid];
		Node& mid = nodes_[midid];
		Node& child = nodes_[childid];
		mid.parent_ = nodeid;
		child.parent_ = midid;
		node.out_[outc] = midid;
		child.off_ += depth;
		child.len_ -= depth;
		mid.out_[t_[child.off_]+1] = childid;
		return midid;
	}
	
	/**
	 * Walk down from the root, following a path corresponding to string q.
	 * Return coordinates of where walk ended.
	 */
	bool fromRoot(
		TOff off,            // in: offset of first character in T for walk
		TOff d,              // in: # characters in T for walk
		TNodeId& node_o,     // out: node walked to
		TChar& c_o,          // out: depth below node
		TOff& depth_o) const // out:
	{
		TOff i = 0;
		TNodeId curid = 0; // start at root
		while(i < d) {
			TChar c = t_[off + i]+1; // next character
			assert_range(0, 5, c);
			if(nodes_[curid].out_[c] == INVALID_NODE) { // child with that character?
				return false; // no such child!
			}
			TNodeId childid = nodes_[curid].out_[c];
			assert_eq(c, t_[nodes_[childid].off_]+1);
			i++;
			// match rest of characters past off + i, stopping when we exhaust
			// the query substring or when we've exhausted the characters
			// labeling the edge
			TOff j = 1;
			TOff edgeLen = nodes_[childid].len_, edgeOff = nodes_[childid].off_;
			while(j < edgeLen && i < d) {
				if(t_[off + i] != t_[edgeOff + j]) {
					return false; // failed to match character along edge!
				}
				j++; i++; // advance along query and edge
			}
			if(j == edgeLen) {
				curid = childid; // exhausted edge, descend to child
			} else {
				assert_eq(i, d); // exhausted query
				node_o = curid;
				c_o = c;
				depth_o = j;
				return true;
			}
		}
		// exhausted query and edge at the same time
		node_o = curid;
		depth_o = 0;
		return true;
	}
	
	/**
	 * Use skip-count trick to walk from an internal node labeled cx (c is a
	 * char, x is a string) to the location in the tree labeled x.
	 */
	void skipCount(
		TNodeId v,           // in:
		TOff d,              // in:
		TOff off,            // in:
		TNodeId& node_o,     // out:
		TChar& c_o,          // out:
		TOff& depth_o) const // out:
	{
		if(!nodes_[v].isRoot() && !nodes_[v].hasSuffixLink()) {
			// at an internal node without a suffix link; must walk up by one
			// more node to get to either (a) the root or (b) an internal node
			// with a suffix link
			off = nodes_[v].off_; // adjust off accordingly
			d += nodes_[v].len_;  // adjust d accordingly
			v = nodes_[v].parent_;
			assert(nodes_[v].isRoot() || nodes_[v].hasSuffixLink());
		}
		// now v is either the root or an internal node with a suffix link
		if(nodes_[v].isRoot()) { // v is the root
			assert_gt(d, 0);
			d--; off++; // adjust d and off; we didn't traverse a suffix link
			if(d > 0) {
				// walk down from the root
				ASSERT_ONLY(bool ret = ) fromRoot(off, d, node_o, c_o, depth_o);
				assert(ret);
				return;
			} else {
				node_o = v; // we stopped at the root
				depth_o = 0;
			}
		} else { // v is an internal node with a suffix link
			assert(nodes_[v].hasSuffixLink());
			TNodeId sv = nodes_[v].slink_;
			TNodeId curid = sv;
			assert(!nodes_[sv].isLeaf());
			if(d == 0) {
				// journey ends at sv
				node_o = sv;
				depth_o = 0;
				return;
			}
			// walk down from sv
			TOff i = off;
			while(d > 0) {
				int c = t_[i] + 1;
				TNodeId nextid = nodes_[curid].out_[c];
				if(nodes_[nextid].len_ < d) {
					d -= nodes_[nextid].len_;
					i += nodes_[nextid].len_;
					curid = nextid;
				} else if(nodes_[nextid].len_ == d) {
					node_o = nextid;
					depth_o = 0;
					return;
				} else {
					node_o = curid;
					c_o = c;
					depth_o = d;
					return;
				}
			}
		}
	}

	/**
	 * Descend into all the leaf nodes below the given locus and append their
	 * ids to the sufs list.
	 */
	void suffixesBelow(Locus& l, EList<TOff>& sufs) {
		assert(l.repOk());
		if(l.below_ > 0) {
			l.na_ = l.nb_;
			l.nb_ = INVALID_NODE;
			l.below_ = 0;
			l.oc_ = -1;
		}
		if(nodes_[l.na_].isLeaf()) {
			sufs.push_back(nodes_[l.na_].id_);
		} else {
			// Recursively look for suffixes below children
			for(size_t i = 0; i < 6; i++) {
				if(nodes_[l.na_].out_[i] != INVALID_NODE) {
					Locus lc(nodes_[l.na_].out_[i], INVALID_NODE, -1, 0);
					suffixesBelow(lc, sufs);
				}
			}
		}
	}

	/**
	 * Ukkonen's linear-time-and-space algorithm for building a suffix tree.
	 */
	void ukkonen();
	
	TNodeId rootid_;
	TNodeId bottomid_;
	const char* t_;
	TOff tlen_;
	EFactory<Node> nodes_;
	
	EList<TOff> suftmp_;
};

#endif /* defined(STREE_H_) */
