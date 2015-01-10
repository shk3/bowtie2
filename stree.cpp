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

#include "mask.h"
#include "stree.h"

typedef SuffixTree::TOff TOff;
typedef SuffixTree::TNodeId TNodeId;
const TNodeId SuffixTree::INVALID_NODE = std::numeric_limits<TNodeId>::max();

/**
 * Check for internal consistency.
 */
bool SuffixTree::repOk() const {
	TOff nleaves = 0;
	for(size_t i = 0; i < nodes_.size(); i++) {
		if(nodes_[i].isLeaf()) {
			nleaves++;
		} else if(nodes_[i].isInternal()) {
			assert_gt(nodes_[i].numChildren(), 0);
		}
		assert_neq(1, nodes_[i].numChildren());
	}
	assert_eq(nleaves, tlen_+1);
	// Make sure we can walk to each of the leaves using fromRoot.  Make
	// sure there is a one to one map between suffixes and leaves.
	EList<bool> seen;
	seen.resize(tlen_+1);
	seen.fill(false);
	for(TOff i = 0; i < tlen_+1; i++) {
		TNodeId nodeid = rootid_;
		TChar c = -1;
		TOff depth = 0;
		if(i < tlen_) {
			ASSERT_ONLY(bool ret = ) fromRoot(i, tlen_ - i, nodeid, c, depth);
			assert(ret);
		}
		// This must be an internal node with a child
		if(depth == 0) {
			assert(nodes_[nodeid].isInternal() || nodes_[nodeid].isRoot());
			assert_neq(std::numeric_limits<TNodeId>::max(), nodes_[nodeid].out_[0]);
			nodeid = nodes_[nodeid].out_[0];
			assert_neq(nodeid, std::numeric_limits<TNodeId>::max());
			assert(nodes_[nodeid].isLeaf());
			assert(!seen[nodes_[nodeid].id_]);
			seen[nodes_[nodeid].id_] = true;
		} else {
			TNodeId childid = nodes_[nodeid].out_[c];
			assert_eq(depth+1, nodes_[childid].len_);
			assert(nodes_[childid].isLeaf());
			assert(!seen[nodes_[childid].id_]);
			seen[nodes_[childid].id_] = true;
		}
	}
	return true;
}

/**
 * Ukkonen's linear-time-and-space algorithm for building a suffix tree.
 */
void SuffixTree::ukkonen() {
	// Some extension positions are "rooted", i.e. all subsequent
	// extensions of that suffix will fall into "Case 1" whereby the root
	// edge length is incremented by 1
	TOff skipext = 1;
	TNodeId lastleaf = bottomid_;
	TOff tlen = tlen_ + 1;
	for(TOff i = 0; i < tlen - 1; i++) { // for each phase
		TNodeId lastInternal = std::numeric_limits<TNodeId>::max();
		TChar c;
		if(i+1 == tlen_) {
			c = 0;
		} else {
			c = t_[i+1] + 1; // character we're adding this phase
		}
		// lastn, lastdepth, lastoff are the node, depth and off we
		// skipCount from
		TNodeId lastn = nodes_[lastleaf].parent_;
		TOff lastdepth = i + 1 - nodes_[lastleaf].off_;
		TOff lastoff = nodes_[lastleaf].off_;
		for(TOff j = skipext; j < i+2; j++) { // for each extension
			// do skip-count
			TNodeId nodeid = std::numeric_limits<TNodeId>::max();
			TChar outc = -1;
			TOff depth = 0;
			assert_lt(lastn, tlen * 2);
			assert_leq(lastdepth, tlen);
			assert_lt(lastoff, tlen);
			skipCount(lastn, lastdepth, lastoff, nodeid, outc, depth);
			assert(nodeid != std::numeric_limits<TNodeId>::max());
			if(depth == 0) { // ended in a node
				if(lastInternal != std::numeric_limits<TNodeId>::max()) {
					// let this be the destination of the suffix link
					// outgoing from the most recently added internal node
					assert(!nodes_[lastInternal].hasSuffixLink());
					nodes_[lastInternal].slink_ = nodeid;
					lastInternal = std::numeric_limits<TNodeId>::max();
				}
				if(nodes_[nodeid].out_[c] == std::numeric_limits<TNodeId>::max()) {
					// no outgoing edge on c currently, so create one
					lastn = nodeid;
					lastdepth = 0;
					lastoff = i+1;
					lastleaf = newNode(i+1, tlen-i-1, j);
					nodes_[nodeid].out_[c] = lastleaf;
					nodes_[lastleaf].parent_ = nodeid;
				} else {
					break; // Rule 3: Already there
				}
			} else {
				// skip-count ended in edge
				TNodeId child = nodes_[nodeid].out_[outc];
				if(t_[nodes_[child].off_ + depth] + 1 != c) {
					// Rule 2: Not already there; split edge
					TNodeId midid = lastn = splitEdge(nodeid, outc, depth);
					lastdepth = 0;
					lastoff = i+1;
					// mid is a new internal node, in need of a suffix link
					lastleaf = newNode(i+1, tlen-i-1, j);
					nodes_[midid].out_[c] = lastleaf;
					nodes_[lastleaf].parent_ = midid;
					if(lastInternal != std::numeric_limits<TNodeId>::max()) {
						assert(!nodes_[lastInternal].hasSuffixLink());
						nodes_[lastInternal].slink_ = midid;
					}
					lastInternal = midid;
				} else {
					assert(lastInternal == std::numeric_limits<TNodeId>::max());
					break; // Rule 3: Already there
				}
			}
			skipext = max(skipext, j+1);
		}
		assert(lastInternal == std::numeric_limits<TNodeId>::max());
	}
}

/**
 * Find maximal exact matches (MEMs) between the indexed string t and the
 * given pattern p.  Only report MEMs that are at least l characters long.
 */
void SuffixTree::mems(
	const char *p,                         // query
	TOff plen,                             // query length
	EList<Triple<TOff, TOff, TOff> >& mem, // len, t-off, q-off
	TOff l)                                // min MEM len
{
	TNodeId node = rootid_, child = INVALID_NODE;
	TOff below = 0, depth = 0;
	suftmp_.clear();
	for(TOff i = 0; i < plen; i++) {
		int cp = firsts5[(int)p[i]] + 1;
		assert_range(1, 5, cp);
		// How do we know when we've reached the end of a MEM?  Basically,
		// whenever we set of suffixes below the current pointer changes.
		bool first = true;
		assert_leq(depth, i);
		while(true) {
			if(below == 0) {
				assert(child == INVALID_NODE);
				if(nodes_[node].out_[cp] != INVALID_NODE) {
					if(first && depth >= l) {
						// Gather suffixes below the other children,
						// besides the one we're investigating
						for(int ic = 0; ic < 6; ic++) {
							if(ic != cp && nodes_[node].out_[ic] != INVALID_NODE) {
								Locus ln(nodes_[node].out_[ic], INVALID_NODE, -1, 0);
								suffixesBelow(ln, suftmp_);
								for(size_t j = 0; j < suftmp_.size(); j++) {
									mem.expand();
									mem.back().a = depth;
									mem.back().b = suftmp_[j];
									mem.back().c = i - depth;
								}
								suftmp_.clear();
							}
						}
					}
					below++; depth++;
					child = nodes_[node].out_[cp];
					if(below == nodes_[child].len_) {
						node = child;
						child = INVALID_NODE;
						below = 0;
					}
					break; // from while(true)
				} else if(nodes_[node].isRoot()) {
					break;
				} else {
					// No outgoing edge on next pattern character
					if(first && depth >= l) {
						// Gather suffixes below
						Locus ln(node, INVALID_NODE, -1, 0);
						suffixesBelow(ln, suftmp_);
						for(size_t j = 0; j < suftmp_.size(); j++) {
							mem.expand();
							mem.back().a = depth;
							mem.back().b = suftmp_[j];
							mem.back().c = i - depth;
						}
						suftmp_.clear();
					}
					int c = -1;
					TNodeId newnode = 0;
					TOff newbelow = 0;
					skipCount(node, below, 0, newnode, c, newbelow);
					node = newnode;
					below = newbelow;
					depth--;
					if(below > 0) {
						child = nodes_[node].out_[c];
						assert_lt(below, nodes_[child].len_);
					} else {
						child = INVALID_NODE;
					}
				}
			} else {
				assert(child != INVALID_NODE);
				assert_lt(below, nodes_[child].len_);
				assert_gt(below, 0);
				TOff off = nodes_[child].off_ + below;
				assert_leq(off, tlen_);
				if(off == tlen_ || t_[off]+1 != cp) {
					if(first && depth >= l) {
						// Gather suffixes below
						Locus ln(child, INVALID_NODE, -1, 0);
						suffixesBelow(ln, suftmp_);
						for(size_t j = 0; j < suftmp_.size(); j++) {
							mem.expand();
							mem.back().a = depth;
							mem.back().b = suftmp_[j];
							mem.back().c = i - depth;
						}
						suftmp_.clear();
					}
					assert(nodes_[node].slink_ != INVALID_NODE || nodes_[node].isRoot());
					int c = -1;
					TNodeId newnode = 0;
					TOff newbelow = 0;
					skipCount(node, below, nodes_[child].off_, newnode, c, newbelow);
					node = newnode;
					below = newbelow;
					depth--;
					if(below > 0) {
						assert_range(1, 5, c);
						child = nodes_[node].out_[c];
						assert_lt(below, nodes_[child].len_);
					} else {
						child = INVALID_NODE;
					}
				} else {
					below++; depth++;
					if(below == nodes_[child].len_) {
						node = child;
						child = INVALID_NODE;
						below = 0;
					}
					break; // from while(true)
				}
			}
			first = false;
		}
	}
	if(depth >= l) {
		// Gather suffixes below
		if(below > 0) {
			assert(child != INVALID_NODE);
			node = child;
		}
		// Gather suffixes below
		Locus ln(node, INVALID_NODE, -1, 0);
		suffixesBelow(ln, suftmp_);
		for(size_t j = 0; j < suftmp_.size(); j++) {
			mem.expand();
			mem.back().a = depth;
			mem.back().b = suftmp_[j];
			mem.back().c = plen - depth;
		}
		suftmp_.clear();
	}
}

#ifdef MAIN_STREE

int main(int argc, const char **argv) {
	{
		BTDnaString s("ATAACA", true);
		SuffixTree stree;
		stree.init(s);
		assert(stree.repOk());
	}

	{
		BTDnaString s("ggtatctatatgcgctagc", true);
		SuffixTree stree;
		stree.init(s);
		assert(stree.repOk());
	}
	
	{
		SuffixTree stree;
		BTDnaString t("agaaga", true);
		BTDnaString p("aga", true);
		stree.init(t);
		EList<Triple<TOff, TOff, TOff>> mem;
		stree.mems(p, mem, 4);
		assert(mem.empty());
	}
	
	{
		SuffixTree stree;
		BTDnaString t("agaaga", true);
		BTDnaString p("aga", true);
		stree.init(t);
		EList<Triple<TOff, TOff, TOff>> mem;
		stree.mems(p, mem, 3);
		assert(!mem.empty());
		mem.sort();
		assert_eq(2, mem.size());
		assert_eq(3, mem[0].a);
		assert_eq(0, mem[0].b);
		assert_eq(0, mem[0].c);
		assert_eq(3, mem[1].a);
		assert_eq(3, mem[1].b);
		assert_eq(0, mem[1].c);
	}
	
	{
		SuffixTree stree;
		BTDnaString t("agaaga", true);
		BTDnaString p("agagaga", true);
		stree.init(t);
		EList<Triple<TOff, TOff, TOff>> mem;
		stree.mems(p, mem, 4);
		assert(mem.empty());
	}
	
	{
		SuffixTree stree;
		BTDnaString t("agaaga", true);
		BTDnaString p("agagaga", true);
		stree.init(t);
		EList<Triple<TOff, TOff, TOff>> mem;
		stree.mems(p, mem, 3);
		assert(!mem.empty());
		mem.sort();
		
		assert_eq(3, mem[0].a);
		assert_eq(0, mem[0].b);
		assert_eq(0, mem[0].c);
		
		assert_eq(3, mem[1].a);
		assert_eq(0, mem[1].b);
		assert_eq(2, mem[1].c);
		
		assert_eq(3, mem[2].a);
		assert_eq(0, mem[2].b);
		assert_eq(4, mem[2].c);
		
		assert_eq(3, mem[3].a);
		assert_eq(3, mem[3].b);
		assert_eq(0, mem[3].c);
		
		assert_eq(3, mem[4].a);
		assert_eq(3, mem[4].b);
		assert_eq(2, mem[4].c);
		
		assert_eq(3, mem[5].a);
		assert_eq(3, mem[5].b);
		assert_eq(4, mem[5].c);
	}

	{
		SuffixTree stree;
		BTDnaString t("ACG", true);
		BTDnaString p("ACG", true);
		stree.init(t);
		EList<Triple<TOff, TOff, TOff>> mem;
		stree.mems(p, mem, 3);
		assert(!mem.empty());
		mem.sort();
		assert_eq(1, mem.size());
	}

	{
		SuffixTree stree;
		BTDnaString t("ACG", true);
		BTDnaString p("ACGACG", true);
		stree.init(t);
		EList<Triple<TOff, TOff, TOff>> mem;
		stree.mems(p, mem, 3);
		assert(!mem.empty());
		mem.sort();
		assert_eq(2, mem.size());
	}

	{
		SuffixTree stree;
		BTDnaString t("A", true);
		BTDnaString p("AC", true);
		stree.init(t);
		EList<Triple<TOff, TOff, TOff>> mem;
		stree.mems(p, mem, 3);
		assert(mem.empty());
	}
}

#endif
