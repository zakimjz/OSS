//! \file lattice_node.h - struct to represent a node on itemset lattice
#ifndef _LATTICE_NODE_H
#define _LATTICE_NODE_H

#include <ext/hash_set>
#include "graph_iso_check.h"
#include <algorithm>
#include "helper_funs.h"
#include "pattern_factory.h"

//! A lattice_node template structure
template <typename PAT >
struct lattice_node
{

  typedef lattice_node<PAT> L_NODE;
  typedef typename PAT::VERTEX_T V_T;
  typedef typename PAT::EDGE_T E_T;
  typedef pair<int, int> EDGE;

  string get_key() {
    const typename PAT::CAN_CODE& cc = _pat->canonical_code();
    std::string min_dfs_cc = cc.to_string();
    return min_dfs_cc;
  }
 
	//! Constructor
  lattice_node(PAT* p) {
    _pat = p;
  }
  
  bool _is_processed; //!< it is true, when we know all neighbors and their status of this pattern
  PAT* _pat;//!< Store a pattern in lattice node
  vector<L_NODE*> _neighbors;//!< Store all the neighbors of a node
  vector<double> _neighbor_prob; 
  int _super_cnt;
};

#endif
