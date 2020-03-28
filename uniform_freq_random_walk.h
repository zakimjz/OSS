//! \file random_walk_manager_freq.h - class to perform the random walk
#ifndef _RANDOM_WALK_MANAGER2_H_
#define _RANDOM_WALK_MANAGER2_H_

#include <algorithm>
#include <ext/hash_set>
#include <ext/hash_map>
#include "helper_funs.h"
#include "pattern_factory.h"
#include "lattice_node.h"
#include "random.h"
#include "time_tracker.h"
#include "functional"

/**
 * This performs random walk on frequent patterns to find uniform sample of frequent
 * patterns. It save the frequent patterns for fast processing
 */
//! A RandomWalkManager_Freq template class.
template<class PAT>
class Uniform_Freq_Random_Walk
{

  public:

  typedef lattice_node<PAT> LATTICE_NODE;
  typedef HASHNS::hash_map<string, LATTICE_NODE*, hash_func<string>, equal_to<string> > NODE_MAP;
  typedef typename NODE_MAP::iterator NS_IT;
  typedef typename NODE_MAP::const_iterator CONST_NS_IT;

  typedef HASHNS::hash_map<string, int, hash_func<string>, equal_to<string> > FREQ_CNT_MAP;
  typedef typename FREQ_CNT_MAP::iterator FC_IT;
  typedef Database<PAT> DATABASE;
  typedef PatternFactory<PAT> PATTERN_FACTORY;

	//! Constructor
  Uniform_Freq_Random_Walk(DATABASE* d, int iter) {
    _last = 0;
    _pf = PATTERN_FACTORY::instance(d);
    _maxiter = iter;
    cout << "Iteration count:" << _maxiter << endl;
  }
	
	//! get PatternFactory object
  PatternFactory<PAT>* get_pat_factory() const {
    return _pf;
  }

  // random walk manager initialize with a frequent pattern node
	/*! \fn LATTICE_NODE* initialize() 
 		*  \brief A member function to initialize the walk in itemset Lattice. Initialization completed
		*		by selecting an size one frequent pattern.
 		*  \return a pointer of LATTICE_NODE type.
 		*/
  LATTICE_NODE* initialize() {
		
    PAT* p =_pf->get_one_random_one_edge_frequent_pattern();
    
    const typename PAT::CAN_CODE& cc = check_isomorphism(p);
    p->set_canonical_code(cc);
    cout << p->get_canonical_code().to_string() << endl;
    if (p->_is_frequent==false) {
      cout << "ERROR:this pattern is infrquent\n";
      exit(1);
    }
    LATTICE_NODE* ln = create_lattice_node(p);
    process_node(ln);
    return ln;
  }

	/*! \fn void walk(LATTICE_NODE* current, int &iter) 
 		*  \brief A member function to do the random walk on itemset Lattice. Walk continues until
		* iteration count hits maximum number of iteration. From a current node this method selects
		* next node to jump based on acceptance probability calculated by Metropolis-Hastings 
		* algorithm.
		*	\param current a pointer of LATTICE_NODE.
		* \param iter a reference of an integer 
	*/
 	bool walk(LATTICE_NODE* current, int &iter) {

    int step = 1;
    while (true) {
      process_node(current);
			if(current->_neighbors.size()==0)// if current has no neighbors return and start over again.
			{
				return 0;						
			}
      PAT* p = current->_pat;
      if (iter >= _maxiter) {
				tt.stop();
				cout<< "time taken="<<tt.print()<<endl;
        cout << "Dist:";
        FC_IT fit;
        for (fit = _freq_cnt.begin(); fit != _freq_cnt.end(); fit++) {
          cout << fit->first << "(" << fit->second << ")" << endl;
        }
        cout << "maxiter:" << _maxiter << " iter:" << iter << " Returing from here" << endl;				
        return 1;
      }

      if (step >=10) {  // assuming after 10 steps, the walk mixes to uniformity
        iter++;
        string cc = p->get_canonical_code().to_string();
        FC_IT fit = _freq_cnt.find(cc);
        if (fit != _freq_cnt.end()) {
          fit->second++;
        }
        else {
          _freq_cnt.insert(make_pair(cc, 1));
          cout << "Iter:" << iter << " Total sampled:" << _freq_cnt.size() << " Total traversed:" << _node_map.size() << endl;
        }
      }
      // cout << "In walk: looking for new node to visit" << endl;
      LATTICE_NODE* next = get_next(current);
      _last = p;
      current = next;
      step++;
    }
  }

	/*! \fn LATTICE_NODE* get_next(LATTICE_NODE* current) const 
 		*  \brief A member function to get next node on itemset lattice to jump from current. 
		*	 Acceptance probability calculation of Metropolis-Hastings 
		* algorithm is implemented here.
		*	\param current a pointer of LATTICE_NODE.
		* \return a pointer of LATTICE_NODE. 
	*/
  LATTICE_NODE* get_next(LATTICE_NODE* current) const {
    int total=current->_neighbor_prob.size();
#ifdef PRINT 
   std::copy(current->_neighbor_prob.begin(), current->_neighbor_prob.end(), ostream_iterator<double>(cout," "));
    cout << endl;
#endif
    vector<double> prob(total+1);
    prob[0]=current->_neighbor_prob[0];
    for (int i=1; i<total; i++) {
      prob[i]=prob[i-1]+current->_neighbor_prob[i];
    }
    assert(prob[total-1]<=1.00001);
    prob[total]=1;
#ifdef PRINT
    std::copy(prob.begin(), prob.end(), ostream_iterator<double>(cout," "));
    cout << endl;
#endif
    int idx;
    do {
      idx = randomWithDiscreteProbability(prob);
    } while (idx == total);
//    cout << "returning with:" << idx << endl; 
    return current->_neighbors[idx];
  }

	/*! \fn LATTICE_NODE* create_lattice_node(PAT*& p) 
 		*  \brief A member function to create a new lattice node. 
		*	 It first check whether the pattern p come as parameter is already a lattice node from its canonical code.
		*	 If not a new lattice node is created.  
		*	\param p a reference of a pointer of PAT.
		* \return a pointer of LATTICE_NODE 
	*/
  LATTICE_NODE* create_lattice_node(PAT*& p) {
    const typename PAT::CAN_CODE& cc = check_isomorphism(p);
    p->set_canonical_code(cc);
    std::string min_dfs_cc = cc.to_string();

    LATTICE_NODE* node = exists(min_dfs_cc);
    if (node == 0) {  // new pattern
      node = new LATTICE_NODE(p);
      node->_is_processed = false;
      insert_lattice_node(min_dfs_cc, node);
    }
    else {
      delete p;
      p = node->_pat;
    }
    return node;
  }

	/*! \fn LATTICE_NODE* exists(string p) 
 		*  \brief A member function to check exixtance of a lattice node. 
		*	\param p a string.
		* \return a pointr of LATTICE_NODE.
	*/
  LATTICE_NODE* exists(string p) {;
    CONST_NS_IT it = _node_map.find(p);
    return (it != _node_map.end())? it->second : 0;
  }

	/*! \fn void insert_lattice_node(string p, LATTICE_NODE* ln) 
 		*  \brief A member function to store newly created lattice node. 
		*	\param p a string.
		* \param ln a pointer of LATTICE_NODE.
	*/
  void insert_lattice_node(string p, LATTICE_NODE* ln) {
    _node_map.insert(make_pair(p, ln));
  }

	/*! \fn void process_node(LATTICE_NODE* n) 
 		*  \brief A member function to process a lattice node.
		* This function generates all frequent super and sub patterns of a processed node (n).	
		* It also perform the degree calculation of n as well as of all minned super and sup patterns. 
		*	\param n a LATTICE_NODE pointer
	*/
  void process_node(LATTICE_NODE* n) {
    if (n->_is_processed) return;
    PAT* p = n->_pat;
    assert(p->get_sup_ok() == 0);
//#ifdef PRINT
    cout << "Current pattern:\n";
    cout << *p;
//#endif
    vector<PAT*> neighbors;
		
    _pf->get_freq_super_patterns(p, neighbors); 
		    
		_pf->get_sub_patterns(p, neighbors); 
#ifdef PRINT
    cout << "Its neighbors:" << endl;   
   cout << "Total neighbors="<< neighbors.size() << endl;
#endif
   for (int i=0; i<neighbors.size(); i++) {
      PAT* one_neighbor = neighbors[i];
#ifdef PRINT
      cout << *one_neighbor;
#endif
      int its_degree=_pf->get_super_degree(one_neighbor)+ _pf->get_sub_degree(one_neighbor);
#ifdef PRINT
      cout << "Its degree:" << its_degree << endl;
#endif
      double prob = 1.0 / (its_degree>neighbors.size()? its_degree : neighbors.size());
      LATTICE_NODE* ln = create_lattice_node(one_neighbor);
      int status;
      n->_neighbors.push_back(ln);
      n->_neighbor_prob.push_back(prob);
        
      const typename PAT::CAN_CODE& cc = check_isomorphism(one_neighbor);
      one_neighbor->set_canonical_code(cc);
    }  
    n->_is_processed=true;
  }

  private:
  FREQ_CNT_MAP _freq_cnt; //!< store all sampled frequent itemset patterns.
  NODE_MAP _node_map;	//!< store all lattice node.
  PatternFactory<PAT>* _pf;//!< a PatternFactory object
  int _maxiter;//!< store maximum number of iteration.
  PAT* _last;//!< store last node of the random walk.
  time_tracker tt;
};

#endif
