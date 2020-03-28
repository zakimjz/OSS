//! \file pattern.h - class to represent a graph pattern
#ifndef _PATTERN_H_
#define _PATTERN_H_

#include "adj_matrix.h"
#include "graph_can_code.h"
#include <algorithm>

// For graph pattern that does not have any edge label
//! A Pattern Template Class For graph pattern that does not have any edge label
template <typename V_T>
class Pattern 
{
  public:
    typedef AdjMatrix<V_T> ST;
    typedef std::pair<V_T, V_T> EDGE;
		//! Constructor
    Pattern(const vector<V_T>& labels) { 
      _matrix = new ST(labels);
    }

  private:
    vector<int> _vat;
    multiset<EDGE> _edges;
    ST* _matrix;
}; 

//!forward declaration for friends
template <typename V_T, typename E_T>
class ExPattern;

template <typename V_T, typename E_T>
ostream& operator<< (ostream& ostr, const ExPattern<V_T, E_T>& );


//!A ExPattern Template Class for graph pattern that has both vertex and edge label
template <typename V_T, typename E_T>
class ExPattern
{

  friend ostream& operator<< <>(ostream& osr, const ExPattern<V_T, E_T>&);

  public:
    typedef ExPattern<V_T, E_T> PAT; //!< Type definition of PAT
    typedef V_T VERTEX_T; //!< Type definition of VERTEX_T
    typedef E_T EDGE_T; //!< Type definition of EDGE_T
    typedef FullLabelAdjMatrix<V_T, E_T> ST; //!< Type definition of ST
    typedef std::pair<pair<V_T, V_T>, E_T> EDGE; //!< Type definition of EDGE
    typedef typename multiset<EDGE>::iterator EDGE_IT;
    typedef typename multiset<EDGE>::const_iterator EDGE_CIT;
    typedef ::canonical_code<V_T, E_T> CAN_CODE; //!< Type definition of CAN_CODE
    typedef typename CAN_CODE::INIT_TYPE CC_INIT_TYPE;
    typedef typename CAN_CODE::COMPARISON_FUNC CC_COMPARISON_FUNC;
    typedef set<pair<V_T, E_T> > NEIGHBORS;//!< Type definition of NEIGHBORS
    typedef typename NEIGHBORS::const_iterator NCIT;
    typedef typename NEIGHBORS::iterator NIT;
    
		/*! \fn ExPattern() 
 		*  \brief A Constructor.
 		*/
    ExPattern() {
      _matrix = 0;
      _sup_ok = -1;
      _is_frequent = false;
      _code_known=false;
      _status_known=false;
    }

		/*! \fn ExPattern(const vector<V_T>& labels) 
 		*  \brief A Constructor.
 		*  \param labels a constant reference of V_T vector.
 		*/
    ExPattern(const vector<V_T>& labels) {
      //cout << "At the pattern constructor " << endl;
      _matrix = new ST(labels);
      _sup_ok = -1;
      _is_frequent = false;
      _code_known = false;
      _status_known=false;
    }

		/*! \fn ~ExPattern() 
 		*  \brief A Destructor.
 		*/
   ~ExPattern() {
     if (_matrix) {
       delete _matrix;
     }
     _matrix=0;
   }

		/*! \fn void add_edge(unsigned int i, unsigned int j, E_T e) 
 		*  \brief A member function to add an edge to a pattern.
 		*  \param i,j an unsigned integer.
		*  \param e a E_T type.
 		*/
   void add_edge(unsigned int i, unsigned int j, E_T e) {
      _matrix->add_edge(i,j,e);
      V_T v1 = _matrix->label(i);
      V_T v2 = _matrix->label(j);
      pair<V_T, V_T> v_pair = (v1 < v2) ? make_pair(v1, v2) : make_pair(v2,v1);
      EDGE edge = make_pair(v_pair, e);
      _edges.insert(edge);
      _sup_ok = 1;  // super-set of exact VAT
    }

    PAT* make_null_pattern(int max_id_in_vat) const {
      PAT* clone = new PAT();
      clone->_vat.resize(max_id_in_vat); 
      for (int i=0; i<=clone->_vat.size();i++)
        clone->_vat[i]=i;
      clone->_sup_ok=0;
      clone->_is_frequent=true;
      clone->_matrix=0;
      return clone;
    }

		/*! \fn PAT* clone() const 
 		*  \brief A member function to clone (making copy) a pattern.
		*		\return a PAT type pointer
 		*/
    PAT* clone() const {
#ifdef PRINT
      cout << "In clone:\n";
#endif
      PAT* clone = new PAT();
      clone->_matrix = this->_matrix->clone();
      clone->_vat = _vat; 
      clone->_edges = _edges;
      clone->_is_frequent= _is_frequent;
      clone->_removable_edge_known = false;
      clone->_sup_ok = _sup_ok;
      clone->_status_known=false;
#ifdef PRINT
      cout << "Got the new pattern, returing from clone ..." << endl;
#endif
      return clone;
    }
    
		/*! \fn E_T get_edge_label(unsigned int i, unsigned int j)
 		*  \brief A member function returning edge label of an edge define by two end vectex i and j.
		* \param i, j an unsigned integer
		*		\return a E_T
 		*/
    E_T get_edge_label(unsigned int i, unsigned int j) const {
      return _matrix->get_edge_label(i, j);
    }
 
		/*! \fn V_T label(unsigned int i) const
 		*  \brief A member function returning vertex label of a vertex.
		* \param i an unsigned integer
		*		\return a V_T
 		*/
    V_T label(unsigned int i) const {
      return _matrix->label(i);
    }

		/*! \fn size_t size() const 
 		*  \brief A member function return a pattern size/length.
		*		\return a size_t type value.
 		*/
    size_t size() const {
      if (_matrix)
        return _matrix->size();
      else
        return 0;
    }

		/*! \fn void get_vids_for_this_label(const V_T& v, vector<int>& ret_val) const
 		*  \brief A member function returning vertex ids with the same label as vertex v has.
		* \param v a constant reference of V_T.
		* \param ret_val a reference of integer vector.
 		*/
    void get_vids_for_this_label(const V_T& v, vector<int>& ret_val) const {
      _matrix->get_vid_for_this_label(v, ret_val);
    }

		/*! \fn EDGE remove_edge(const int&a, const int& b) 
 		*  \brief A member function to remove an edge from a pattern.
 		*  \param a,b a constant integer.
		*  \return EDGE
 		*/
    EDGE remove_edge(const int&a, const int& b) {  // returning the removed edge
      V_T v1 = this->label(a);
      V_T v2 = this->label(b);
      E_T e = _matrix->get_edge_label(a,b);

      _matrix->remove_edge(a,b);
      pair<V_T, V_T> v_pair = (v1 < v2) ? make_pair(v1, v2) : make_pair(v2,v1);
      EDGE edge = make_pair(v_pair, e);
      EDGE_IT eit = _edges.find(edge);
      if (eit == _edges.end()) {
        cout << "ERROR: request for non-existent edge removal" << endl;
        exit(1);
      }
      _edges.erase(eit);
      _sup_ok = -1; //current VAT is sub-set of exact VAT
      return edge;
    }

/*    PAT* remove_multiple_edges(const set<EDGE>& edge_set) {
      vector<pair<pair<int,int>, E_T > > to_be_removed;
      AdjIterator iter(_matrix);
      for (iter.first(); !iter.is_done(); iter.next()) {
        pair<int, int> vid_p = iter.current();
        E_T e=this->get_edge_labels(vid_p.first,vid_p.second);
        V_T v1 = _matrix->label(vid_p.first);
        V_T v2 = _matrix->label(vid_p.second);
        pair<V_T, V_T> v_pair = (v1 < v2) ? make_pair(v1, v2) : make_pair(v2,v1);
        EDGE edge = make_pair(v_pair, e);
        if (edge_set.find(edge) == edge_set.end()) {
          to_be_removed.push_back(make_pair(v_pair,e));
        }
      }
      PAT* ret_val=new PAT();
      ret_val->_matrix=_matrix->remove_multiple_edges(to_be_removed);
      if (ret_val->_matrix==0) return 0;
      typedef typename PAT::EDGE_LABEL_CIT EDGE_LABEL_CIT;
      EDGE_LABEL_CIT cit=ret_val->_matrix->_elabel.begin();
      for (;cit != ret_val->_matrix->_elabel.end();cit++) {
        pair<int,int> vid_pairs=cit->first;
        E_T elabel= cit->second;
        V_T v1=ret_val->_matrix->_vlabel[vid_pairs.first];
        V_T v2=ret_val->_matrix->_vlabel[vid_pairs.second];
        EDGE edge = make_pair(make_pair(v1,v2), elabel);
        ret_val->_edges.insert(edge);
      }
      ret_val->_sup_ok = -1; // current-vat is subset of exact VAT
      return ret_val;
    }
  */ 
		
		/*! \fn unsigned int find_removable_edges() 
 		*  \brief A member function to find possible number of removable an edge from a pattern.
		*  \return an unsigned integer
 		*/
    unsigned int find_removable_edges() {
      if (_removable_edge_known) return _removable_edges.size();
      AdjIterator iter(_matrix);
      for (iter.first(); !iter.is_done(); iter.next()) {
        pair<int, int> edge = iter.current();
        if (! _matrix->essential_edge(edge.first, edge.second)) {
          // cout << "edge (" << edge.first << "," << edge.second << ") is removable\n";
          _removable_edges.push_back(edge); 
        }
      }
      _removable_edge_known = true;
      return _removable_edges.size();
    }

		/*! \fn const vector<pair<int, int> >& get_removable_edges() const 
 		*  \brief A member function to get possible removable an edges from a pattern.
		*  \return a constant vector of <int,int> pair
 		*/
    const vector<pair<int, int> >& get_removable_edges() const {
      return _removable_edges;
    } 

		/*! \fn const ST* get_adj_matrix() const  
 		*  \brief A member function to get _matrix private variable.
		*		\return a constant pointer of ST.
 		*/
    const ST* get_adj_matrix() const {
      return _matrix;
    }

		/*! \fn void set_vat(const vector<int>& v) 
 		*  \brief A member function to set VAT of a pattern.
		*		\param v a constant reference of integer type vector.
 		*/
    void set_vat(const vector<int>& v) {
      _vat = v;
    }

		/*! \fn void set_sup_ok() 
 		*  \brief A member function to set  private variable _sup_ok with i.
 		*/
    void set_sup_status(int i ) {
      _sup_ok=i;
    }

		/*! \fn void set_sup_ok() 
 		*  \brief A member function to set _sup_ok private variable.
 		*/
    int get_sup_ok() const {
      return _sup_ok;
    }

		/*! \fn bool edge_exist(int i, int j) const 
 		*  \brief A member function to check edge existance.
		*  \param i,j integer
		* \return boolean
 		*/
    bool edge_exist(int i, int j) const {
      return _matrix->at(i,j);
    }

		/*! \fn int add_vertex(V_T v) 
 		*  \brief A member function to add vertex in Pattern.
		*  \param v a V_T
		* \return integer
 		*/
    int add_vertex(V_T v) {
      return _matrix->add_vertex(v);
    }

		/*! \fn bool is_frequent() const 
 		*  \brief A member function to get _is_frequent private variable.
		*		\return a boolean
 		*/
    bool is_frequent() const {
      assert(_sup_ok == true);
      return _is_frequent;
    }
    
		/*! \fn void set_canonical_code(const CAN_CODE& c) 
 		*  \brief A member function to set canonical code of a pattern.
		*		\param c a constant reference of CAN_CODE
 		*/
    void set_canonical_code(const CAN_CODE& c) {
      _canonical_code = c; 
    }

		/*! \fn const multiset<EDGE>& get_edgeset() 
 		*  \brief A member function to get edge set of a pattern.
		*		\return a constant reference of EDGE multiset.
 		*/
    const multiset<EDGE>& get_edgeset() const {
      return _edges;
    }

		/*! \fn bool is_super_pattern(PAT* pat) 
 		*  \brief A member function to check whether a pattern pat is super pattern of current pattern.
		* \param pat, a pointer of PAT
		*		\return a boolean.
 		*/
    bool is_super_pattern(PAT* pat) {
      const multiset<EDGE>& mset = pat->get_edgeset();
      if (includes(_edges.begin(), _edges.end(), mset.begin(), mset.end())) {
        return true;
      }
      else {
        return false;
      }
    }

		/*! \fn const vector<int>& get_vat() const 
 		*  \brief A member function to get _vat private variable.
		*		\return a constant reference of integer type vector.
 		*/
    const vector<int>& get_vat() const {
      return _vat;
    }
 
		
    int get_vat_size() const {
      return _vat.size();
    }

		/*! \fn void print_vat() const 
 		*  \brief A member function to print VAT of a pattern.
 		*/
    void print_vat() const {
      cout <<"\nVAT:Size(" << _vat.size() << ")\n";
      std::copy(_vat.begin(), _vat.end(), ostream_iterator<int>(cout," "));
      cout << endl;
    }

		/*! \fn void join_vat(PAT* p) 
 		*  \brief A member function to join VAT of a pattern to its super/sub pattern using set_intersection STL algorithm.
		*		\param p a PAT type pointer.
 		*/
    void join_vat(PAT* p) {
#ifdef PRINT
      cout << "Before join: vat size:" << _vat.size() << endl;
#endif
      vector<int> out_vector;
      set_intersection(_vat.begin(), _vat.end(), p->_vat.begin(), p->_vat.end(),
                back_inserter(out_vector));
      _vat = out_vector;
#ifdef PRINT
      cout << "After join: vat size:" << _vat.size() << endl;
#endif
    }

    void set_freq() {
      _sup_ok=0;
      _is_frequent = true;
    }

		/*! \fn std::string canonical_code() const 
 		*  \brief A member function to get canonical code of a pattern.
		*		\return CAN_CODE reference.
 		*/
    const CAN_CODE& get_canonical_code() const { return _canonical_code;} 
 
		/*! \fn void add_tid_to_vat(int tid) 
 		*  \brief A member function to add graph id to a vat of a pattern.
		*		\param tid an integer.
 		*/
    void add_tid_to_vat(int tid) {
      vector<int>::iterator it;
      it = lower_bound(_vat.begin(), _vat.end(), tid);
      _vat.insert(it, tid);
    }

    //void add_tid_to_vat(vector<int>::iterator it, int tid) const {
      //_vat.insert(it, tid);
    //}

		/*! \fn int edge_counter(const EDGE& e) const 
 		*  \brief A member function to count existance of an edge e in a pattern.
		*		\param e a constant reference of EDGE.
		* \return an integer
 		*/
    int edge_counter(const EDGE& e) const {
      pair<EDGE_IT, EDGE_IT> eit_p = equal_range(_edges.begin(), _edges.end(), e);
      int cnt=0;
      for (;eit_p.first != eit_p.second;eit_p.first++) cnt++;
      return cnt;
    }

    int sub_pat_cnt() const {
      return _removable_edges.size();
    }

    int edge_cnt() const {
      return _edges.size();
    }
   
    bool is_code_known() const {return _code_known;} //!< getter of _code_known
    void set_code_known() {_code_known=true;}//!< setter of _code_known

    bool is_status_known() const {return _status_known;} //!< getter of _status_known
    void set_status_known() {_status_known=true;}//!< setter of _status_known

    int get_good_for_class() const {return _good_for_class;} //!< getter of _good_for_class
    void set_good_for_class(int val) {//!< getter of _good_for_class
      _good_for_class = val;
    }

    void set_neighbor_cnt(int x) {_neighbor_count=x;}//!< setter of _neighbor_count
    int get_neighbor_cnt() { return _neighbor_count;}//!< getter of _neighbor_count

    bool _is_frequent;
  private:

    vector<int> _vat;           //!< contains the graph ids that has at least these edges
    multiset<EDGE> _edges; //!< store edges of a pattern
    vector<pair<int, int> > _removable_edges;  //!< removing these edges keeps the pattern connected
    ST* _matrix;
    int _sup_ok;               //!< 0=exact, +1: Super-Set of exact VAT, -1: Sub-Set of exact VAT
    CAN_CODE _canonical_code;
    bool _removable_edge_known;
    bool _code_known;
    bool _status_known;
    int _good_for_class;  
    int _neighbor_count;
};

//! << operator overloading template 
template<typename V_TYPE, typename E_TYPE>
ostream& operator<< (ostream& ostr, const ExPattern<V_TYPE, E_TYPE>& p){

  typedef typename ExPattern<V_TYPE, E_TYPE>::EDGE EDGE;
  if (p.size() == 0) {
    cout << "NULL" << endl;
    cout << "Support:" << p._vat.size() << endl;
    return ostr;
  }
  ostr<< "SIZE:" <<  p.size() << endl;
  ostr<< "EDGE COUNT:" << p._edges.size() << endl;
  ostr<< "Vlabel & ADJ LIST:\n";
  ostr<< *(p._matrix);
  if (p._sup_ok == 0) { 
    ostr<< "Support available, value:" << p._vat.size() << endl;
  }
  else if (p._sup_ok==1){ 
    ostr<< "Support NOT available yet, status:+1" << endl; 
  }
  else { 
    ostr<< "Support NOT available yet, status:-1" << endl; 
  }
  return ostr;
}

#endif
