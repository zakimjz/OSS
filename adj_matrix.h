//! \file adj_matrix.h - class for adjacency matrix implementation
#ifndef _ADJ_MATIX_H_
#define _ADJ_MATIX_H_

#include <map>
#include <ext/hash_map>
#include <set>
#include <algorithm>
#include <iterator>
#include "helper_funs.h"
#include "matrix_base.h"

using namespace std;

//! forward declarations for friend functions
template <typename V_T>
class AdjMatrix;

// filling out the matching matrix to match a as a 
// subgraph of b, the matching matrix should have 
// size of (size(a), size(b))
template <typename V_T>
void matcher(const AdjMatrix<V_T>& a, const AdjMatrix<V_T>& b, Matrix&);


//! A AdjMatrix Class for graph with labeled vertex
template <typename V_T = int>
class AdjMatrix : public SqrSymMatrix {

  friend void matcher<>(const AdjMatrix<V_T>&, const AdjMatrix<V_T>&, Matrix&);
  
  public:
    typedef AdjMatrix<V_T> MAT_V_LABEL;
		//! Default Constructor
    AdjMatrix() : SqrSymMatrix() { }

		//!Constructor
    AdjMatrix(size_t n) : SqrSymMatrix(n) { 
      _vlabel.resize(n);
      _size = n;
    }
 		//! Constructor
    AdjMatrix(const vector<V_T>& labels) : SqrSymMatrix(labels.size()) {
      std::copy(labels.begin(), labels.end(), back_inserter(_vlabel));
      _size = _vlabel.size();
    }
		
		/*! \fn void neighbors(const unsigned int& i, vector<int>& ret_val) const 
 		*  \brief A member function to handle neighbors. 
		*	\param i a constant reference of unsigned integer.
		* \param ret_val a reference of integer vector.
	*/
    void neighbors(const unsigned int& i, vector<int>& ret_val) const {
      if (i > _size) {
        cout << "ERROR: Matrix::Requesting degree of non-existent vertex:" << i << endl;
        exit(1);
      }
      for (size_t pos = _data[i].find_first(); pos != BITVECTOR::npos; pos = _data[i].find_next(pos))
        ret_val.push_back(pos);
    }

		/*! \fn unsigned int degree(const unsigned int& i) const 
 		*  \brief A member function to get degree of a veretx. 
		*	\param i a constant reference of integer.
		* \return unsigned integer.
	*/
    unsigned int degree(const unsigned int& i) const {
      if (i > _size) {
        cout << "ERROR: AdjMatrix::Requesting degree of non-existent vertex " << i << endl;
        exit(1);
      }
      return _data[i].count(); 
    }

		/*! \fn void set_vlabel(unsigned int i, V_T l) 
 		*  \brief A member function to set label of a vertex. 
		*	\param i an unsigned integer.
		*	\param l a V_T.
		* \return boolean.
	*/
    void set_vlabel(unsigned int i, V_T l) {
      if (_size > i) 
        _vlabel[i] = l;  
      else
        cout << "ERROR: AdjMatrix::set_vlabel--> vertex label assignment on invalid position" << endl;
    }

		/*! \fn const V_T& label(unsigned int i) const 
 		*  \brief A member function to get a label of a vertex. 
		*	\param i an unsigned integer.
		* \return constant reference of V_T.
	*/
    const V_T& label(unsigned int i) const {
      if (_size > i)
        return _vlabel[i]; 
      else {
        cout << "ERROR: AdjMatrix::label--> requested vertex label from invalid position" << endl;
        cout << "POS:=" << i << " SIZE:=" << _size << endl;
        exit(1);
      }
    }

		/*! \fn void add_edge(unsigned int i, unsigned int j) 
 		*  \brief A member function to add a new edge. 
		*	\param i,j an unsigned integer.
	*/
    void add_edge(unsigned int i, unsigned int j) {
      this->set(i,j,1);
      this->set(j,i,1);
    }


		/*! \fn void get_vid_for_this_label(V_T v, vector<int>& vec) 
 		*  \brief A member function to get vertex ids for a label. 
		*	\param v a V_T.
		*	\param vec a reference of integer vector.
	*/
    void get_vid_for_this_label(V_T v, vector<int>& vec) {
      for (int i=0; i<_size; i++) {
        if (_vlabel[i] == v) vec.push_back(i);
      }
    }

    inline unsigned int get_edge_count() const {return _edge_cnt;};

		/*! \fn int add_vertex(V_T l) 
 		*  \brief A member function to add a new vertex. 
		*	\param l a V_T.
		* \return integer.
	*/
    int add_vertex(V_T l) {
      int ret_val = SqrSymMatrix::add_vertex();
      _vlabel.push_back(l);
      return ret_val;
    }
  protected:

    vector<V_T> _vlabel;
    size_t _edge_cnt;
};
/*! A Template function matcher.
	* \brief This function build "M" matrix used in Ullman's subgraph isomorphism algorithm. Function applies label invariant and degree invariant to make less 
number of one(1) in M maintaining the restriction.
*/
template<typename V_T>
void matcher(AdjMatrix<V_T>& a, AdjMatrix<V_T>& b, Matrix& M) {
  int size1 = a.size(); int size2 = b.size();
	cout << size1 <<size2<<endl;
  if (size1 > size2) return;
	cout<< "reached here" <<endl;
  // using vertex label invariant
  for (int i = 0; i < size1; i++) {
    V_T v1 = a.label(i);
    for (int j = 0; j < size2; j++) {
      V_T v2 = b.label(j);
      if (v1 == v2) M.set(i,j,1);
    }
  }
  cout << "After vertex label invariant\n";
  cout << M << endl;
  // using own degree invariant
  vector<pair<int, int> > match_pair;
  for (int i = 0; i < size1; i++) {
    if (M.rowset_cnt(i)  == 0) return;
    if (M.rowset_cnt(i) == 1) continue;
    unsigned int d1 = a.degree(i);
    for (int j = 0; j < size2; j++) {
      if (M.at(i,j) == 0) continue;
      unsigned int d2 = b.degree(j);
      if (d1 > d2)
        M.set(i,j,0);
      else 
        match_pair.push_back(make_pair(i,j)); 
    }
  }

  cout << "After own degree invariant\n";
  cout << M << endl;
  // using neighbor-degree and label invariant
  for (int idx=0; idx < match_pair.size(); idx++) {
    int i = match_pair[idx].first;
    int j = match_pair[idx].second;
    vector<int> neighbors; 
    a.neighbors(i, neighbors);
    vector<int> neigh_degree1(neighbors.size());
    multiset<V_T> neigh_labels1;
 
    for (int k=0; k < neighbors.size(); k++) {
      unsigned int d_n = a.degree(neighbors[k]);
      neigh_labels1.insert(a.label(neighbors[k]));
      neigh_degree1.push_back(d_n); 
    }
    
    vector<int> neighbors2; 
    b.neighbors(j, neighbors2);
    vector<int> neigh_degree2(neighbors2.size());
    multiset<V_T> neigh_labels2;
    for (int k=0; k < neighbors2.size(); k++) {
      unsigned int d_n = b.degree(neighbors2[k]);
      neigh_labels2.insert(b.label(neighbors2[k]));
      neigh_degree2.push_back(d_n); 
    }
    //std::copy(neigh_labels1.begin(), neigh_labels1.end(), ostream_iterator<char>(cout, " "));
    //cout << endl;
    //std::copy(neigh_labels2.begin(), neigh_labels2.end(), ostream_iterator<char>(cout, " "));
    //cout << endl;

    // first applying neighbor-label invariant
    if (!includes (neigh_labels2.begin(), neigh_labels2.end(), neigh_labels1.begin(),
        neigh_labels1.end())) {
      cout << "setting " << i << " " << j << " to 0" << endl;
      M.set(i,j,0);
      continue;
    }

    // now applying neighbor-degree invariant
    sort(neigh_degree1.begin(), neigh_degree1.end(), greater<int>());
    sort(neigh_degree2.begin(), neigh_degree2.end(), greater<int>());

    int k = 0, l = 0;
    for (; k < neigh_degree1.size(), l < neigh_degree2.size(); ) {
      if (neigh_degree1[k] <= neigh_degree2[l]) {k++; l++;}
      else 
        break;
    }
    if (k < neigh_degree1.size()) M.set(i, j, 0); 
  }
  cout << "at end\n";
  cout << M << endl;
}

// forward declarations for friend functions
template <typename V_T, typename E_T>
class FullLabelAdjMatrix;

// filling out the matching matrix to match a as a 
// subgraph of b, the matching matrix should have 
// size of (size(a), size(b))
template <typename V_T, typename E_T>
void matcher(const FullLabelAdjMatrix<V_T, E_T>& a, const FullLabelAdjMatrix<V_T, E_T>& b, Matrix&);

template <typename V_T, typename E_T>
ostream& operator<< (ostream& ostr, const FullLabelAdjMatrix<V_T, E_T>& M);


//! A FullLabelAdjMatrix Class for graph with labeled vertex end edge
template <typename V_T = int, typename E_T = int>
class FullLabelAdjMatrix : public AdjMatrix<V_T> {
  
  friend ostream& operator<< <>(ostream& osr, const FullLabelAdjMatrix&);
  friend void matcher<>(const FullLabelAdjMatrix<V_T, E_T>&, const FullLabelAdjMatrix<V_T, E_T>&, Matrix&);

  public:
    typedef V_T VERTEX_T;
    typedef E_T EDGE_T;
    typedef map<pair<int, int>, E_T> EDGE_LABEL;
    typedef typename EDGE_LABEL::iterator EDGE_LABEL_IT;
    typedef typename EDGE_LABEL::const_iterator EDGE_LABEL_CIT;
    typedef FullLabelAdjMatrix<V_T, E_T> MAT_V_E_LABEL;

		//!Default Constructor    
		FullLabelAdjMatrix() : AdjMatrix<V_T>() { }

		//! Constructor
    FullLabelAdjMatrix(size_t n) : AdjMatrix<V_T>(n) { }

		//! Constructor
    FullLabelAdjMatrix(const vector<V_T>& labels) : AdjMatrix<V_T>(labels) {}

		//! Constructor
    FullLabelAdjMatrix(const FullLabelAdjMatrix<V_T, E_T>& rhs) {
      this->_data = rhs._data;
      this->_vlabel = rhs._vlabel;
      this->_size = rhs._size;
      _elabel = rhs._elabel; 
    }

		//! Constructor
    FullLabelAdjMatrix<V_T, E_T>& operator=(const FullLabelAdjMatrix<V_T, E_T> & m) {
      this->_data = m._data;
      this->_vlabel=m._vlabel;
      _elabel.insert(m._elabel.begin(), m._elabel.end());
      return *this;
    }

		/*! \fn FullLabelAdjMatrix<V_T, E_T>* clone() const 
 		*  \brief A member function to make a copy/clone of a matrix. 
		* \return Pointer of FullLabelAdjMatrix
	*/
    FullLabelAdjMatrix<V_T, E_T>* clone() const {
      
      FullLabelAdjMatrix<V_T, E_T>* m_new = new FullLabelAdjMatrix(*this);
      EDGE_LABEL_IT it = m_new->_elabel.begin();
      return m_new;
    }

		/*! \fn void set_edge_label(unsigned int i, unsigned int j, E_T e)
 		*  \brief A member function to add a set edge label. 
		*	\param i,j an unsigned integer.
		*	\param e a E_T.
	*/
    void set_edge_label(unsigned int i, unsigned int j, E_T e) {
       if (i ==j) return;  // no self loop allowed
       pair<EDGE_LABEL_IT, bool> ret_type;
       pair<int, int> key = (i<j) ? make_pair(i,j) :make_pair(j,i);
       ret_type = _elabel.insert(make_pair(key, e));
       if (ret_type.second == false) {
         cout << "Edge label insertion unsuccessful (1)" << endl; 
         exit(1);
       }
    }

		/*! \fn const E_T& get_edge_label(unsigned int i, unsigned int j) const 
 		*  \brief A member function to get the label of an edge. 
		*	\param i,j an unsigned integer.
		* \return constant reference of E_T
	*/
    const E_T& get_edge_label(unsigned int i, unsigned int j) const {
      
      if ((i == j) || (i>=this->_size) || (j>=this->_size)) {
        cout << "ERROR: AdjMatrix::get_edge_label--> requested edge label from invalid position" << endl;
        cout <<"i=" << i << " j:=" << j << " size:=" << this->_size << endl;
        exit(1);
      }
      EDGE_LABEL_CIT it;
      pair<int,int> key = (i < j)? make_pair(i,j) : make_pair(j,i);
      it = _elabel.find(key);
      if (it != _elabel.end()) 
        return it->second;
      else {
        EDGE_LABEL_CIT it2;
        for (it2=_elabel.begin(); it2 != _elabel.end(); it2++) {
          cout << it2->first.first << " " << it2->first.second << " " << it2->second << endl;
        }
        cout <<"i=" << i << " j:=" << j << " size:=" << this->_size << " could not find in edge-map" << endl;
        cout << "ERROR: AdjMatrix::get_edge_label--> requested edge label from invalid position" << endl;
        exit(1);
      }
    }

		/*! \fn void remove_edge(const int&i, const int& j) 
 		*  \brief A member function to remove an edge by removing both of end vertices. 
		*	\param i,j a constant reference of integer.
	*/
    void remove_edge(const int&i, const int& j) {
      this->set(i,j,0);
      this->set(j,i,0);
     
      // both vertex got isolated, happens only when we
      // are deleting an edge from an one-edge pattern
      if (this->_data[i].none() && this->_data[j].none()) {  
        assert(this->_size == 2); // we have only two vertices
        this->_size = 0;
        this->_data.clear();
        this->_vlabel.clear();
        this->_elabel.clear();
        return;
      }
      bool v_removed=false;

      if (this->_data[i].none()) {
        remove_vertex(i);
        v_removed=true;
      }
      else if (this->_data[j].none()) {
        remove_vertex(j);
        v_removed=true;
      }
      if (v_removed==false) {
        pair<int,int> key = (i < j)? make_pair(i,j) : make_pair(j,i);
        _elabel.erase(key);
      }
    }
	
	
		/*! \fn int add_vertex(V_T l) 
 		*  \brief A member function to add a new vertex. 
		*	\param l a V_T.
		* \return integer.
	*/
/*    FullLabelAdjMatrix<V_T, E_T>* remove_multiple_edges(const vector<pair<int, int> >& d_edges ) {
    
      vector<int> to_be_degree(this->_size,-1);  // -1 means remain unchanged, 0 means to be deleted
      for (int k=0; k<d_edges.size(); k++) {
        int i=d_edges[k].first;
        to_be_degree[i]=this->degree(i)-1;
        i=d_edges[k].second;
        to_be_degree[i]=this->degree(i)-1;
      }
      map<int, int> d_map;
      vector<int> rev_map(this->_size);
      int vid_num=0;
      for (int k=0;k<this->_size;k++) {
        if (to_be_degree[k]!=0) {
          d_map.insert(make_pair(vid_num,k));
          rev_map[k]=vid_num;
          vid_num++;
        }
      }
      if (vid_num==0) return 0;
      vector<V_T> labels(vid_num);
      for (int k=0;k<vid_num;k++) {
        labels[k] = this->_vlabel[d_map[k]];
      }
      FullLabelAdjMatrix<V_T,E_T>* ret_val=new FullLabelAdjMatrix<V_T,E_T>(labels);
      AdjIterator iter(this->_matrix);
      for (iter.first(); !iter.is_done(); iter.next()) {
        pair<int, int> edge = iter.current();
        int first = rev_map[edge.first];
        int second = rev_map[edge.second];
        if (first==-1 || second==-1) continue;
        ret_val->add_edge(first,second);
        E_T elabel=get_edge_label(edge.first, edge.second);
        ret_val->_elabel.insert(make_pair(make_pair(first,second),elabel));
      }
    }*/
		
		/*! \fn  void add_edge(unsigned int i, unsigned int j, E_T e) 
 		*  \brief A member function to add a new edge. 
		*	\param i,j a unsigned integer.
		*	\param e a E_T.
	*/
    void add_edge(unsigned int i, unsigned int j, E_T e) {
      this->set(i,j,1);
      this->set(j,i,1);
      set_edge_label(i,j,e);
    }

  private:

		/*! \fn void remove_vertex(size_t i) 
 		*  \brief A member function to remove a vertex from adjacency matrix. 
		*	\param i a size_t.
		
	*/
    void remove_vertex(size_t i) {
      int its_degree = this->degree(i);
      assert(its_degree == 0);
      vector<pair<int, int> > v_pairs;
      v_pairs.reserve(_elabel.size());
      vector<E_T> elabels;
      elabels.reserve(_elabel.size());

      EDGE_LABEL_CIT it = _elabel.begin();
      for(; it != _elabel.end();it++) {
        const pair<int, int>& key = it->first;
        const E_T& elabel = it->second;
        if (key.first < i && key.second < i) {
          v_pairs.push_back(key); 
          elabels.push_back(elabel); 
        }
        else if (key.first < i && key.second > i) {
          v_pairs.push_back(make_pair(key.first, key.second-1));
          elabels.push_back(elabel); 
        }
        else if (key.first > i && key.second > i) {
          v_pairs.push_back(make_pair(key.first-1, key.second-1));
          elabels.push_back(elabel); 
        }
        else if (key.first == i || key.second == i) {
          continue;
        }
        else {
          cout << "ERROR: program should not come to this line" << endl;
          exit(1);
        }
      }
      _elabel.clear();
      for (int j=0; j<v_pairs.size(); j++) {
        _elabel.insert(make_pair(v_pairs[j], elabels[j]));
      }
      this->_vlabel.erase(this->_vlabel.begin()+i);
      SqrSymMatrix::change_adj_matrix(this->_size-1, v_pairs);
    }

    EDGE_LABEL _elabel;
};

/*! A Template function matcher.
	* \brief This function build "M" matrix used in Ullman's subgraph isomorphism algorithm. Function applies label invariant and degree invariant to make less 
number of one(1) in M maintaining the restriction.
*/
template<typename V_T, typename E_T>
void matcher(const FullLabelAdjMatrix<V_T, E_T> & a,const  FullLabelAdjMatrix<V_T, E_T> & b, Matrix& M) {

//  cout << "a:" << endl;
//  cout << a << endl;
//  cout << "b:" << endl;
//  cout << b << endl; 
  int size1 = a.size(); int size2 = b.size();
//	cout <<size1<<size2<<endl;
  if (size1 > size2) return;

  // using vertex label invariant
  for (int i = 0; i < size1; i++) {
    V_T v1 = a.label(i);
    for (int j = 0; j < size2; j++) {
      V_T v2 = b.label(j);
      if (v1 == v2) M.set(i,j,1);
    }
  }
//  cout << "After vertex label invariant\n";
//  cout << M << endl;
  // using own degree invariant
  vector<pair<int, int> > match_pair;
  for (int i = 0; i < size1; i++) {
    if (M.rowset_cnt(i)  == 0) return;
    if (M.rowset_cnt(i) == 1) continue;
    unsigned int d1 = a.degree(i);
    for (int j = 0; j < size2; j++) {
      if (M.at(i,j) == 0) continue;
      unsigned int d2 = b.degree(j);
      if (d1 > d2)
        M.set(i,j,0);
      else 
        match_pair.push_back(make_pair(i,j)); 
    }
  }

//  cout << "After own degree invariant\n";
//  cout << M << endl;
  // using neighbor-degree and label invariant
  for (int idx=0; idx < match_pair.size(); idx++) {
    int i = match_pair[idx].first;
    int j = match_pair[idx].second;
    vector<int> neighbors; 
    a.neighbors(i, neighbors);
    vector<int> neigh_degree1;
    multiset<pair<V_T, E_T> > neigh_labels1;
 
    for (int k=0; k < neighbors.size(); k++) {
      unsigned int d_n = a.degree(neighbors[k]);
      E_T e_l = a.get_edge_label(i,neighbors[k]);
      neigh_labels1.insert(make_pair(a.label(neighbors[k]), e_l));
      neigh_degree1.push_back(d_n); 
    }
    
    vector<int> neighbors2; 
    b.neighbors(j, neighbors2);
    //cout << j << ", its neighbors:";
    //std::copy(neighbors2.begin(), neighbors2.end(), ostream_iterator<int>(cout, " ")); cout << endl;
  
    vector<int> neigh_degree2;
    multiset<pair<V_T, E_T> > neigh_labels2;
    for (int k=0; k < neighbors2.size(); k++) {
      unsigned int d_n = b.degree(neighbors2[k]);
      E_T e_l = b.get_edge_label(j,neighbors2[k]);
      neigh_labels2.insert(make_pair(b.label(neighbors2[k]), e_l));
      neigh_degree2.push_back(d_n); 
    }
    /*
    typename multiset<pair<V_T, E_T> >::iterator it1;
    for (it1= neigh_labels1.begin(); it1 != neigh_labels1.end(); it1++) {
      cout << it1->first << " " << it1->second << ",";
    }
    cout << endl;
    for (it1= neigh_labels2.begin(); it1 != neigh_labels2.end(); it1++) {
      cout << it1->first << " " << it1->second << ",";
    }
    cout << endl;
    */
    // first applying neighbor-label invariant
    if (!includes (neigh_labels2.begin(), neigh_labels2.end(), neigh_labels1.begin(),
        neigh_labels1.end())) {
      // cout << "(1) setting " << i << " " << j << " to 0\n";
      M.set(i,j,0);
      continue;
    }

    // now applying neighbor-degree invariant
    sort(neigh_degree1.begin(), neigh_degree1.end(), greater<int>());
    sort(neigh_degree2.begin(), neigh_degree2.end(), greater<int>());
    //cout << "DEGREE:\n";
    //std::copy(neigh_degree1.begin(), neigh_degree1.end(), ostream_iterator<int>(cout, " ")); cout << endl;
    //std::copy(neigh_degree2.begin(), neigh_degree2.end(), ostream_iterator<int>(cout, " ")); cout << endl;
    int k = 0, l = 0;
    while (k<neigh_degree1.size() && l<neigh_degree2.size()) {
      if (neigh_degree1[k] <= neigh_degree2[l]) {k++; l++;}
      else 
        break;
    }
    if (k < neigh_degree1.size()) {
      // cout << "(2) setting " << i << " " << j << " to 0\n";
      M.set(i, j, 0); 
    }
  }
//  cout << "at end\n";
//  cout << "Getting out from matcher with M:=\n";
// 	cout << M << endl;
}
/*! \fn ostream& operator<< (ostream& ostr, const FullLabelAdjMatrix<V_T, E_T>& M) 
 		*  \brief A overloading function of operator << to print matrix. 
		*	\param oste a reference of ostream
		*	\param M a constant reference of FullLabelAdjMatrix
		* \return reference of ostream
	*/
template<typename V_T, typename E_T>
ostream& operator<< (ostream& ostr, const FullLabelAdjMatrix<V_T, E_T>& M){
  for (int i = 0; i < M._size; i++){
    ostr << i << "(" << M.label(i) << ")  ";
  }
  ostr << endl;
  ostr <<"-------------------" << endl;
  AdjIterator iter(&M);
  for (iter.first(); !iter.is_done(); iter.next()) {
    pair<int, int> edge = iter.current();
    int first = edge.first;
    int second = edge.second;
    ostr << first << " " << second << " " << M.label(first) << " " 
         << M.get_edge_label(first, second) << " " << M.label(second)<< endl;
  }
 return ostr;
}

#endif
