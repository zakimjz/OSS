//! \file matrix_base.cpp 
#include "matrix_base.h"

using namespace std;


//! default constructor
Matrix::Matrix() {
  _row = _col = 0;
} 

//! constructor that allocates memory, and initialize all cell of the matrix as 0
Matrix::Matrix(size_t r, size_t c) {
  _row = r; _col = c;
  BITVECTOR one_row(c);
  _data.insert(_data.end(), r, one_row);
}

//! Destructor
Matrix::~Matrix() {
  //cout << "In this destructor" << endl;
  //cout << "row:" << _row << " and col:" << _col << endl;
}

const boost::dynamic_bitset<>& Matrix::operator[](size_t i) const {
  return _data[i];
}


/*! \fn void Matrix::allocate(size_t r, size_t c) 
 		*  \brief allocate memory for a matrix of size r X c initialize all entry with 0. 
		*	\param r,c  a size_t.
		
	*/
void Matrix::allocate(size_t r, size_t c) {
  Matrix(r, c);
}

//! 
/*! \fn bool Matrix::at(const unsigned int& r, const unsigned int& c) const 
 		*  \brief Return an element of the matrix [r][c]. 
		*	\param r,c a constant reference of unsigned integer.
		* \return boolean.
	*/
bool Matrix::at(const unsigned int& r, const unsigned int& c) const {
  if ((r > _row) || (c > _col)) return 0;
  return _data[r][c];
}

/*! \fn void Matrix::reset(const unsigned int& r) 
 		*  \brief A member function to reset the matrix. 
		*	\param r a constant reference of unsigned integer.
	*/
void Matrix::reset(const unsigned int& r) {
  _data[r].reset();
}
   
// 
/*! \fn void Matrix::set(const unsigned int& r, const unsigned int& c, bool val) 
 		*  \brief set the value of the element [r][c]. 
		*	\param r,c a constant reference of unsigned integer.
		* \param val a boolean.
		*/
void Matrix::set(const unsigned int& r, const unsigned int& c, bool val) {
  _data[r].set(c,val);
}

/*! \fn size_t Matrix::rowset_cnt(unsigned int i) 
 		*  \brief A member function to count non zero element in a row. 
		*	\param i an unsigned integer.
		* \return size_t.
	*/
size_t Matrix::rowset_cnt(unsigned int i){
  return _data[i].count();
}

/*! \fn bool Matrix::rowset_empty(unsigned int i) const 
 		*  \brief A member function to check empty row of a matrix. 
		*	\param i an unsigned integer.
		* \return boolean.
	*/
bool Matrix::rowset_empty(unsigned int i) const {
  if (_data[i].none()) return true;
  return false;
}

/*! \fn void Matrix::neighbors(const unsigned int& i, vector<unsigned int>& ret_val) const 
 		*  \brief A member function to find neighbors. 
		*	\param i a constant reference of unsigned integer.
		* \param ret_val a reference of unsigned integer vector.
	*/
void Matrix::neighbors(const unsigned int& i, vector<unsigned int>& ret_val) const {
  if (i > _row) {
    cout << "ERROR: Matrix::Requesting degree of non-existent vertex" << endl;
    exit(1);
  }
  for (size_t pos = _data[i].find_first(); pos != BITVECTOR::npos; pos = _data[i].find_next(pos))
    ret_val.push_back(pos);
}

/*! \fn Matrix operator*(const Matrix& A, const Matrix& B) 
 		*  \brief A member function to overload * operator for matrix multiplication. 
		*	\param A,B a constant reference of Matrix.
		* \return Matrix
	*/
Matrix operator*(const Matrix& A, const Matrix& B) { 
  assert(A.col() == B.row());
  unsigned int r = A.row();
  unsigned int c = B.col();

  Matrix C(A.row(), B.col());
  for (int i = 0; i < r; i++) 
    for (int j = 0; j < c; j++) 
      for (int k= 0; k < A.col(); k++)
        C._data[i][j] = C._data[i][j] + A.at(i,k) * B.at(k,j);
  return C; 
} 

/*! \fn Matrix Transpose(const Matrix& A) 
 		*  \brief A member function to get transpose of a Matrix. 
		*	\param A a constant reference of Matrix.
		* \return Matrix.
	*/
Matrix Transpose(const Matrix& A) { 
  Matrix B(A.col(), A.row());
  for (int i = 0; i < B.row(); i++)
    for (int j = 0; j < B.col(); j++)
      B._data[i][j] = A.at(j,i);
  return B;
}

/*! \fn ostream& operator<< (ostream& ostr, const Matrix& M) 
 		*  \brief A member function to overload << operator for matrix print. 
	*/
ostream& operator<< (ostream& ostr, const Matrix& M){

  for (int i = 0; i < M.row(); i++)
    for (int j = 0; j < M.col(); j++) {
      ostr << (M.at(i,j)? '1' : '0');
      if (j < M.col()-1)
        ostr << " ";
      else
        ostr << "\n";
    }
  return ostr;
}

SqrSymMatrix::SqrSymMatrix() : Matrix() { }

SqrSymMatrix::SqrSymMatrix(size_t n) : Matrix(n,n) { 
  _size = n;
}

bool SqrSymMatrix::at(const unsigned int& r, const unsigned int& c) const {
  if ((r > _size) || (c > _size)) return 0;
  return _data[r][c];
}

void SqrSymMatrix::neighbors(const unsigned int& i, vector<unsigned int>& ret_val) const {
  if (i > _size) {
    cout << "ERROR: Matrix::Requesting neighbors of non-existent vertex" << endl;
    exit(1);
  }
  for (size_t pos = _data[i].find_first(); pos != BITVECTOR::npos; pos = _data[i].find_next(pos))
    ret_val.push_back(pos);
}

/*! \fn int SqrSymMatrix::add_vertex() 
 		*  \brief A member function to increase row size because of new vertex addition. 
		
		* \return an integer.
	*/
int SqrSymMatrix::add_vertex() {
  for (int i=0; i<_size;i++) {
    _data[i].resize(_size+1);
  }
  BITVECTOR one_row(_size+1);
  _data.push_back(one_row);
  _size = _size+1;
  return _size-1; // returning vertex-id
}

/*! \fn void SqrSymMatrix::change_adj_matrix(int size, const vector<pair<int, int> >& adj_list) 
 		*  \brief A member function to update adjacency list of a graph patten when it is extended or reduced in size. 
		*	\param size an integer.
		*	\param adj_list a constant reference to pair<int,int> vector.
	*/
void SqrSymMatrix::change_adj_matrix(int size, const vector<pair<int, int> >& adj_list) {
  _data.clear();
  BITVECTOR one_row(size);
  _data.insert(_data.end(), size, one_row);
  for (int i=0; i< adj_list.size(); i++) {
    pair<int, int> one_edge = adj_list[i];
    set(one_edge.first, one_edge.second, 1);
    set(one_edge.second, one_edge.first, 1);
  }
  _size = size;
}

/*! \fn bool SqrSymMatrix::essential_edge(int src, int dst) const 
 		*  \brief A member function to check exixtance of a path from scr to dst vetex of a graph pattern. 
		*	\param src,dst an integer.
		* \return boolean.
	*/
bool SqrSymMatrix::essential_edge(int src, int dst) const {
  BITVECTOR visited(_size);
  if (_data[src].count() == 1 || _data[dst].count() == 1) return false;
  visited[src] = true;
  bool can_reach_dst = false;
  for (size_t pos = _data[src].find_first(); pos != BITVECTOR::npos;
        pos = _data[src].find_next(pos)) {
    if (pos != dst) {
      can_reach_dst = dfs_visit(pos, visited, dst);
      if (can_reach_dst) return false;  // this edge is not essential
    }
  }
  assert(visited[dst] == false);
  return true;
}
 
/*! \fn bool SqrSymMatrix::dfs_visit(int src, BITVECTOR& visited, int dst) const 
 		*  \brief A member function to make a dfs visit starting from scr. 
		*	\param src,dst an integer.
	*	\param visited a reference of BITVECTOR.
		* \return a boolean.
	*/
bool SqrSymMatrix::dfs_visit(int src, BITVECTOR& visited, int dst) const {
  visited[src] = true;
  for (size_t pos = _data[src].find_first(); pos != BITVECTOR::npos;
       pos = _data[src].find_next(pos)) {
    if (visited[pos] == true) continue;
    if (pos == dst) return true;
    bool ret_val = dfs_visit(pos, visited, dst);
    if (ret_val == 1) return true;
  }
  return false;
}

//! Constructor
AdjIterator::AdjIterator(const SqrSymMatrix* m) {
  _m = m;
  _i = 0;
  _j = 0;
  _is_done = false;
  _max_i = m->_data.size()-1;
}

/*! \fn void AdjIterator::first() 
 		*  \brief A member function to set iterator at first position. 
		
	*/
void AdjIterator::first() {
  _i = 1;
  _j = -1;
  next();
}

/*! \fn void AdjIterator::next()  
 		*  \brief A member function to get next one in _data BITVECTOR. 
		
		
	*/
void AdjIterator::next() {
  while (true) {
    size_t p = (_j == -1) ? _m->_data[_i].find_first() :
               _m->_data[_i].find_next(_j);
    if (p > _i) { 
      _i++;
      _j = -1;
      if (_i > _max_i) {
        _is_done = true;
        break;
      } 
    }
    else if (p == boost::dynamic_bitset<>::npos) {
      _is_done = true;
      break;
    }
    else {
      _j = p;
      break;
    }
  } 
}

/*! \fn pair<size_t, size_t> AdjIterator::current() const 
 		*  \brief A member function to return the current position of iterator. 
		
		* \return pair<size_t,size_t>.
	*/
pair<size_t, size_t> AdjIterator::current() const {
  if (_is_done) {
    cout << "ERROR:Iterator bound error" << endl;
    exit(1);
  }
  return make_pair(_j, _i);
}

bool AdjIterator::is_done() const {
  return _is_done;
}

/*! \fn NodeAdjIterator::NodeAdjIterator(const SqrSymMatrix* m, size_t i) 
 		*  \brief Constructor. 
		*	\param m a constant pointer of SqrSymMatrix.
		*	\param i a size_t.
		
	*/
NodeAdjIterator::NodeAdjIterator(const SqrSymMatrix* m, size_t i) {
  _m = m;
  _i = i;
  _j = 0;
  _is_done = false;
  if (i > m->_data.size()-1) {
    cout << "ERROR: index out of vector limit" << endl;
    exit(1);
  }
}

/*! \fn void AdjIterator::first() 
 		*  \brief A member function to set iterator at first position. 
		
	*/
void NodeAdjIterator::first() {
  _j = -1;
  next();
}

/*! \fn void AdjIterator::next()  
 		*  \brief A member function to get next one in _data BITVECTOR. 
		
		
	*/
void NodeAdjIterator::next() {
  size_t p = (_j == -1) ? _m->_data[_i].find_first() :
               _m->_data[_i].find_next(_j);
  if (p == boost::dynamic_bitset<>::npos)
    _is_done = true;
  else
    _j = p;
}

/*! \fn pair<size_t, size_t> AdjIterator::current() const 
 		*  \brief A member function to return the current position of iterator. 
		
		* \return  a size_t.
	*/
size_t NodeAdjIterator::current() const {
  if (_is_done) {
    cout << "ERROR:Iterator bound error" << endl;
    exit(1);
  }
  return _j;
}

bool NodeAdjIterator::is_done() const {
  return _is_done;
}
