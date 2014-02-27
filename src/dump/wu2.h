#ifndef _MAIN_H
#define _MAIN_H

//#define _median linasm::median
#define _median partition_median 
#define _lowerquartile partition_lowerquartile
#define _upperquartile partition_upperquartile
#define _interquartilerange partition_interquartilerange

#include <exception>
#include <sstream>
using namespace std;

template<class T>
inline void SWAP(T &a, T &b) {T dum=a; a=b; b=dum;}
template<class T>
inline void Swap(T &a, T &b) {T dum=a; a=b; b=dum;}

// Vector and Matrix Classes

template <class T>
class Array {
 private:
  int nn;	// size of Array. upper index is nn-1
  int nnn;	// capacity of Array. nnn>=nn
  T *v;
 public:
  Array();
  explicit Array(int n);		// Zero-based Array
  Array(int n, const T &a);	        // Initialize to constant value
  Array(int n, const T *a);	        // Initialize to Array
  Array(const Array &rhs);	        // Copy constructor
  Array & operator=(const Array &rhs);	// Assignment
  typedef T value_type;                 // Make T available externally
  inline T & operator[](const int i);	// i'th element
  inline const T & operator[](const int i) const;
  inline int size() const;
  inline int capacity() const;
  void resize(int newn);               // Resize (contents not preserved)
  void reuse(int newn);               // Resize (contents preserved)
  void assign(int newn, const T &a);   // Resize and assign a constant value
  void assign(const T &a);             // assign all elements a constant value
  void assign(int newn, const T *a);   // Resize and assign a constant value
  ~Array();
};

// Array definitions
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

template <class T>
Array<T>::Array() : nn(0), nnn(0), v(NULL) {}

template <class T>
Array<T>::Array(int n) : nn(n), nnn(n), v(n>0 ? new T[n] : NULL) {}

template <class T>
Array<T>::Array(int n, const T& a) : nn(n), nnn(n), v(n>0 ? new T[n] : NULL)
{ for(int i=0; i<n; i++) v[i] = a; }

template <class T>
Array<T>::Array(int n, const T *a) : nn(n), nnn(n), v(n>0 ? new T[n] : NULL)
{ for(int i=0; i<n; i++) v[i] = *a++; }

template <class T>
Array<T>::Array(const Array<T> &rhs) : nn(rhs.nn), nnn(rhs.nnn), v(nnn>0 ? new T[nnn] : NULL)
{  for(int i=0; i<nnn; i++) v[i] = rhs[i]; }

template <class T>
Array<T> & Array<T>::operator=(const Array<T> &rhs)
{    
  if (this != &rhs)
    {
      if (nnn != rhs.nnn) {
	if (v != NULL) delete [] (v);
	nn=rhs.nn;
	nnn=rhs.nnn;
	v= nnn>0 ? new T[nnn] : NULL;
      }
      for (int i=0; i<nnn; i++) v[i]=rhs[i];
    }
  return *this;
}
  
template <class T>
inline T & Array<T>::operator[](const int i)	//subscripting
{
  if (i<0 || i>=nn) { throw("Array subscript out of bounds"); }
  return v[i];
}

template <class T>
inline const T & Array<T>::operator[](const int i) const	//subscripting
{
  if (i<0 || i>=nn) {  throw("Array subscript out of bounds"); }
  return v[i];
}

template <class T>
inline int Array<T>::size() const
{  return nn; }

template <class T>
inline int Array<T>::capacity() const
{  return nnn; }

template <class T>
void Array<T>::resize(int newn)
{
  if (newn != nnn ) {
    if (v != NULL) delete[] (v);
    nnn = nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
}

template <class T>
void Array<T>::reuse(int newn)
{
  if (newn>0 && newn <= nnn) nn = newn;
  else {
    delete[] (v);
    nnn = nn = 0;
    v = NULL;
  }
}

template <class T>
void Array<T>::assign(int newn, const T& a)
{
  if (newn > nnn) {
    if (v != NULL) delete[] (v);
    nnn = nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
  else nn=newn;
  for (int i=0;i<nn;i++) v[i] = a;
}

template <class T>
void Array<T>::assign(const T& a)
{ for (int i=0;i<nn;i++) v[i] = a; }

template <class T>
void Array<T>::assign(int newn, const T* a)
{
  if (newn > nnn) {
    if (v != NULL) delete[] (v);
    nnn = nn = newn;
    v = nn > 0 ? new T[nn] : NULL;
  }
  else nn=newn;
  for (int i=0;i<nn;i++) v[i] = a[i];
}

template <class T>
Array<T>::~Array()
{  if (v != NULL) delete[] (v); }

// end of Array definitions
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

// 2D matrix definations
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

template <class T>
class Mat {
 private:
  int nn;
  int mm;
  T **v;
 public:
  Mat();
  Mat(int n, int m);			// Zero-based Array
  Mat(int n, int m, const T &a);	// Initialize to constant
  Mat(int n, int m, const T *a);	// Initialize to Array
  Mat(const Mat &rhs);		// Copy constructor
  Mat & operator=(const Mat &rhs);	// Assignment
  typedef T value_type;                 // Make T available externally
  inline T* operator[](const int i);	// Subscripting: pointer to row i
  inline const T* operator[](const int i) const;
  inline int nrows() const;
  inline int ncols() const;
  void resize(int newn, int newm);      // Resize (contents not preserved)
  void assign(int newn, int newm, const T &a); // resize and assign a constant value
  ~Mat();
};

template <class T>
Mat<T>::Mat() : nn(0), mm(0), v(NULL) {}

template <class T>
Mat<T>::Mat(int n, int m) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1;i<n;i++) v[i] = v[i-1] + m;
}

template <class T>
Mat<T>::Mat(int n, int m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
}

template <class T>
Mat<T>::Mat(int n, int m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : NULL)
{
  int i,j,nel=m*n;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< n; i++) v[i] = v[i-1] + m;
  for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
}

template <class T>
Mat<T>::Mat(const Mat &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : NULL)
{
  int i,j,nel=mm*nn;
  if (v) v[0] = nel>0 ? new T[nel] : NULL;
  for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
}

template <class T>
Mat<T> & Mat<T>::operator=(const Mat<T> &rhs)
{
  // postcondition: normal assignment via copying has been performed;
  //		if matrix and rhs were different sizes, matrix
  //		has been resized to match the size of rhs
  if (this != &rhs) {
    int i,j,nel;
    if (nn != rhs.nn || mm != rhs.mm) {
      if (v != NULL) {
	delete[] (v[0]);
	delete[] (v);
      }
      nn=rhs.nn;
      mm=rhs.mm;
      v = nn>0 ? new T*[nn] : NULL;
      nel = mm*nn;
      if (v) v[0] = nel>0 ? new T[nel] : NULL;
      for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
    }
    for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
  }
  return *this;
}
  
template <class T>
inline T* Mat<T>::operator[](const int i)  //subscripting: pointer to row i
{
  if (i<0 || i>=nn) { throw("Mat subscript out of bounds"); }
  return v[i];
}

template <class T>
inline const T* Mat<T>::operator[](const int i) const
{
  if (i<0 || i>=nn) { throw("Mat subscript out of bounds"); }
  return v[i];
}

template <class T>
inline int Mat<T>::nrows() const { return nn; }

template <class T>
inline int Mat<T>::ncols() const { return mm; }

template <class T>
void Mat<T>::resize(int newn, int newm)
{
  int i,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
}

template <class T>
void Mat<T>::assign(int newn, int newm, const T& a)
{
  int i,j,nel;
  if (newn != nn || newm != mm) {
    if (v != NULL) {
      delete[] (v[0]);
      delete[] (v);
    }
    nn = newn;
    mm = newm;
    v = nn>0 ? new T*[nn] : NULL;
    nel = mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : NULL;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
  }
  for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
}

template <class T>
Mat<T>::~Mat()
{
  if (v != NULL) {
    delete[] (v[0]);
    delete[] (v);
  }
}

// 2D matrix definations
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

// 3D matrix definations
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

template <class T>
class Mat3D {
 private:
  int nn;
  int mm;
  int kk;
  T ***v;
 public:
  Mat3D();
  Mat3D(int n, int m, int k);
  inline T** operator[](const int i);	//subscripting: pointer to row i
  inline const T* const * operator[](const int i) const;
  inline int dim1() const;
  inline int dim2() const;
  inline int dim3() const;
  ~Mat3D();
};

template <class T>
Mat3D<T>::Mat3D(): nn(0), mm(0), kk(0), v(NULL) {}

template <class T>
Mat3D<T>::Mat3D(int n, int m, int k) : nn(n), mm(m), kk(k), v(new T**[n])
{
  int i,j;
  v[0] = new T*[n*m];
  v[0][0] = new T[n*m*k];
  for(j=1; j<m; j++) v[0][j] = v[0][j-1] + k;
  for(i=1; i<n; i++) {
    v[i] = v[i-1] + m;
    v[i][0] = v[i-1][0] + m*k;
    for(j=1; j<m; j++) v[i][j] = v[i][j-1] + k;
  }
}

template <class T>
inline T** Mat3D<T>::operator[](const int i) //subscripting: pointer to row i
{ return v[i]; }

template <class T>
inline const T* const * Mat3D<T>::operator[](const int i) const
{ return v[i]; }

template <class T>
inline int Mat3D<T>::dim1() const { return nn; }

template <class T>
inline int Mat3D<T>::dim2() const { return mm; }

template <class T>
inline int Mat3D<T>::dim3() const { return kk; }

template <class T>
Mat3D<T>::~Mat3D()
{
  if (v != NULL) {
    delete[] (v[0][0]);
    delete[] (v[0]);
    delete[] (v);
  }
}
// end of 3D matrix definations
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

/* convert number or anything to string */
template <class T>
inline std::string to_string (const T& t)
{
  std::stringstream ss;
  ss << t;
  return ss.str();
}



// http://www.cplusplus.com/forum/general/20554/
// http://stackoverflow.com/questions/1134388/stdendl-is-of-unknown-type-when-overloading-operator
/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
class dual_stream {
 public:
 dual_stream(std::ostream& os1, std::ostream& os2) : 
  os1(os1), os2(os2) {};
  
  template<class T>
    dual_stream& operator<<(const T& x) {
    os1 << x;
    os2 << x;
    return *this;
  }
  
  // function that takes a custom stream, and returns it
  typedef dual_stream& (*dual_streamManipulator)(dual_stream&);
  
  // take in a function with the custom signature
  dual_stream& operator<<(dual_streamManipulator manip)
    {
      // call the function, and return it's value
      return manip(*this);
    }
  
  // define the custom endl for this stream.
  // note how it matches the `dual_streamManipulator`
  // function signature
  static dual_stream& endl(dual_stream& stream)
  {
    // print a new line
    std::cout << std::endl;
    
    // do other stuff with the stream
    // std::cout, for example, will flush the stream
    stream << "Called dual_stream::endl!" << std::endl;
    
    return stream;
  }
  
  // this is the type of std::cout
  typedef std::basic_ostream<char, std::char_traits<char> > CoutType;
  
  // this is the function signature of std::endl
  typedef CoutType& (*StandardEndLine)(CoutType&);
  
  // define an operator<< to take in std::endl
  dual_stream& operator<<(StandardEndLine manip)
    {
      // call the function, but we cannot return it's value
      // manip(std::cout);
      manip(os1);
      manip(os2);
      
      return *this;
    }
  
 private:
  std::ostream& os1;
  std::ostream& os2;
};

/*ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
#endif
