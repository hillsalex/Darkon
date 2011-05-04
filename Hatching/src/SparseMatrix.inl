/* ---------------------------------------------- 
   file: SparseMatrix.inl
   auth: Travis Fischer
   acct: tfischer
   date: Spring 2008

   Inline function definitions for a sparse matrix
   ---------------------------------------------- */

#ifndef MATRIX_SPARSE_INL_
#define MATRIX_SPARSE_INL_

#include <string.h>
#include <assert.h>
// Constructors
// ------------

// Constructs an m by n sparse matrix
inline SparseMatrix::SparseMatrix(int m, int n)
   : m_elements(), m_m(m), m_n(n)
{ }

inline SparseMatrix::SparseMatrix(int m, int n, const SparseArray &elements) 
   : m_elements(elements), m_m(m), m_n(n)
{ }

// Copy Constructor
inline SparseMatrix::SparseMatrix(const SparseMatrix &m) {
   m_elements = m.m_elements;
   m_m = m.m_m;
   m_n = m.m_n;
}


// Accessor Operators
// ------------------

// Returns the number of rows in this SparseMatrix
inline int SparseMatrix::getM() const {
   return m_m;
}

// Returns the number of cols in this SparseMatrix
inline int SparseMatrix::getN() const {
   return m_n;
}

inline SparseArray &SparseMatrix::getData() {
   return m_elements;
}

inline const SparseArray &SparseMatrix::getData() const {
   return m_elements;
}


// Equality Operators
// ------------------
inline       bool       SparseMatrix::operator==(const SparseMatrix &m) const {
   return (m_m == m.m_m && m_n == m.m_n && m_elements == m.m_elements);
}

inline       bool       SparseMatrix::operator!=(const SparseMatrix &m) const {
   return !((*this) == m);
}


// Mutator Operators
// -----------------
inline       SparseMatrix &SparseMatrix::operator =(const SparseMatrix &m) {
   m_m = m.m_m;
   m_n = m.m_n;
   m_elements = m.m_elements;
   
   return *this;
}

inline       SparseMatrix &SparseMatrix::operator*=(const SparseMatrix &m) {
   return (*this = (*this) * m);
}

inline       SparseMatrix &SparseMatrix::operator+=(const SparseMatrix &m) {
   return (*this = (*this) + m);
}

inline       SparseMatrix &SparseMatrix::operator-=(const SparseMatrix &m) {
   return (*this = (*this) - m);
}


// Arithmetic Operators
// -----------------
inline       SparseMatrix  SparseMatrix::operator* (const SparseMatrix &rhs) const {

   SparseColArrayConstIterator j, p;
   SparseArrayConstIterator i;
   
   SparseArray a;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i) {
      int row = i->first;
      SparseColArray &aRow = a[row];
      
      for(j = i->second.begin(); j != i->second.end(); ++j) {
         int col = j->first;
         
         if (rhs.m_elements.count(col) > 0) {
            const SparseColArray &o = rhs.m_elements.find(col)->second;
            
            for(p = o.begin(); p != o.end(); ++p)
               aRow[p->first] += p->second * j->second;
         }
      }
   }
  
   return SparseMatrix(m_m, rhs.m_n, a);
}

inline       SparseMatrix  SparseMatrix::operator+ (const SparseMatrix &rhs) const {

   SparseColArrayConstIterator j, p;
   SparseArrayConstIterator i;
   
   SparseArray a;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i) {
      int row = i->first;
      SparseColArray &aRow = a[row];
      
      for(j = i->second.begin(); j != i->second.end(); ++j) {
         int col = j->first;
         
         aRow[col] = j->second;
      }
   }
   
   for(i = rhs.m_elements.begin(); i != rhs.m_elements.end(); ++i) {
      int row = i->first;
      SparseColArray &aRow = a[row];
      
      for(j = i->second.begin(); j != i->second.end(); ++j) {
         int col = j->first;
         
         aRow[col] += j->second;
      }
   }
  
   return SparseMatrix(m_m, m_n, a);
}

inline       SparseMatrix  SparseMatrix::operator- (const SparseMatrix &rhs) const {

   SparseColArrayConstIterator j, p;
   SparseArrayConstIterator i;
   
   SparseArray a;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i) {
      int row = i->first;
      SparseColArray &aRow = a[row];
      
      for(j = i->second.begin(); j != i->second.end(); ++j) {
         int col = j->first;
         
         aRow[col] = j->second;
      }
   }
   
   for(i = rhs.m_elements.begin(); i != rhs.m_elements.end(); ++i) {
      int row = i->first;
      SparseColArray &aRow = a[row];
      
      for(j = i->second.begin(); j != i->second.end(); ++j) {
         int col = j->first;
         
         aRow[col] -= j->second;
      }
   }
  
   return SparseMatrix(m_m, m_n, a);
}

inline       SparseMatrix  SparseMatrix::operator* (const double scale) const {
   SparseArray a(m_elements);
   SparseColArrayIterator j;
   SparseArrayIterator i;
   
   for(i = a.begin(); i != a.end(); ++i)
      for(j = i->second.begin(); j != i->second.end(); ++j)
         j->second *= scale;

   return SparseMatrix(m_m, m_n, a);
}

inline       SparseMatrix  SparseMatrix::operator/ (const double scale) const {
   SparseArray a(m_elements);
   SparseColArrayIterator j;
   SparseArrayIterator i;
   

   
   for(i = a.begin(); i != a.end(); ++i)
      for(j = i->second.begin(); j != i->second.end(); ++j)
         j->second *= scale;

   return SparseMatrix(m_m, m_n, a);
}


// More Complex Functionality
// --------------------------

// Returns whether or not this SparseMatrix is square
inline bool SparseMatrix::isSquare() const {
   return (m_m == m_n);
}

// Returns whether or not this SparseMatrix is diagonal
inline bool SparseMatrix::isDiagonal() const {
   SparseColArrayConstIterator j;
   SparseArrayConstIterator i;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i)
      for(j = i->second.begin(); j != i->second.end(); ++j)
         if (i->first != j->second)
            return false;

   return true;
}


// Returns the number of elements currently in this SparseMatrix
inline int SparseMatrix::getSize() const {
   SparseColArrayConstIterator j;
   SparseArrayConstIterator i;
   int size = 0;
   
   //cerr << "SparseMatrix::getSize()" << endl;
   for(i = m_elements.begin(); i != m_elements.end(); ++i)
      size += i->second.size();
   
   return size;
}

// Returns the full size in bytes of this SparseMatrix if it were to be 
// expanded in non-sparse form
inline int SparseMatrix::fullSize() const {
   return m_m * m_n * sizeof(double);
}

// Stores a (row-major) full version of this SparseMatrix in paramater out, 
// which should be a preallocated buffer of fullSize() bytes
inline void SparseMatrix::toFull(double *out) const {
   memset(out, 0, fullSize());
   
   SparseColArrayConstIterator j;
   SparseArrayConstIterator i;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i)
      for(j = i->second.begin(); j != i->second.end(); ++j)
         out[i->first * m_n + j->first] = j->second;
}

// Adds the value to this SparseMatrix
inline SparseMatrix &SparseMatrix::addValue(int row, int col, double value) {

   m_elements[row][col] = value;
   return *this;
}

inline SparseMatrix &SparseMatrix::incValue(int row, int col, double value) {
    if (hasValue(row,col))
    {
        m_elements[row][col] = value + getValue(row,col);
    }
    else setValue(row,col,value);
    return *this;
}

// Adds the value to this SparseMatrix
inline SparseMatrix &SparseMatrix::setValue(int row, int col, double value) {
   return addValue(row, col, value);
}

// Returns whether or not this SparseMatrix contains a value for the given cell
inline bool   SparseMatrix::hasValue(int row, int col) {

   return m_elements.count(row) > 0 && m_elements[row].count(col) > 0;
}

// Returns the value stored at (row, col) or zero if none exists
inline double SparseMatrix::getValue(int row, int col) {
    assert(row<getM() && col < getN());
   if (!hasValue(row, col))
      return 0;
   
   return m_elements[row][col];
}

inline void SparseMatrix::setIdentityColumn(int col) {
   SparseArrayIterator i;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i) {
      if (i->second.count(col) > 0)
         i->second[col] = 0;
   }
   
   setValue(col, col, 1.0);
}

// Cleans up a matrix (0's out entries that are less than epsilon)
inline void SparseMatrix::cleanup() {
   SparseColArrayIterator j;
   SparseArrayIterator i;
   
   for(i = m_elements.begin(); i != m_elements.end(); ++i) {
      for(j = i->second.begin(); j != i->second.end(); ++j) {
         if (j->second == 0)
            j->second = 0;
      }
   }
}

inline void SparseMatrix::clearRow(int row) {

   if (m_elements.find(row) != m_elements.end())
      m_elements.erase(m_elements.find(row));
}

// Extra operators where SparseMatrix is on right-hand side
// -----------------------------------------------------
inline       SparseMatrix  operator* (const double scale, const SparseMatrix &rhs) {
   return rhs * scale;
}

#endif // MATRIX_SPARSE_INL_

