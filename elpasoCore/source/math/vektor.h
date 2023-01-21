/* Copyright (c) 2023. Authors listed in AUTHORS.md

 * This file is part of elPaSo-Core.

 * elPaSo-Core is free software: you can redistribute it and/or modify it
 * under the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or (at your option)
 * any later version.

 * elPaSo-Core is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.

 * You should have received a copy of the GNU Lesser General Public License along
 * with elPaSo-Core (COPYING.txt and COPYING.LESSER.txt). If not, see
 * <https://www.gnu.org/licenses/>. 
 */

#ifndef INFAM_VEKTOR_H
#define INFAM_VEKTOR_H

#include <cmath>
#include <iomanip>
#include <iostream>


namespace infam {
//! @brief vector class
//! @author Dirk Clasen
//! @date 22.06.2005
template <class TFLOAT, class TUINT>
class tVector {
 private:
  TUINT m_Size;    //! number of entries
  TFLOAT *m_Data;  //! vector entries

  void merge_sort(tVector<TFLOAT, TUINT> &temp, int left, int right);
  void merge(tVector<TFLOAT, TUINT> &temp, int left, int mid, int right);

 public:
  tVector();
  tVector(TUINT Size);
  tVector(const tVector<TFLOAT, TUINT> &other);
  ~tVector();

  //! return the size of this vector
  inline TUINT size(void) const { return m_Size; }

  //! set new size of vector
  void resize(TUINT Size);

  //! access a single value of the vector (read/write)
  inline TFLOAT &operator[](TUINT pos) { return m_Data[pos]; };

  //! access a single value of the vector (read only)
  inline TFLOAT operator[](TUINT pos) const { return m_Data[pos]; };

  //! v_i /= value
  tVector<TFLOAT, TUINT> &operator/=(const TFLOAT &value);

  //! allocation operator " = " for class
  tVector<TFLOAT, TUINT> &operator=(const tVector<TFLOAT, TUINT> &other);

  //! assign a value to each cell of the vector
  void setValue(const TFLOAT &value);

  //! |v|^2
  TFLOAT abs2() const;

  //! abs of the vector
  inline TFLOAT abs() const { return std::sqrt(abs2()); }

  //! dot product this . v
  TFLOAT dot(const tVector<TFLOAT, TUINT> &v) const;

  //! vector will be sorted
  void sort();

  //! returns minimum value
  TFLOAT getMin() const;

  //! returns maximum value
  TFLOAT getMax() const;

  //! write this vector to a stream
  std::ostream &write(std::ostream &os) const;
};

//! @note unit-tested
template <class TFLOAT, class TUINT>
tVector<TFLOAT, TUINT>::tVector() {
  m_Size = 0;
  m_Data = new TFLOAT[m_Size];
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
tVector<TFLOAT, TUINT>::tVector(TUINT Size) {
  m_Size = Size;
  m_Data = new TFLOAT[m_Size];
  for (TUINT k = 0; k < Size; k++) m_Data[k] = 0.;
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
tVector<TFLOAT, TUINT>::tVector(const tVector<TFLOAT, TUINT> &other) {
  m_Size = other.size();
  m_Data = new TFLOAT[m_Size];
  for (TUINT k = 0; k < m_Size; k++) m_Data[k] = other[k];
}

template <class TFLOAT, class TUINT>
tVector<TFLOAT, TUINT>::~tVector() {
  delete[] m_Data;
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
inline void tVector<TFLOAT, TUINT>::resize(TUINT Size) {
  if (m_Size != Size) {
    m_Size = Size;
    delete[] m_Data;
    m_Data = new TFLOAT[m_Size];
    for (TUINT k = 0; k < Size; k++) m_Data[k] = 0.;
  }
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
inline void tVector<TFLOAT, TUINT>::setValue(const TFLOAT &value) {
  for (TUINT k = 0; k < m_Size; k++) m_Data[k] = value;
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
tVector<TFLOAT, TUINT> &tVector<TFLOAT, TUINT>::operator/=(
    const TFLOAT &value) {
  for (TUINT k = 0; k < m_Size; k++) m_Data[k] /= value;

  return *this;
}

//! @note unit-tested
//! allocation operator " = " for class
template <class TFLOAT, class TUINT>
tVector<TFLOAT, TUINT> &tVector<TFLOAT, TUINT>::operator=(
    const tVector<TFLOAT, TUINT> &other) {
  if (m_Size != other.size()) {
    m_Size = other.size();
    delete[] m_Data;
    m_Data = new TFLOAT[m_Size];
  }

  for (TUINT k = 0; k < m_Size; k++) m_Data[k] = other[k];

  return *this;
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
TFLOAT tVector<TFLOAT, TUINT>::abs2() const {
  TFLOAT res = TFLOAT(0.0);

  for (TUINT i = 0; i < size(); i++) res += m_Data[i] * m_Data[i];

  return res;
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
inline TFLOAT tVector<TFLOAT, TUINT>::dot(
    const tVector<TFLOAT, TUINT> &v) const {
  TFLOAT res = TFLOAT(0.0);
  if (m_Size == v.size()) {
    for (TUINT i = 0; i < size(); i++) res += m_Data[i] * v[i];
  } else {
    throw cException("illegal call dot : improper size of vectors!", __FILE__,
                     __LINE__);
  }

  return res;
}

// ---------------------------------------------------------------------------
//   crossproduct  c = a x b
// ---------------------------------------------------------------------------
//! @note unit-tested
template <class TFLOAT, class TUINT>
inline void cross_product(const tVector<TFLOAT, TUINT> &a,
                          const tVector<TFLOAT, TUINT> &b,
                          tVector<TFLOAT, TUINT> &c) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = -a[0] * b[2] + a[2] * b[0];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

// ---------------------------------------------------------------------------
//   scale all cells of the vector v by c
//   @note unit-tested
// ---------------------------------------------------------------------------
template <class TFLOAT1, class TFLOAT2, class TUINT>
inline void scale(tVector<TFLOAT1, TUINT> &V, const TFLOAT2 &c) {
  for (TUINT z = 0; z < V.size(); z++) V[z] *= c;
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
void tVector<TFLOAT, TUINT>::sort() {
  tVector<TFLOAT, TUINT> temp(size());
  merge_sort(temp, 0, size() - 1);
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
void tVector<TFLOAT, TUINT>::merge_sort(tVector<TFLOAT, TUINT> &temp, int left,
                                        int right) {
  int mid = 0;
  if (right > left) {
    mid = (right + left) / 2;
    merge_sort(temp, left, mid);
    merge_sort(temp, mid + 1, right);
    merge(temp, left, mid + 1, right);
  }
}

//! @note unit-tested
template <class TFLOAT, class TUINT>
void tVector<TFLOAT, TUINT>::merge(tVector<TFLOAT, TUINT> &temp, int left,
                                   int mid, int right) {
  int i, left_end, num_elements, tmp_pos;

  left_end = mid - 1;
  tmp_pos = left;
  num_elements = right - left + 1;

  while ((left <= left_end) && (mid <= right)) {
    if (m_Data[left] <= m_Data[mid]) {
      temp[tmp_pos] = m_Data[left];
      tmp_pos = tmp_pos + 1;
      left = left + 1;
    } else {
      temp[tmp_pos] = m_Data[mid];
      tmp_pos = tmp_pos + 1;
      mid = mid + 1;
    }
  }
  while (left <= left_end) {
    temp[tmp_pos] = m_Data[left];
    left = left + 1;
    tmp_pos = tmp_pos + 1;
  }
  while (mid <= right) {
    temp[tmp_pos] = m_Data[mid];
    mid = mid + 1;
    tmp_pos = tmp_pos + 1;
  }
  for (i = 0; i < num_elements; i++) {
    m_Data[right] = temp[right];
    right = right - 1;
  }
}

template <class TFLOAT, class TUINT>
inline TFLOAT tVector<TFLOAT, TUINT>::getMin() const {
  TFLOAT res = TFLOAT(m_Data[0]);
  for (int i = 1; i < m_Size; ++i)
    if (m_Data[i] < res) res = m_Data[i];
  return res;
}

template <class TFLOAT, class TUINT>
inline TFLOAT tVector<TFLOAT, TUINT>::getMax() const {
  TFLOAT res = TFLOAT(m_Data[0]);
  for (int i = 1; i < m_Size; ++i)
    if (m_Data[i] > res) res = m_Data[i];
  return res;
}

/*BEGIN_NO_COVERAGE*/
template <class TFLOAT, class TUINT>
inline std::ostream &tVector<TFLOAT, TUINT>::write(std::ostream &os) const {
  int w = os.width();
  int p = os.precision();

  os.setf(std::ios::scientific);
  for (TUINT i = 0; i < size(); i++)
    os << std::setw(15) << m_Data[i] << std::endl;
  os.unsetf(std::ios::scientific);

  os.width(w);
  os.precision(p);
  return os;
}

// ---------------------------------------------------------------------------
//   overloaded outputoperator <<
// ---------------------------------------------------------------------------
template <class TFLOAT, class TUINT>
inline std::ostream &operator<<(std::ostream &os,
                                const tVector<TFLOAT, TUINT> &v) {
  return v.write(os);
}
/*END_NO_COVERAGE*/

}  // namespace infam

#endif
