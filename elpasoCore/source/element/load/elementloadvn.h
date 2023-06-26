
#ifndef INFAM_ELEMENT_LOAD_VN
#define INFAM_ELEMENT_LOAD_VN

#include "elementload.h"


/**
 * specify the type of the boundary condition
 * normal to fluid's face
 */
enum eTypeLoad
{
  undefined,
  flux,
  vn
};


/**
 * @brief normal velocity vn or flux q acting normal to surface of fluid element
 * @author Dirk Clasen
 * @date 16.08.2005
 */
class cElementLoadVn :
  public cElementLoad,
  private tCounter<cElementLoadVn>
{
private:
  eTypeLoad m_TypeLoad; ///< specifies if m_Value is flux or normal velocity
  PetscReal m_Value;    ///< value for \f$ v_n = \frac{\partial p}{\partial n} \f$ or flux q

public :
  cElementLoadVn(const eTypeLoad &MyType = undefined);
  cElementLoadVn(const cElementLoadVn &other);
  ~cElementLoadVn();

  //! number of instances of this object
  static size_t howMany(void) { return tCounter<cElementLoadVn>::howMany(); }

  //! get the type of the load - either normal velocity or flux
  eTypeLoad getType(void) const { return m_TypeLoad; }

  //! set a new type to this element load
  void setType(const eTypeLoad &NewType) { m_TypeLoad = NewType; }

  //! return the value for the normal velocity or flux
  inline PetscReal getValue(void) const { return m_Value; }

  //! write this object to a stream
  //! @param os outputstream
  //! @return modified outputstream
  std::ostream& write(std::ostream &os) const;

  //! read data to this object
  //! @param is inputstream
  //! @return modified inputstream
  std::istream& read(std::istream &is);

  //! write this object to XML stream
  //! @param os outputstream
  //! @return der modified outputstream
  std::ostream& writeXml(std::ostream &os) const;
};


//! overloaded outputoperator
inline std::ostream& operator<<(std::ostream &os, const cElementLoadVn &other)
{
  return other.write(os);
}


//! overloaded inputoperator
inline std::istream& operator>>(std::istream &is, cElementLoadVn &other)
{
  return other.read(is);
}

#endif
