
#include "elementloadvn.h"


cElementLoadVn::cElementLoadVn(const eTypeLoad &MyType)
{
  m_Value = 0.;
  setType( MyType );
}


cElementLoadVn::cElementLoadVn(const cElementLoadVn &other) :
  cElementLoad(other)
{
  m_Value = other.getValue();
  setType( other.getType( ) );
}


cElementLoadVn::~cElementLoadVn()
{
  // empty
}


std::ostream& cElementLoadVn::write(std::ostream &os) const
{
  if (m_TypeLoad == flux)
  {
    os << "flux boundary" << std::endl;
    os << "  Id = " << getId() << std::endl;
    os << "  q  = " << getValue() << std::endl;
  }
  else if (m_TypeLoad == vn)
  {
    os << "normal velocity" << std::endl;
    os << "  Id = " << getId() << std::endl;
    os << "  vn  = " << getValue() << std::endl;
  }
  else
  {
    throw cException("no type specified for cElementLoadVn", __FILE__, __LINE__);
  }

  return os;
}


std::istream& cElementLoadVn::read(std::istream &is)
{
  cId::read(is);
  is >> m_Value;
  return is;
}


std::ostream& cElementLoadVn::writeXml(std::ostream &os) const
{
  if (m_TypeLoad == flux)
  {
    os << "<ElemLoad Type=\"flux\">";
    os << "<Id>" << getId() << "</Id>";
    os << "<q>" << getValue() << "</q>";
    os << "</ElemLoad>";
  }
  else if (m_TypeLoad == vn)
  {
    os << "<ElemLoad Type=\"vn\">";
    os << "<Id>" << getId() << "</Id>";
    os << "<vn>" << getValue() << "</vn>";
    os << "</ElemLoad>";
  }
  else
  {
    throw cException("no type specified for cElementLoadVn", __FILE__, __LINE__);
  }

  return os;
}
