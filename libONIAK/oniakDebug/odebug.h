#ifndef __ONIAK_ODEBUG_H__
#define __ONIAK_ODEBUG_H__

#include <vector>

namespace ONIAK
{

template<typename T>
void DEBUG_VECTOR_EQUAL(const std::vector<T>& a, const std::vector<T>& b)
{
#ifdef ONIAK_DEBUG
  assert(a.size() == b.size());
  for (size_t i = 0; i < a.size(); i++)
  {
    assert(a[i] == b[i]);
  }
#endif
}

}

#endif