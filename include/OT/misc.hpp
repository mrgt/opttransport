#ifndef OT_MISC_HPP
#define OT_MISC_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iterator>

namespace OT
{
  namespace ublas = boost::numeric::ublas;
  typedef ublas::vector<double> uvector;
  typedef ublas::matrix<double, ublas::column_major> umatrix;

  struct null_output_iterator : 
    std::iterator< std::output_iterator_tag,
                   null_output_iterator > {
    /* no-op assignment */
    template<typename T>
    void operator=(T const&) { }

    null_output_iterator & operator++() { 
        return *this; 
    }

    null_output_iterator operator++(int) { 
        null_output_iterator it(*this); 
        ++*this;
        return it;
    }

    null_output_iterator & operator*() { return *this; }
  };

  bool is_null_output_iterator(const null_output_iterator &oi)
  {
    return true;
  }

  template <class OutputIterator>
  bool is_null_output_iterator(const OutputIterator &oi)
  {
    return false;
  }
}

#define OT_DEBUG_SHOW(x) std::cerr << "[" << __FUNCTION__ << "] " <<	\
      #x << " = " << x << "\n";

#undef NDEBUG
#include <assert.h>
#define OT_ASSERT(x) assert(x)

#endif
