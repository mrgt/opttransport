#ifndef OT_SERIALIZATION_HPP
#define OT_SERIALIZATION_HPP

#include <boost/archive/tmpdir.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/assume_abstract.hpp>

namespace boost
{
  namespace serialization {
    
    template<class Archive, class K>
    void save(Archive & ar, const typename CGAL::Point_2<K> &p,
	      const unsigned int version)
    {
      ar & p.x();
      ar & p.y();      
    }

    template<class Archive, class K>
    void load(Archive & ar, typename CGAL::Point_2<K> &p,
	      const unsigned int version)
    {
      double x, y;
      ar & x;
      ar & y;
      p = CGAL::Point_2<K> (x,y);
    }

  } // namespace serialization
}

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

BOOST_SERIALIZATION_SPLIT_FREE
  (CGAL::Point_2<CGAL::Exact_predicates_inexact_constructions_kernel>);

#endif
