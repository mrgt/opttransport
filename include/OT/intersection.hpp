#ifndef OT_CONVEX_INTERSECTION_HPP
#define OT_CONVEX_INTERSECTION_HPP

//#define OT_USE_CGAL_POLYGON_INTERSECTION
#ifdef OT_USE_CGAL_POLYGON_INTERSECTION

#include <CGAL/intersections.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/Boolean_set_operations_2.h>

namespace OT
{
  template <class K>
  bool convex_intersection (const typename CGAL::Polygon_2<K> &A,
			    const typename CGAL::Polygon_2<K> &B,
			    typename CGAL::Polygon_2<K> &out)
  {
    typedef CGAL::Polygon_2<K> Polygon_2;
    typedef CGAL::Polygon_with_holes_2<K> Polygon_with_holes_2;
    typedef std::list<Polygon_with_holes_2> Pwh_list_2;

    Pwh_list_2 intR;
    typename Pwh_list_2::const_iterator it;

    CGAL::intersection (A, B, std::back_inserter(intR));

    if (intR.size() != 1)
      return false;

    const Polygon_with_holes_2 &poly = *intR.begin();

    if (poly.is_unbounded())
      return false;

    out = poly.outer_boundary();
    return true;
  }
}

#else
#include <OT/intersection_orourke.hpp>
#endif


#endif

