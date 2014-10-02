#ifndef CGAL_RAY_2_POLYGON_2_INTERSECTION_H
#define CGAL_RAY_2_POLYGON_2_INTERSECTION_H

#include <CGAL/Segment_2.h>
#include <CGAL/Ray_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/number_utils.h>
#include <CGAL/Straight_2.h>



#include <CGAL/Object.h>

CGAL_BEGIN_NAMESPACE

namespace internal 
{
  template <class K, class StraightType>
  class Straight_2_Polygon_2_pair
  {
    typedef typename CGAL::Polygon_2<K> Polygon;
  public:
    enum Intersection_results {NO_INTERSECTION, POINT, SEGMENT};
    Straight_2_Polygon_2_pair(StraightType const  *str,
			      Polygon const *poly)
      : _str(str), _poly(poly), _known(false) {}

    Intersection_results intersection_type() const
    {
      typedef typename K::Line_2  Line_2;
      if (_known)
	return _result;
      _known = true;
      Straight_2_<K> straight(*_str);
      
      // FIXME: we assume that the polygon is convex and
      // oriented counterclockwise
      size_t N = _poly->size();
      for (size_t i = 0; i < N; ++i)
	{
	  straight.cut_right_off(Line_2(_poly->vertex(i),
					_poly->vertex((i+1)%N)));
	}
      switch (straight.current_state()) {
      case Straight_2_<K>::EMPTY:
        _result = NO_INTERSECTION;
        return _result;
      case Straight_2_<K>::POINT: {
        straight.current(_intersection_point);
        _result = POINT;
        return _result;
      }
      case Straight_2_<K>::SEGMENT: {
        typename K::Segment_2 seg;
        straight.current(seg);
        _intersection_point = seg.source();
        _other_point = seg.target();
        _result = SEGMENT;
        return _result;
      }
      default:  // should not happen.
        CGAL_kernel_assertion_msg(false, "Internal CGAL error.");
        _result = NO_INTERSECTION;
        return _result;
      }
    }

    typename K::Point_2    intersection_point() const
    {
      if (!_known)
	intersection_type();
      CGAL_kernel_assertion(_result == POINT);
      return _intersection_point;
    }

    typename K::Segment_2  intersection_segment() const
    {
      typedef typename K::Segment_2 Segment_2;
      if (!_known)
	intersection_type();
      CGAL_kernel_assertion(_result == SEGMENT);
      return Segment_2(_intersection_point, _other_point);
    }

  protected:
    StraightType const * _str;
    Polygon const      * _poly;
    mutable bool                    _known;
    mutable Intersection_results     _result;
    mutable typename K::Point_2         _intersection_point;
    mutable typename K::Point_2         _other_point;
  };

  template <class K>
  class Ray_2_Polygon_2_pair : 
    public Straight_2_Polygon_2_pair<K, typename K::Ray_2>
  {
    typedef typename
    internal::Straight_2_Polygon_2_pair <K, typename K::Ray_2> Parent;

  public:
    Ray_2_Polygon_2_pair(typename K::Ray_2 const *ray,
			 typename CGAL::Polygon_2<K> const *poly) : 
      Parent (ray, poly)
    {}
  };

  template <class K>
  class Segment_2_Polygon_2_pair : 
    public Straight_2_Polygon_2_pair<K, typename K::Segment_2> 
  {
    typedef typename
    internal::Straight_2_Polygon_2_pair <K, typename K::Segment_2> Parent;

  public:
    Segment_2_Polygon_2_pair(typename K::Segment_2 const *seg,
			     typename CGAL::Polygon_2<K> const *poly):
      Parent (seg, poly)
    {}
  };  
  
  template <class K>
  Object
  intersection(const typename CGAL::Ray_2<K> &ray, 
	       const typename CGAL::Polygon_2<K> &tr,
	       const K&)
  {
    typedef Ray_2_Polygon_2_pair<K> is_t;
    is_t ispair(&ray, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
      return Object();
    case is_t::POINT:
      return make_object(ispair.intersection_segment());
    case is_t::SEGMENT:
      return make_object(ispair.intersection_segment());
    }
  }

  template <class K>
  Object
  intersection(const typename CGAL::Segment_2<K> &seg, 
	       const typename CGAL::Polygon_2<K> &tr,
	       const K&)
  {
    typedef Segment_2_Polygon_2_pair<K> is_t;
    is_t ispair(&seg, &tr);
    switch (ispair.intersection_type()) {
    case is_t::NO_INTERSECTION:
    default:
      return Object();
    case is_t::POINT:
      return make_object(ispair.intersection_segment());
    case is_t::SEGMENT:
      return make_object(ispair.intersection_segment());
    }
  }

  template <class K>
  Object
  intersection(const typename CGAL::Polygon_2<K> &tr,
	       const typename CGAL::Ray_2<K> &ray, 
	       const K& k)
  {
    return internal::intersection(ray, tr, k);
  }

  template <class K>
  Object
  intersection(const typename CGAL::Polygon_2<K> &tr,
	       const typename CGAL::Segment_2<K> &seg, 
	       const K& k)
  {
    return internal::intersection(seg, tr, k);
  }


  template <class K>
  inline bool do_intersect(const typename K::Ray_2 &p1,
			   const typename CGAL::Polygon_2<K> &p2,
			   const K&)
  {
    typedef Ray_2_Polygon_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
  }

  template <class K>
  inline bool do_intersect(const typename K::Segment_2 &p1,
			   const typename CGAL::Polygon_2<K> &p2,
			   const K&)
  {
    typedef Segment_2_Polygon_2_pair<K> pair_t;
    pair_t pair(&p1, &p2);
    return pair.intersection_type() != pair_t::NO_INTERSECTION;
  }

  template <class K>
  inline bool do_intersect(const typename CGAL::Polygon_2<K> &p1,
			   const typename K::Ray_2 &p2,
			   const K&)
  {
    return internal::do_intersect(p2,p1);
  }

  template <class K>
  inline bool do_intersect(const typename CGAL::Polygon_2<K> &p2,
			   const typename K::Segment_2 &p1,
			   const K&)
  {
    return internal::do_intersect(p2, p1);
  }
} // namespace internal


template <class K>
inline bool do_intersect(const CGAL::Polygon_2<K> &p,
			 const Ray_2<K> &ray)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(ray, p);
}


template <class K>
inline bool do_intersect(const CGAL::Polygon_2<K> &p,
			 const Segment_2<K> &seg)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(seg, p);
}

template <class K>
inline bool do_intersect(const Ray_2<K> &ray,
			 const CGAL::Polygon_2<K> &p)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(ray, p);
}

template <class K>
inline bool do_intersect(const Segment_2<K> &seg,
			 const CGAL::Polygon_2<K> &p)
{
  typedef typename K::Do_intersect_2 Do_intersect;
  return Do_intersect()(seg, p);
}

template <class K>
inline Object
intersection(const Ray_2<K> &ray, const CGAL::Polygon_2<K> &p)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(ray, p);
}

template <class K>
inline Object
intersection(const CGAL::Polygon_2<K> &p, const Ray_2<K> &ray)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(ray, p);
}


template <class K>
inline Object
intersection(const Segment_2<K> &seg, const CGAL::Polygon_2<K> &p)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(seg, p);
}

template <class K>
inline Object
intersection(const CGAL::Polygon_2<K> &p, const Segment_2<K> &seg)
{
  typedef typename K::Intersect_2 Intersect;
  return Intersect()(seg, p);
}


CGAL_END_NAMESPACE

#endif
