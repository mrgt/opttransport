#ifndef OT_DENSITY_BASE_2_HPP
#define OT_DENSITY_BASE_2_HPP

#include <CGAL/Polygon_2.h>
#include <CGAL/intersection_2.h>
//#include <CGAL/Ray_2_Polygon_2_intersection.h>
#include <OT/random.hpp>
#include <OT/intersection_orourke.hpp>

namespace OT
{
  namespace details
  {
    template <class K>
    void output_polygon_ps(std::ostream &os,
			   const typename CGAL::Polygon_2<K> &p)
    {
      size_t N = p.size();
      if (N == 0)
	return;
      
      //std::cout << "is_convex?" << p.is_convex() << "\n";
      os << "newpath\n";
      os << p[N-1].x() << "\t" <<  p[N-1].y() << " moveto\n";
      for (size_t i = 0; i < p.size(); ++i)
	os << p[i].x() << "\t" <<  p[i].y() << " lineto\n";

      std::cerr << "closepath stroke\n";
    }

    template <class BBox>
    double compute_clamp_radius(const BBox &bb)
    {
      double x = std::sqrt(pow(bb.xmin(), 2.0) + pow(bb.ymin(), 2.0));
      double y = std::sqrt(pow(bb.xmax(), 2.0) + pow(bb.ymax(), 2.0));
      return 2.0 * std::max(x,y);
    }
  }

  template <class K>
  class Density_base_2
  {  
  public:
    typedef typename CGAL::Polygon_2<K> Polygon;
    typedef typename K::Segment_2 Segment;
    typedef typename K::Ray_2 Ray;
    typedef typename K::Point_2 Point;
    typedef typename K::Vector_2 Vector;
    typedef typename K::FT FT;
    typedef typename K::Iso_rectangle_2 BBox;

  protected:
    Polygon _bounding_poly;
    FT _r;

  public:
    Density_base_2(): _r(0.0) {}

    Density_base_2(const Polygon &bpoly) :
      _bounding_poly(bpoly),
      _r(details::compute_clamp_radius(bpoly.bbox()))
    {}

    const Polygon &
    bounding_poly() const
    {
      return _bounding_poly;
    }

    FT
    clamp_radius() const
    {
      return _r;
    }
  };

  template<class K>
    bool clip (const typename CGAL::Polygon_2<K> &bound,
	       const typename CGAL::Polygon_2<K> &poly, 
	       typename CGAL::Polygon_2<K> & clipped_poly)
  {
    return OT::convex_intersection (bound, poly, clipped_poly);
  }
#if 0
  template<class K>
    bool clip (const typename CGAL::Polygon_2<K> &bound,
	       const typename K::Segment_2 &segment,
	       typename K::Segment_2 &clipped_segment) 
  {
    CGAL::Object isect_object =
      CGAL::internal::intersection(segment, bound, K());
    return CGAL::assign(clipped_segment, isect_object);
  }
  
  template <class K>
    bool clip (const typename CGAL::Polygon_2<K> &bound,
	       const typename K::Ray_2 &ray,
	       typename K::Segment_2 &clipped_segment) 
  {
    CGAL::Object isect_object =
      CGAL::internal::intersection(ray, bound, K());
    return CGAL::assign(clipped_segment, isect_object);
  }
#endif

  template <class K, class Density>
  double
  mass (const Density &mu,
	const CGAL::Polygon_2<K> &upoly,
	bool clipped = false)
  {
    if (clipped)
      return mu.mass_(upoly);
    else 
      {
	typename CGAL::Polygon_2<K> poly;
	if (!clip(mu.bounding_poly(), upoly, poly))
	  return 0;
	return mu.mass_(poly);
      }
  }

  template <class Density, class K>
  double
  mass (const Density &mu,
	const typename CGAL::Segment_2<K> &useg,
	bool clipped = false)
  {
    if (clipped)
      mu.mass_(useg);
    else 
      {
	typename K::Segment_2 seg;
	if (!clip(mu.bounding_poly(), useg, seg))
	  return 0;
	return mu.mass_(seg);
      }
  }

  template <class K, class Density>
  double
  mass (const Density &mu,
	const typename CGAL::Ray_2<K> &uray)
  {    
    typename K::Segment_2 seg;
    if (!clip(mu.bounding_poly(), uray, seg))
      return 0;
    return mu.mass_(seg);
  }


  template <class K, class Density>
  bool
  centroid (const Density &mu,
	    const CGAL::Polygon_2<K> &upoly,
	    std::pair<typename K::Point_2, double> &p,
	    bool clipped = false)
  {
    if (clipped)
      return mu.centroid_(upoly, p);
    else 
      {
	typename CGAL::Polygon_2<K> poly;
	if (!clip(mu.bounding_poly(), upoly, poly))
	  return false;
	return mu.centroid_(poly, p);
      }
  }


  template <class K, class Density, class Function>
  double
  integrate (const Density &mu,
	     const Function &f,
	     const typename CGAL::Polygon_2<K> &upoly,
	     bool clipped = false) 
  {
    if (clipped)
      mu.integrate_(upoly, f);
    else 
      {
	typename CGAL::Polygon_2<K> poly;
	if (!clip(mu.bounding_poly(), upoly, poly))
	  return 0;
	return mu.integrate_(poly, f);
      }
  }

#if 0
  template <class K, class Density, class Function, class Linear>
  double
    integrate_linear_ (const Density &mu,
			const Function &f,
			const Linear &useg,
			bool clipped = false)
  {
    if (clipped)
      mu.integrate_(f, useg);
    else 
      {
	typename K::Segment_2 seg;
	if (!clip(mu.bounding_poly(), useg, seg))
	  return 0;
	return mu.integrate_(seg);
      }
  }

  template <class K, class Density, class Function>
  double
    integrate (const Density &mu,
	       const Function &f,
	       const typename K::Segment_2 &useg,
	       bool clipped = false)
  {
    return integrate_linear_(mu, f, useg, clipped);
  }

  template <class K, class Density, class Function>
  double
  integrate (const Density &mu,
	     const Function &f,
	     const typename K::Ray_2 &uray,
	     bool clipped = false)
  {
    return integrate_linear_(mu, f, uray, clipped);
  }
#endif

  template <class Density, class Engine>
  typename Density::Point
  random_point (const Density &mu, Engine &eng)
  {
    return mu.random_point_(eng);
  }


}


#endif
