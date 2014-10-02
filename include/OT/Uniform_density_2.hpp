#ifndef OT_UNIFORM_DENSITY_HPP
#define OT_UNIFORM_DENSITY_HPP

#include <CGAL/Bbox_2.h>
#include <OT/intersection.hpp>
#include <OT/Density_base_2.hpp>

namespace OT
{
  template <class K>
  class Uniform_density_2 : public Density_base_2<K>
  {  
  public:
    typedef typename OT::Density_base_2<K> Parent;
    typedef typename CGAL::Polygon_2<K> Polygon;
    typedef typename K::Segment_2 Segment;
    typedef typename K::Ray_2 Ray;
    typedef typename K::Point_2 Point;
    typedef typename K::Vector_2 Vector;
    typedef typename K::FT FT;
    typedef typename K::Iso_rectangle_2 BBox;
    typedef OT::Random_point_in_polygon<K> Random_point_in_polygon;

  private:
    double _scale;
    Random_point_in_polygon _rpp;

    double compute_scale()
    {
      double area = Parent::bounding_poly().area();
      OT_ASSERT(area > 0.0);
      return 1.0/area;
    }

  public:
    Uniform_density_2() : _scale(0.0) {}

    Uniform_density_2(const Polygon &bpoly) :
      Density_base_2<K>(bpoly),
      _scale(compute_scale()),
      _rpp (Parent::bounding_poly())
    {}
   
    inline bool 
    centroid_ (const Polygon &poly,
	       std::pair <Point, FT> &p) const
    {
      size_t N = poly.size();
      Vector vcentroid (0.0, 0.0);
      double weight = 0.0;

      for (size_t i = 2; i < N; ++i)
	{
	  // triangle 0, i-1, i
	  FT area = fabs(CGAL::area(poly[0], poly[i-1], poly[i]));
	  vcentroid = vcentroid + 
	    (area/3.0) * ( (poly[0] - CGAL::ORIGIN) + 
			   (poly[i-1] - CGAL::ORIGIN) + 
			   (poly[i] - CGAL::ORIGIN));
	  weight += area;
	}

      if (CGAL::is_zero(weight))
	return false;
      
      p = std::make_pair
	(Point ( CGAL::ORIGIN + ((1.0/weight) * vcentroid) ),
	 weight * _scale);
      return true;
    }
    
    double
    mass_ (const Polygon &poly) const
    {
      return poly.area() * _scale;
    }

    double 
    mass_ (const Segment &s) const
    {
      return sqrt(CGAL::squared_distance (s.source(), s.target())) * _scale;
    }

    template <class Function>
    double 
    integrate_ (const Polygon &poly, 
		const Function &f) const
    {
      double res = 0.0;
      for (size_t i = 2; i < poly.size(); ++i)
	res += f(poly[0], poly[i-1], poly[i]);
      return res * _scale;
    }

    template <class Function>
    double 
    integrate_ (const Segment &segment, 
		const Function &f) const
    {
      return f(segment.source(), segment.target());
    }

    template <class Engine>
    Point
    random_point_(Engine &eng) const
    {
      return _rpp(eng);
    }
  };

}

#endif
