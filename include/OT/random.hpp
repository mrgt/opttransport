#ifndef OT_RANDOM_HPP
#define OT_RANDOM_HPP

#include <CGAL/Polygon_2.h>
#include <CGAL/Origin.h>
#include <CGAL/Kernel/global_functions_2.h>

#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>

#include <OT/misc.hpp>

namespace OT
{ 
  template <class K>
  class Random_point_in_triangle
  {
  public:
    typedef typename K::Point_2 Point;
    typedef typename K::Vector_2 Vector;
    Vector _a, _b, _c;

  public:
    Random_point_in_triangle (const Point &a,
			      const Point &b,
			      const Point &c) : 
      _a(a - CGAL::ORIGIN),
      _b(b - CGAL::ORIGIN),
      _c(c - CGAL::ORIGIN)
    {}

    template <class Engine>
    Point operator () (Engine &eng) const
    {
      boost::uniform_real<> u01;
      double ta(u01(eng)), tb(u01(eng));
      
      if (ta+tb > 1.0)
	{
	  ta = 1.0 - ta;
	  tb = 1.0 - tb;
	}
      double tc = 1.0 - ta - tb;
      
      return CGAL::ORIGIN + (ta * _a + tb * _b + tc * _c);
    }
  };

  template <class K>
  class Random_point_in_polygon
  {
  public:
    typedef typename K::Point_2 Point;
    typedef typename CGAL::Polygon_2<K> Polygon;
    typedef typename OT::Random_point_in_triangle<K> RP_in_triangle;
    
    std::vector<RP_in_triangle> _rptri;
    std::vector<double> _cum;

  public:
    Random_point_in_polygon ()
    {}

    Random_point_in_polygon (const Polygon &poly)
    {
      double cum = 0.0;

      for (size_t i = 2; i < poly.size(); ++i)
	{
	  const Point &a = poly[0], &b = poly[i-1], &c = poly[i];
	  _rptri.push_back(RP_in_triangle(a,b,c));
	  double area = fabs(CGAL::area(a,b,c));
	  cum += area;
	  _cum.push_back(cum);
	}

      if (cum == 0.0)
	return;

      // normalize _cum between 0 and 1
      for (size_t i = 0; i < _cum.size(); ++i)
	  _cum[i] /= cum;
    }

    template <class Engine>
    Point operator () (Engine &eng) const
    {
      OT_ASSERT(_cum.size() > 0);

      boost::uniform_real<> u01(0,1);
      double t = u01(eng);

      for (size_t i = 0; i < _cum.size(); ++i)
	{
	  if (t <= _cum[i])
	    return (_rptri[i])(eng);
	}
      return Point(0.0,0.0);
    }
  };
}

#endif
