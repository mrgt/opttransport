#ifndef OT_INTEGRATE_HPP
#define OT_INTEGRATE_HPP

#include <CGAL/Cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Simple_cartesian.h>
#include <utility>

#include <OT/tessellation.hpp>
#include <OT/intersection.hpp>
#include <OT/Density_base_2.hpp>

namespace OT
{
  template <class DT, class Density>
  bool
  centroid(const Density &density,
	   const DT &dt, typename DT::Vertex_handle v,	   
	   std::pair<typename Density::Point, double> &p) 
  {
    typedef typename DT::Geom_traits::Kernel K;
    typename CGAL::Polygon_2<K> poly;

    if (!tessellate_voronoi_cell(dt, v, std::back_inserter(poly),
				density.clamp_radius()))
      return false;

    return OT::centroid(density, poly, p);
  }

  template <class DT, class Density>
  double
  mass(const Density &density,
       const DT &dt, typename DT::Vertex_handle v) 
  {
    typedef typename DT::Geom_traits::Kernel K;
    typename CGAL::Polygon_2<K> poly;

    if (!tessellate_voronoi_cell(dt, v, std::back_inserter(poly),
				density.clamp_radius()))
      return 0.0;

    return OT::mass(density, poly);
  }

  template <class DT, class Density, class Function>
  double
  integrate(const Density &density,
	    const DT &dt, typename DT::Vertex_handle v,
	    Function f)
  {
    typedef typename DT::Geom_traits::Kernel K;
    typename CGAL::Polygon_2<K> poly;

    if (!tessellate_voronoi_cell(dt, v, std::back_inserter(poly),
				density.clamp_radius()))
      return 0.0;

    return OT::integrate(density, f, poly);
  }
}
#endif
