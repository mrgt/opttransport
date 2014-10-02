#ifndef OT_TESSELLATION_HPP
#define OT_TESSELLATION_HPP

#include <OT/misc.hpp>

namespace OT
{
  namespace details
  {
    template <class RT, class V>
    bool is_hidden(const RT &rt, const V &v)
    {
      return v->is_hidden();
    }
    
    template <class V, class K>
    bool is_hidden(const CGAL::Delaunay_triangulation_2<K> &dt,
		   const V &v)
    {
      return false;
    }
  }
  
  template <class DT, class OutputIterator>
  bool
  tessellate_voronoi_cell(const DT &dt,
			  typename DT::Vertex_handle v,
			  OutputIterator it,
			  double R)
  {
    typedef typename DT::Edge_circulator Edge_circulator;
    typedef typename DT::Point Point;
    typedef typename DT::Segment Segment;
    typedef typename DT::Geom_traits::Kernel K;
    typedef typename DT::Geom_traits::Ray_2 Ray;
    typedef typename DT::Geom_traits::Vector_2 Vector;

    OT_ASSERT(R > 0.0);

    if (OT::details::is_hidden(dt, v))
      return false;

    // First step: determine whether the Voronoi cell is infinite.
    // If this is the case, we want to start with c = infinite edge.
    Edge_circulator c = dt.incident_edges (v), cm = c--, done(c);    
    bool adjacent_to_infinite = false;
    do
      {
	if (dt.is_infinite(c))
	  {
	    adjacent_to_infinite = true;
	    break;
	  }
      } while (++c != done);
    done = c;
    
    bool after_infinite = true;
    //#define SIMPLE_DEBUG

#ifdef SIMPLE_DEBUG
    CGAL::Polygon_2<K> P;
#endif
    do
      {
	if (dt.is_infinite(c))
	  ++c;

	CGAL::Object o = dt.dual(c);

	Segment s;
	if (CGAL::assign (s, o) )
	  {
#ifdef SIMPLE_DEBUG
	    P.push_back(s.source());
	    std::cerr << "segment ";
	    //OT_DEBUG_SHOW(c->first);
	    OT_DEBUG_SHOW(c->second);
	    OT_DEBUG_SHOW(s.source());
	    OT_DEBUG_SHOW(s.target());
#endif
	    if (s.source() != s.target())
	      *it ++ = s.source();
	    continue;
	  }

	Ray r;
	if (!CGAL::assign(r,o))
	  {
	    // FIXME only possible case left: dual is a line. we
	    // currently neglect this case
	    return false;
	  }

	Vector to = r.to_vector();
	to = (R/sqrt(to.squared_length())) * to;
	
	if (after_infinite)
	  {
	    *it++ = r.source() + to;
#ifdef SIMPLE_DEBUG
	    P.push_back(r.source() + to);
	    std::cerr << "after infinite ";
	    OT_DEBUG_SHOW(r.source() + to);
#endif
	    after_infinite = false;
	  }
	else
	  {
#ifdef SIMPLE_DEBUG
	    P.push_back(r.source());
	    P.push_back(r.source() + to);
	    std::cerr << "before infinite ";
	    OT_DEBUG_SHOW(r.source());
	    OT_DEBUG_SHOW(r.source() + to);
#endif
	    *it++ = r.source();
	    *it++ = r.source() + to;
	  }
      }
    while (++c != done);

#ifdef SIMPLE_DEBUG
    OT_DEBUG_SHOW(P.is_simple());
#endif

    return true;
  }
}

#endif
