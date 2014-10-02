#ifndef LLOYD_HPP
#define LLOYD_HPP

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <OT/Discrete_measure_2.hpp>

#include <OT/random.hpp>
#include <OT/misc.hpp>
#include <OT/integration.hpp>

namespace OT
{

  template <class K, class Engine>
  void
  lloyd_init(CGAL::Delaunay_triangulation_2<K> &dt,
	     const Discrete_measure_2<K> &meas,
	     size_t k,
	     Engine &engine)
  {
    typedef typename CGAL::Delaunay_triangulation_2<K> DT;
    typedef typename DT::Vertex_handle Vertex_handle;
    typedef typename DT::Finite_vertices_iterator Iterator;
    typedef typename K::Point_2 Point;

    //std::vector<Point> positions;
    do
      {
	Point p = OT::random_point(meas, engine);
	dt.insert(p);
      } while (dt.number_of_vertices() != k);

    // dt.insert(positions.begin(), positions.end());
  }


  template <class K>
  double
  lloyd_energy(CGAL::Delaunay_triangulation_2<K> &dt,
	       const Discrete_measure_2<K> &meas)
  {
    typedef typename CGAL::Delaunay_triangulation_2<K> DT;
    typedef typename DT::Vertex_handle Vertex_handle;

    double energy (0);
    for (size_t i = 0; i < meas.size(); ++i)
      {
	Vertex_handle v = dt.nearest_vertex(meas.position(i));
	energy += meas.mass(i)*CGAL::squared_distance(v->point(), meas.position(i));
      }
    return energy;
  }

  template <class K>
  double
  lloyd_step(CGAL::Delaunay_triangulation_2<K> &dt,
	     const Discrete_measure_2<K> &meas,
	     bool cleaning_move = false)
  {
    typedef typename CGAL::Delaunay_triangulation_2<K> DT;
    typedef typename DT::Vertex_handle Vertex_handle;
    typedef typename DT::Finite_vertices_iterator Iterator;
    typedef typename K::Point_2 Point;
    typedef typename K::Vector_2 Vector;

    size_t k = dt.number_of_vertices();

    std::map<Vertex_handle, size_t> indices;
    std::vector<Vertex_handle> vertices;

    size_t id = 0;
    for (Iterator it = dt.finite_vertices_begin();
	 it != dt.finite_vertices_end(); ++it, ++id)
      {
	indices[it] = id;
	vertices.push_back (it);
      }
	
    std::vector<Vector> translations (k, Vector(0,0));
    std::vector<double> masses (k, 0.0);

    for (size_t i = 0; i < meas.size(); ++i)
      {
	Vertex_handle v = dt.nearest_vertex(meas.position(i));
	size_t id = indices[v];

	masses[id] += meas.mass(i);
	translations[id] = translations[id] + 
	  meas.mass(i) * (meas.position(i) - CGAL::ORIGIN);
      }

    double max_move = 0.0;
    double energy = 0.0;
    std::vector<Point> new_positions;

    for (size_t id = 0; id < k; ++id)
      {
	Point position;

	//OT_DEBUG_SHOW(masses[id]);
	if (masses[id] > 0)
	  {
	    position = CGAL::ORIGIN + (1.0/masses[id]) * translations[id];
	    new_positions.push_back(position);
	    //OT_DEBUG_SHOW(position - vertices[id]->point());
	  }
	else 
	  {
	    position = vertices[id]->point();
	    
	    if (!cleaning_move)
	      new_positions.push_back(position);
	  }

	max_move = std::max(max_move,
			    CGAL::squared_distance(position,
						   vertices[id]->point()));
      }
    
    dt.clear();
    dt.insert(new_positions.begin(), new_positions.end());
    
    max_move = sqrt(max_move);
    //OT_DEBUG_SHOW(k);
    //OT_DEBUG_SHOW(max_move);

    return max_move;
  }

  template <class RT, class Density, class Engine>
  double
  lloyd_step(RT &dt,
	     Density &density,
	     Engine &engine)
  {
    typedef typename RT::Point Point;

    std::vector<Point> new_positions;    
    double max_move = 0.0;

    typename RT::Finite_vertices_iterator it;
    size_t i = 0;
    for (it = dt.finite_vertices_begin(); it != dt.finite_vertices_end(); ++it)
      {
	std::pair<Point, double> p;
	Point oldp = it->point();
	
	if (!OT::centroid(density, dt, it, p))
	  {
	    //new_positions.push_back(random_point(density, engine));
	    new_positions.push_back(oldp);
	    //std::cerr << "lloyd_step: no centroid\n";
	    continue;
	  }
	
	new_positions.push_back(p.first);
	max_move =
	  std::max(max_move, sqrt(CGAL::squared_distance(oldp, p.first)));
      }
    //OT_DEBUG_SHOW(max_move);
      
    dt.clear();
    dt.insert(new_positions.begin(), new_positions.end());

    return max_move;
  }

  template <class K, class RT, class Density, class Engine,
	    class OutputPositionsIterator,
	    class OutputMassesIterator,
	    class OutputPolygonsIterator>
  size_t 
  k_means (RT &rt,
	   Density density,
	   Engine &engine,
	   OutputPositionsIterator out_pos,
	   OutputMassesIterator out_mass,
	   OutputPolygonsIterator out_poly,
	   double stop_threshold_pct = 0.001,
	   size_t max_iter = 100)
  {
    double stop_threshold = stop_threshold_pct * density.clamp_radius();
    size_t N = 0;

    while (lloyd_step(rt, density, engine) >= stop_threshold &&
	   N < max_iter)
      ++N;

    typename RT::Finite_vertices_iterator it;

    for (it = rt.finite_vertices_begin();
	 it != rt.finite_vertices_end(); ++it, ++N)
      {
	if (!is_null_output_iterator(out_pos))
	  *out_pos ++ = it->point();
	if (!is_null_output_iterator(out_mass) ||
	    !is_null_output_iterator(out_poly))
	  {
	    typename CGAL::Polygon_2<K> poly, cpoly;	      
	    OT::tessellate_voronoi_cell(rt, it,
					std::back_inserter(poly),
					density.clamp_radius());
	    OT::convex_intersection(poly, density.bounding_poly(), cpoly);	

	    if (!is_null_output_iterator(out_mass))
	      *out_mass ++ = OT::mass(density, poly);
	    if (!is_null_output_iterator(out_poly))
	      *out_poly++ = cpoly;
	  }
      }

    return N;
  }

  // stop_threshold_pct is a proportion of density.clamp_radius();
  template <class K, class Density,
	    class OutputPositionsIterator,
	    class OutputMassesIterator,
	    class OutputPolygonsIterator>
  size_t 
  k_means (size_t k,
	   Density density,
	   OutputPositionsIterator out_pos,
	   OutputMassesIterator out_mass,
	   OutputPolygonsIterator out_poly,
	   double stop_threshold_pct = 0.001,
	   size_t max_iter = 100)
  {
    typedef typename CGAL::Delaunay_triangulation_2<K> RT;
    boost::mt19937 engine;

    RT rt;
    for (size_t i = 0; i < k; ++i)
      rt.insert(OT::random_point(density, engine));

    return OT::k_means<K>(rt, density, engine, out_pos, out_mass, out_poly,
			  stop_threshold_pct, max_iter);
  }

  template <class K, class Density,
	    class OutputPositionsIterator>
  size_t 
  k_means (size_t k,
	   Density density,
	   OutputPositionsIterator out_pos,
	   double stop_threshold_pct = 0.001)
  {
    return OT::k_means<K>(k, density, out_pos, 
			  null_output_iterator(),
			  null_output_iterator(),
			  stop_threshold_pct);
  }

  // stop_threshold_pct is a proportion of density.clamp_radius();
  template <class K, class Density,
	    class InputPositionsIterator,
	    class OutputPositionsIterator,
	    class OutputMassesIterator,
	    class OutputPolygonsIterator>
  size_t 
  k_means (InputPositionsIterator begin,
	   InputPositionsIterator end,		   
	   Density density,
	   OutputPositionsIterator out_pos,
	   OutputMassesIterator  out_mass,
	   OutputPolygonsIterator out_poly,
	   double stop_threshold_pct = 0.001)
  {
    typedef typename CGAL::Delaunay_triangulation_2<K> RT;
    boost::mt19937 engine;

    RT rt;
    rt.insert(begin, end);

    return OT::k_means<K>(rt, density, engine, out_pos, out_mass, out_poly,
			  stop_threshold_pct);
  }


}

#endif
