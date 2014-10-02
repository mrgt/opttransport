#ifndef OT_OBJECTIVE_FUNCTION_HPP
#define OT_OBJECTIVE_FUNCTION_HPP

#include <OT/random.hpp>
#include <CGAL/bounding_box.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>

namespace OT
{

  // return the average squared distance from a point in the triangle
  // abc to the point p, multiplied by the area of abc.
  inline double
  trace_covariance_2(double a1, double a2, 
		     double b1, double b2,
		     double c1, double c2,
		     double p1, double p2)
  {
    double u1 = b1 - a1, u2 = b2 - a2;
    double v1 = c1 - a1, v2 = c2 - a2;
    double area = fabs(u1 * v2 - u2 * v1)/2.0;

    double s1 = a1 + b1 + c1;
    double s2 = a2 + b2 + c2;
    double tr0 = (area/12.0) * ( (a1*a1 + a2*a2) + (b1*b1 + b2*b2) +
				 (c1*c1 + c2*c2) + (s1*s1 + s2*s2) );

    return tr0 + ( area * p1 * (p1 - 2.0/3.0 * s1) +
		   area * p2 * (p2 - 2.0/3.0 * s2) );
  }	       
  
  
  // Class defining the function f: x -> || x - p ||^2 
  template <class Point>
  class Distance_squared_2
  {
    Point p;
    
  public:
    Distance_squared_2 (const Point &point) : p(point) {}
    
    double operator () (int i, int j) const
    {
      double dx = p.x() - (double)i;
      double dy = p.y() - (double)j;
      return dx * dx + dy * dy;
    }
    
    // returns the integral of our function over the triangle abc.
    double operator () (const Point &a, const Point &b, const Point &c) const
    {
      return trace_covariance_2(a.x(), a.y(), b.x(), b.y(),
				c.x(), c.y(), p.x(), p.y());
    }

    // returns the integral of our function over the segment ab.
    //
    // | (1-t) a + t b - p |^2
    // = | t (b - a) + a - p |^2
    // = t^2 |b-a|^2 + 2 t <b-a><a-p> + |a-p|^2
    // => with \int_0^1 t^2 = 1/3 and \int_0^1 t = 1/2
    //    the integral is length(ab) * (1/3 |b-a|^2 + <b - a><a - p> + |a-p|^2)
    double operator () (const Point &a, const Point &b) const
    {
      double ddab = CGAL::squared_distance(a,b) ;
      double avg = 1.0/3.0 * ddab + (b - a) * (a - p)
	+ CGAL::squared_distance(a,p);
      return  sqrt(ddab) * avg;
    }
  };
  
  template <class RT>
  void
  rt_copy_vertices_in_order
      (const RT &rt,
       const std::map<typename RT::Bare_point, size_t> &indices,
       std::vector<typename RT::Vertex_handle> &vertices)
  {
    vertices.resize(indices.size(), 0);

    typedef typename RT::Finite_vertices_iterator Finite_vertices_iterator;
    typedef typename RT::Hidden_vertices_iterator Hidden_vertices_iterator;
    typedef typename RT::Vertex_handle Vertex_handle;
    typedef typename std::map<typename RT::Bare_point, size_t> Map;
    
    for (Finite_vertices_iterator it = rt.finite_vertices_begin();
	 it != rt.finite_vertices_end(); ++it)
      {
	typename Map::const_iterator jt = indices.find(it->point().point());
	OT_ASSERT(jt != indices.end());
	vertices[jt->second] = it;
      }
    for (Hidden_vertices_iterator it = rt.hidden_vertices_begin();
	 it != rt.hidden_vertices_end(); ++it)
      {
	typename Map::const_iterator jt = indices.find(it->point().point());
	OT_ASSERT(jt != indices.end());
	vertices[jt->second] = it;
      }
  }

  template <class RT, class PositionsIterator, class WeightsIterator>
  void
  rt_build (RT &rt,
	    PositionsIterator pbegin, PositionsIterator pend,
	    WeightsIterator wbegin, WeightsIterator wend)
  {
    typedef typename RT::Point Point;
    std::vector<Point> points;

    for (; (pbegin != pend) && (wbegin != wend); ++wbegin, ++pbegin)
      points.push_back(Point (*pbegin, *wbegin));
    
    rt.clear();
    rt.insert(points.begin(), points.end());
  }

  template <class K, class PointIterator, class MassIterator>
  double bound_wasserstein_distance_greedy (PointIterator pbegin, 
					    PointIterator pend,
					    MassIterator mbegin,
					    MassIterator mend)
  {
    typedef CGAL::Triangulation_face_base_2<K> Fb;
    typedef CGAL::Triangulation_vertex_base_with_info_2<double,K> Vb;
    typedef CGAL::Triangulation_data_structure_2<Vb,Fb> Tds;
    typedef typename CGAL::Delaunay_triangulation_2<K, Tds> DT;
    typedef typename DT::Vertex_handle Vertex_handle;
    typedef typename DT::Point Point;

    DT dt;

    PointIterator pit = pbegin;
    MassIterator mit = mbegin;

    std::vector< std::pair<Point, double> > source_points;

    double positive(0), negative(0);
    for (; (pit != pend) && (mit != mend); ++pit, ++mit)
      {
	if (*mit > 0)
	  {
	    source_points.push_back(std::make_pair(*pit, *mit));
	    positive += *mit;
	  }
	else
	  {
	    Vertex_handle v = dt.insert(*pit);
	    v->info() = - (*mit);
	    negative -= *mit;
	  }
      }
    //OT_DEBUG_SHOW(positive);
    //OT_DEBUG_SHOW(negative);

    double cost = 0.0;
    double total = 0.0;
    size_t i = 0;
    while ( (i < source_points.size()) && dt.number_of_vertices() > 0)
      {
	Vertex_handle nn = dt.nearest_vertex (source_points[i].first);

	if (nn->info() > source_points[i].second)
	  {
	    cost += source_points[i].second
	      * CGAL::squared_distance (source_points[i].first, nn->point());
	    total += source_points[i].second;
	    nn->info() -= source_points[i].second;
	    ++i;
	  }
	else
	  {
	    cost += nn->info()
	      * CGAL::squared_distance (source_points[i].first, nn->point());	    
	    total += nn->info();
	    source_points[i].second -= nn->info();
	    dt.remove(nn);
	  }
      }
    //OT_DEBUG_SHOW(total);

    return sqrt(cost);
  }

  template <class RT, class Density> class Objective_function;

  template <class Function>
  struct Objective_function_workspace
  {
    typedef typename Function::Regular_triangulation RT;
    typedef typename RT::Geom_traits::Kernel K;
    typedef typename RT::Vertex_handle Vertex_handle;
    typedef typename RT::Bare_point Bare_point;
    typedef typename RT::Point Point;

  protected:
    RT _rt;
    const Function *_f;
    std::vector<Vertex_handle> _vertices;
    std::vector<double> _diff_masses, _gradient;
    size_t _num_active;

    mutable bool _has_wasserstein;
    mutable double _wasserstein;

  public:
    typedef std::vector<double> Gradient_type;
    
  public:
    Objective_function_workspace ()
    {}

    template <class InputIterator>
    void prepare(const Function *f, InputIterator begin, InputIterator end);

    const Gradient_type &gradient() const
    {
      return _gradient;
    }

    double value() const
    {
      double value = squared_wasserstein();

      typename std::vector<Vertex_handle>::const_iterator jt =
	_vertices.begin();

      for (size_t i = 0; jt != _vertices.end(); ++jt, ++i)
	value += - _diff_masses[i] * (*jt)->point().weight();

      return -value;
    }

    template <class Matrix>
    bool hessian(Matrix &m, std::vector<size_t> &active) const
    {
      // Inverse map
      std::map<Vertex_handle, size_t> indices;
      for (size_t i = 0; i < _vertices.size(); ++i)
	  indices[_vertices[i]] = i;

      for (size_t i = 0; i < _vertices.size(); ++i)
	{
	  Vertex_handle v = _vertices[i];

	  if (v->is_hidden())
	    continue;

	  typename RT::Edge_circulator c = _rt.incident_edges (v), done(c);
	  std::map<size_t, double> entries;

	  double total = 0.0;

	  do
	    {
	      Vertex_handle w = c->first->vertex(_rt.ccw(c->second));
	      size_t iw = indices[w];

	      OT_ASSERT(c->first->vertex(_rt.cw(c->second)) == v);
	      typename RT::Segment s;
	      typename RT::Geom_traits::Ray_2 r;

	      if (_rt.is_infinite(c))
		continue;
	     
	      CGAL::Object o = _rt.dual(c); 
	      double mass = 0.0;
	      if (CGAL::assign (s, o))
		  mass = OT::mass (_f->_density, s);
	      else if (CGAL::assign(r,o))
		  mass = OT::mass (_f->_density, r);
	      else
		std::cerr << "[vw] has no dual ?\n";

	      // std::cerr << "(" << i << ", " << iw << ") -> "
	      // 		<< mass << "\n";

	      double sqd = CGAL::squared_distance(w->point(), v->point());
	      double d2vw = mass / (2 * sqrt(sqd));
	      entries[iw] = -d2vw;
	      total += d2vw;
	    }
	  while (++c != done);
	  
	  if (total == 0.0)
	    continue;

	  entries[i] = total;
	  active.push_back(i);
	  std::map<size_t, double>::iterator it;
	  for (it = entries.begin(); it != entries.end(); ++it)
	    {
	      m.push_back(it->first, i,  it->second);
	    }
	}

      return true;
    }

    size_t num_active() const
    {
      return _num_active;
    }

    double squared_wasserstein() const
    {
      if (_has_wasserstein)
	return _wasserstein;

      double value = 0.0;

      typename std::vector<Vertex_handle>::const_iterator jt =
	_vertices.begin();

      for (size_t i = 0; jt != _vertices.end(); ++jt, ++i)
	{
	  Vertex_handle v = *jt;
	  Distance_squared_2<Bare_point> sqdist (v->point().point());
	  value += OT::integrate (_f->_density, _rt, v, sqdist);	  
	}

      _wasserstein = value;
      _has_wasserstein = true;
      
      return _wasserstein;
    }

    double wasserstein() const
    {
      return sqrt(squared_wasserstein());
    }

    double wasserstein_error() const
    {
      double greedy_bound =  OT::bound_wasserstein_distance_greedy<K>
	(_f->_positions.begin(), _f->_positions.end(),
	 _diff_masses.begin(), _diff_masses.end());
      double trivial_bound = sqrt(misaffected_mass()) * _f->diameter();
      //OT_DEBUG_SHOW(misaffected_mass() * _f->diameter());
      return std::min(greedy_bound, trivial_bound);
    }

    double misaffected_mass() const
    {
      // compute erroneously affected mass
      double m = 0.0;
      for (size_t i = 0; i < _diff_masses.size(); ++i)
	m += fabs(_diff_masses[i]);

      return m / 2.0;
    }
  };
  
  template <class RT, class D>
  class Objective_function
  {
  public:
    typedef D Density;
    typedef RT Regular_triangulation;
    typedef typename RT::Geom_traits::Kernel K;
    typedef typename K::Iso_rectangle_2 Iso_rectangle;
    typedef typename RT::Point Point;
    typedef typename RT::Bare_point Bare_point;
    typedef typename OT::Objective_function<RT, Density> Self;
    typedef typename OT::Objective_function_workspace<Self> Workspace;
    friend class OT::Objective_function_workspace<Self>;

  protected:    
    const Density &_density;
    const std::vector<Bare_point> &_positions;
    const std::vector<double> &_masses;
    std::map<Bare_point, size_t> _indices;
    Iso_rectangle _bbox;
    double _quantization_error;
    
  public:
    Objective_function(const Density &density,		       
		       const std::vector<Bare_point> &positions,
		       const std::vector<double> &masses) :
      _density(density), _positions(positions), _masses(masses),
      _quantization_error(0)
    {
      OT_ASSERT(positions.size() == masses.size());

      for (size_t i = 0; i < _positions.size(); ++i)
	_indices[_positions[i]] = i;

      _bbox = CGAL::bounding_box(_positions.begin(), _positions.end());
    }

    template <class InputIterator>
    bool eval (Workspace &wsp,
	       InputIterator begin,
	       InputIterator end) const
    {
      wsp.prepare (this, begin, end);
    }

    size_t number_of_vertices () const
    {
      return _positions.size();
    }

    // returns an upper bound on the diameter of the support of the
    // target measure
    double diameter() const
    {
      return sqrt( pow(_bbox.xmax () - _bbox.xmin (), 2.0) +
		   pow(_bbox.ymax () - _bbox.ymin (), 2.0) );
    }

    void set_quantization_error(double err)
    {
      _quantization_error = err;
    }

    double quantization_error() const
    {
      return _quantization_error;
    }
  };

  template <class Function>
  template <class InputIterator>
  void
  Objective_function_workspace<Function>::prepare
      (const Function *f,
       InputIterator begin,
       InputIterator end)
  {
    _f = f;
    _num_active = 0;
    _has_wasserstein = false;
    _wasserstein = 0;

    size_t N = _f->number_of_vertices();
    size_t Nw = end - begin;

    OT_ASSERT (Nw <= N);
  
    // Recompute regular triangulation with new weights.
    rt_build (_rt,
	      _f->_positions.begin(), _f->_positions.begin() + Nw,
	      begin, end);
    
    for (size_t i = Nw; i < N; ++i)
      _rt.insert (Point (_f->_positions[i], 0.0));
    
    // Extract vertices and reorder
    OT::rt_copy_vertices_in_order(_rt, _f->_indices, _vertices);
    
    // Pre-compute gradient
    typename std::vector<Vertex_handle>::iterator jt = _vertices.begin();

    _diff_masses.resize (N);

    for (size_t i = 0; jt != _vertices.end(); ++jt, ++i)
      {
	Vertex_handle v = *jt;

	const double target_mass = f->_masses[i];
	const double current_mass = OT::mass(f->_density, _rt, v);
	_diff_masses[i] = - (target_mass - current_mass);

	if (current_mass != 0.0)
	  _num_active++;
      }

    _gradient.resize(Nw);
    for (size_t i = 0; i < Nw; ++i)
	_gradient[i] = _diff_masses[i]; // - _diff_masses[N-1];

    //std::cerr << n_active << " active sites\n";
  }
}

#endif
