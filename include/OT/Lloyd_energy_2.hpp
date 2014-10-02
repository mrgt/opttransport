#ifndef OT_LLOYD_ENERGY_2_HPP
#define OT_LLOYD_ENERGY_2_HPP

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/squared_distance_2.h>

namespace OT
{
  template <class RT> class Lloyd_energy;

  template <class Function>
  struct Lloyd_energy_workspace
  {
    typedef typename Function::Delaunay_triangulation DT;
    typedef typename Function::Measure Measure;
    typedef typename DT::Vertex_handle Vertex_handle;
    typedef typename DT::Point Point;

  protected:
    DT _dt;
    const Function *_f;
    std::vector<Point> _positions;
    std::vector<double> _gradient;
    double _energy;
    size_t _size;

  public:
    typedef std::vector<double> Gradient_type;
    
  public:
    Lloyd_energy_workspace ()
    {}

    template <class InputIterator>
    void prepare(const Function *f, InputIterator begin, InputIterator end);

    size_t dimension()
    {
      return _gradient.size();
    }

    const Gradient_type &gradient() const
    {
      return _gradient;
    }

    double value() const
    {
      return _energy;
    }
  };
  
  template <class DT>
  class Lloyd_energy
  {
  public:
    typedef DT Delaunay_triangulation;
    typedef typename DT::Geom_traits::Kernel K;
    typedef typename K::Iso_rectangle_2 Iso_rectangle;
    typedef typename DT::Point Point;
    typedef typename OT::Lloyd_energy<DT> Self;
    typedef typename OT::Lloyd_energy_workspace<Self> Workspace;
    friend class OT::Lloyd_energy_workspace<Self>;
    typedef typename OT::Discrete_measure_2<K> Measure;

  protected:    
    const Measure &_measure;
    
  public:
    Lloyd_energy(const Discrete_measure_2<K> &measure) :
      _measure(measure)
    {}

    template <class InputIterator>
    bool eval (Workspace &wsp,
	       InputIterator begin,
	       InputIterator end) const
    {
      wsp.prepare (this, begin, end);
    }
  };

  template <class Function>
  template <class InputIterator>
  void
  Lloyd_energy_workspace<Function>::prepare
      (const Function *f,
       InputIterator begin,
       InputIterator end)
  {
    std::cerr << "here!!\n";
    _positions.clear();
    for (; begin != end; begin += 2)
	_positions.push_back(Point (*begin, *(begin+1)));

    _dt.clear();
    _dt.insert (_positions.begin(), _positions.end());

    // get a map from Vertex_handle to indices 
    std::map<Vertex_handle, size_t> indices;
    for (size_t i = 0; i < _positions.size(); ++i)
      {
	Vertex_handle nn = _dt.nearest_vertex (_positions[i]);
	indices[nn] = i;
      }

    _size = _dt.number_of_vertices();
    _f = f;
    _gradient.resize(_size * 2);
    std::fill (_gradient.begin(), _gradient.end(), 0);

    const Measure &meas = _f->_measure;
    _energy = 0.0;

    for (size_t i = 0; i < meas.size(); ++i)
      {
	Vertex_handle nn = _dt.nearest_vertex (meas.position(i));
	size_t id = indices[nn];

	_energy += meas.mass(i) * CGAL::squared_distance (nn->point(), meas.position(i));
	_gradient[2 * id] += 2.0 * (nn->point().x() - meas.position(i).x());
	_gradient[2 * id + 1] += 2.0 * (nn->point().y() - meas.position(i).y());
      }

    OT_DEBUG_SHOW (_energy);
  }

}

#include <OT/bfgs_descent.hpp>

namespace OT
{
  template <class DT, class Measure>
  void
  lloyd_bfgs_descent (DT &dt, const Measure &mu)
  {
    typedef typename DT::Point Point;
    typedef typename DT::Finite_vertices_iterator Iterator;

    uvector positions (2 * dt.number_of_vertices());
    size_t i = 0;

    for (Iterator it = dt.finite_vertices_begin();
	 it != dt.finite_vertices_end(); ++it, ++i)
      {
	positions[2*i] = it->point().x();
	positions[2*i+1] = it->point().y();
      }

    OT::bfgs_descent (Lloyd_energy<DT>(mu), 
		      positions, 1e-5);

    dt.clear();
    for (size_t i = 0; i < positions.size(); i += 2)
	dt.insert (Point(positions[i], positions[i+1]));
  }
    
}

#endif
