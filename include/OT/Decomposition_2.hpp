#ifndef DECOMPOSITION_HPP
#define DECOMPOSITION_HPP

#include <OT/Discrete_measure_2.hpp>
#include <OT/Lloyd_energy_2.hpp>
#include <OT/bfgs_descent.hpp>
#include <OT/lloyd.hpp>
#include <OT/misc.hpp>

namespace OT
{
  template <class K>
  class Decomposition_2
  {
  public:
    typedef typename OT::Discrete_measure_2<K> Measure;
    typedef typename std::vector<size_t> Map_to_parent;
    typedef typename K::Point_2 Point;

  private:
    std::vector<Measure> _levels;
    std::vector<Map_to_parent> _maps;

  public:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & _levels;
        ar & _maps;
    }

    void output(std::ostream &os) const
    {
      os << _levels.size() << "\n";
      for(size_t i = 0; i < _levels.size(); ++i)
	_levels[_levels.size()-1-i].output(os);
    }
       
  public:
    Decomposition_2 () {}
    Decomposition_2(const Measure &measure) { reset(measure); }

    void
    reset (const Measure &meas)
    {
      _levels.clear();
      _maps.clear();
      _levels.push_back(meas);
    }

    const Measure &
    level(size_t i) const
    {
      OT_ASSERT(i < _levels.size());
      return _levels[i];
    }

    size_t
    number_of_levels() const
    {
      return _levels.size();
    }

    void precompute_next_level (size_t k, double threshold_pc)
    {    
      typedef typename CGAL::Delaunay_triangulation_2<K> DT;
      typedef typename DT::Vertex_handle Vertex_handle;
      typedef typename DT::Finite_vertices_iterator Iterator;
	
      OT_ASSERT(_levels.size() > 0);
	
      size_t L = _levels.size(), l = L - 1;
      Measure &current_level = _levels[l];

      double best_energy = 1e6;

      DT dt;
      boost::mt19937 engine;
      size_t Ntrials = 5;
      for (size_t i = 0; i < Ntrials; ++i)
	{
	  DT cur_dt;
	  lloyd_init (cur_dt, current_level, level(l).size()/k, engine);
	  
	  double threshold = 0; //threshold_pc * current_level.clamp_radius();
	  while (lloyd_step(cur_dt, current_level) > threshold)
	    ;

	  double energy = lloyd_energy(cur_dt, current_level);
	  if (energy < best_energy)
	    {
	      best_energy = energy;
	      dt = cur_dt;
	    }
	  //OT_DEBUG_SHOW(energy);
	}
      //OT_DEBUG_SHOW(best_energy);
      
      // clean useless points
      lloyd_step (dt, current_level, true);

      // copy back information to _levels[L], _maps[l]
      std::map<Vertex_handle, size_t> indices;

      size_t N = dt.number_of_vertices();
      std::vector<Point> positions (N);
      std::vector<double> masses (N);
      std::vector<size_t> map;

      size_t id = 0;
      for (Iterator it = dt.finite_vertices_begin();
	   it != dt.finite_vertices_end(); ++it, ++id)
	{
	  indices[it] = id;
	  positions[id] = it->point();
	}
      
      for (size_t i = 0; i < current_level.size(); ++i)
      {
	Vertex_handle v = dt.nearest_vertex(current_level.position(i));

	size_t id = indices[v];
	masses[id] += current_level.mass(i);
	map.push_back(id);
      }

      _maps.push_back (map);
      _levels.push_back (Measure (masses, positions));
    }

    void precompute_levels(size_t L,
			   size_t k,
			   double threshold = 0.01)
    {
      OT_ASSERT(_levels.size() > 0);

      for (size_t l = _levels.size() - 1; l < L; ++l)
	precompute_next_level(k, threshold);
    }

    size_t
    ancestor(size_t l, size_t L, size_t i) const
    {
      while (l < L)
	{
	  i = _maps[l][i];
	  l++;
	}

      return i;
    }

    template <class Vector>
    void prepare_weights (size_t l,
			  Vector &to,
			  size_t L,
			  const Vector &from) const
    {
      OT_ASSERT (l < _levels.size());

      size_t N = level(l).size();

      to.resize(N);
      if ( (L >= _maps.size()) || (from.size() != level(L).size()) )
	{
	  std::fill(to.begin(), to.end(), 0.0);
	  return;
	}

      for (size_t i = 0; i < N; ++i)
	to[i] = from[ancestor(l, L, i)];
    }

    double wasserstein_distance(size_t l, size_t L) const
    {
      if (l == L)
	return 0;

      size_t N = level(l).size();

      double w = 0.0;
      for (size_t i = 0; i < N; ++i)
	{
	  size_t j = ancestor(l, L, i);
	  double dd = CGAL::squared_distance(level(l).position(i),
					     level(L).position(j));
	  w += level(l).mass(i) * dd;
	}
      return sqrt(w);
    }
  };
}

#endif
