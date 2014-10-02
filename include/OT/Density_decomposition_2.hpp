#ifndef DENSITY_DECOMPOSITION_HPP
#define DENSITY_DECOMPOSITION_HPP

#include <OT/Density_base_2.hpp>
#include <OT/lloyd.hpp>
#include <OT/misc.hpp>
#include <OT/bfgs_descent.hpp>
//#include <OT/bundle_descent.hpp>
#include <OT/integration.hpp>
#include <OT/newton_descent.hpp>
#include <OT/objective_function.hpp>


namespace OT
{
  template <class K>
  struct Density_decomposition_level
  {
    typedef typename K::Point_2 Point;
    typedef typename CGAL::Polygon_2<K> Polygon;
    
    std::vector<Point> _positions;
    std::vector<double> _masses;
    std::vector<Polygon> _polygons;
    std::vector<size_t> _parents; 

    size_t size() const
    {
      return _positions.size();
    }  
  };


  template<class K>
  typename CGAL::Polygon_2<K>
  enlarge_polygon(const typename CGAL::Polygon_2<K> &p)
  { 
    const CGAL::Bbox_2 bbox = p.bbox();
    typename CGAL::Polygon_2<K> r;
    r.push_back(typename K::Point_2(bbox.xmin() - 1.0, bbox.ymin() - 1.0));
    r.push_back(typename K::Point_2(bbox.xmax() + 1.0, bbox.ymin() - 1.0));
    r.push_back(typename K::Point_2(bbox.xmax() + 1.0, bbox.ymax() + 1.0));
    r.push_back(typename K::Point_2(bbox.xmin() - 1.0, bbox.ymax() + 1.0));
    return r;
  }


  template <class K, class Density>
  class Density_decomposition
  {
  public:
    typedef typename OT::Density_decomposition_level<K> Level;
    typedef typename K::Point_2 Point;
    typedef typename CGAL::Polygon_2<K> Polygon;

  private:
    Density _density;
    std::vector<Level> _levels;
       
  public:
    Density_decomposition (const Density &density = Density()):
      _density(density)
    {}

    const Level &
    level(size_t i)
    {
      OT_ASSERT(i < _levels.size());
      return _levels[i];
    }

    void compute_next_level(size_t k = 10,
			    double threshold = 0.01)
    {
      size_t N = _levels.size();
      _levels.resize(N+1);
      Level &next_level = _levels[N];

      if (N == 0)
	{
	  std::pair<Point, double> c;
	  OT_ASSERT(OT::centroid (_density, _density.bounding_poly(),
				  c, true));
	  next_level._positions.push_back(c.first);
	  next_level._masses.push_back(c.second);
	  next_level._parents.push_back(size_t(-1));
	  next_level._polygons.push_back(_density.bounding_poly());
	  return;
	}

#if 0
      const Level &level = _levels[N-1];
      double total_mass = 0.0;

      for (size_t i = 0; i < level.size(); ++i)
	{
	  size_t Nsteps = OT::k_means<K>
	    (k, Density(_density, level._polygons[i]),
	     std::back_inserter(next_level._positions),
	     std::back_inserter(next_level._masses),
	     std::back_inserter(next_level._polygons));

	  // set correct parent and renormalize masses
	  for (size_t j = 0; j < k; ++j)
	    next_level._parents.push_back(i);

	  for (size_t j = 0; j < k; ++j)
	    {
	      next_level._masses[i * k + j] *= level._masses[i];
	      total_mass += next_level._masses[i * k + j];
	    }
	}
      
#else
      const Level &level = _levels[N-1];


  // std::vector<Point_2> kpos;
  // std::cerr << "begin k_means\n";
  // OT::k_means<K>(1000, density_to, std::back_inserter(kpos));
  // std::cerr << "end k_means\n";
      
      std::vector<Point> positions;      
      for (size_t i = 0; i < level.size(); ++i)
	{
	  size_t Nsteps = OT::k_means<K>
	    (k, Density(_density, level._polygons[i]),
	     std::back_inserter(positions));
	}

      double thresh = 1e-2;
      k_means<K> (positions.begin(), positions.end(),
		  _density,
		  std::back_inserter(next_level._positions),
		  std::back_inserter(next_level._masses),
		  std::back_inserter(next_level._polygons),
		  thresh * _density.clamp_radius());

      CGAL::Delaunay_triangulation_2<K> dt;
      dt.insert(level._positions.begin(),
		level._positions.end());
      std::map<Point,size_t> indices;
      for (size_t i = 0; i < level._positions.size(); ++i)
	indices[level._positions[i]] = i;

      double total_mass = 0.0;
      for (size_t i = 0; i < next_level.size(); ++i)
	{
	  total_mass += next_level._masses[i];

	  Point p = dt.nearest_vertex (next_level._positions[i])->point();
	  next_level._parents.push_back(indices[p]);

	  // std::cerr << "(position, mass, parent)[" << i << "] = ("
	  // 	    << next_level._positions[i] << ", "
	  // 	    << next_level._masses[i] << ", "
	  // 	    << next_level._parents[i] << ")\n";
	}
#endif
      for (size_t i = 0; i < next_level.size(); ++i)
	next_level._masses[i] /= total_mass;

    }

    template <class VectorFrom, class VectorTo>
    void prepare_weights (size_t level,
			  const VectorFrom &from,
			  VectorTo &to)
    {
      OT_ASSERT (level < _levels.size());
      const Level &lvl = _levels[level];
      size_t N = lvl.size();

      to.resize(N);
      if ( (level >= 1) && (_levels[level-1].size() == from.size()) )
	{
	  for (size_t i = 0; i < N; ++i)
	    to[i] = from[lvl._parents[i]];
	}
      else
	std::fill(to.begin(), to.end(), 0.0);
    }

    double wasserstein_distance(size_t l, size_t L)
    {
      if (l == L)
	return 0;

      size_t N = level(L).size();
      std::vector<size_t> from (level(L).size());
      
      for (size_t i = 0; i < N; ++i)
	  from[i] = i;

      for (size_t n = L; n > l; n--)
	{
	  for (size_t i = 0; i < N; ++i)
	    from[i] = level(n)._parents[from[i]];
	}

      double w = 0.0;
      for (size_t i = 0; i < N; ++i)
	{
	  size_t j = from[i];

	  double dd = CGAL::squared_distance(level(l)._positions[j],
					     level(L)._positions[i]);
	  w += level(L)._masses[i] * dd;
	}
      return sqrt(w);
    }


    template<class OtherDensity, class RT>
    void compute_weights (OtherDensity other_density,
			  RT &rt,
			  size_t level,
			  std::vector<double> &result,
			  const std::vector<double> &previous =
			  std::vector<double>(),
			  double target_error = 0.05)
    {
      OT_ASSERT (level < _levels.size());
      const Level &lvl = _levels[level];
      size_t N = lvl.size();

      // fill initial weights
      uvector weights (N, 0.0);
      prepare_weights(level, previous, weights);

      typedef typename OT::Objective_function<RT, OtherDensity> Function;
      Function objective (other_density, lvl._positions, lvl._masses);

      OT::bfgs_descent (objective, weights, target_error);

      rt_build (rt, lvl._positions.begin(), lvl._positions.end(), 
		weights.begin(), weights.end());      
      result.clear();
      std::copy(weights.begin(), weights.end(), std::back_inserter(result));
    }
  };
}

#endif
