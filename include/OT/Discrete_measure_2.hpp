#ifndef OT_DISCRETE_MEASURE_2_HPP
#define OT_DISCRETE_MEASURE_2_HPP

#include <vector>
#include <algorithm>

#include <boost/random/uniform_real.hpp>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Bbox_2.h>

#include <OT/misc.hpp>

namespace OT
{
  class Discrete_random_point_generator_2
  {
    std::vector<double> _cum;
    
  public:
    Discrete_random_point_generator_2 (const std::vector<double> &masses)
    {
      double cur = 0.0;
      for (size_t i = 0; i < masses.size(); ++i)
	{
	  cur += masses[i];
	  _cum.push_back(cur);
	}
    }
    
    Discrete_random_point_generator_2 () {}
    
    template <class Engine>
    size_t operator () (Engine &engine) const
    {
      OT_ASSERT(_cum.size() > 0);

      boost::uniform_real<> u01(0,1);      
      std::vector<double>::const_iterator it = 
	std::lower_bound (_cum.begin(), _cum.end(), u01(engine));

      OT_ASSERT (it != _cum.end());

      return (it - _cum.begin());
    }
  };

  template <class K>
  class Discrete_measure_2
  {
  public:
    typedef typename K::Point_2 Point;
    
  private:
    std::vector<Point> _positions;
    std::vector<double> _masses;    
    Discrete_random_point_generator_2 _rpg;
    CGAL::Bbox_2 _bbox;

  public:
    template<class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & _positions;
        ar & _masses;

	if (Archive::is_loading::value)
	  prepare();
    }

  public:
    Discrete_measure_2() {}

    Discrete_measure_2(const std::vector<double> &masses,
		       const std::vector<Point> &positions) :
      _masses(masses),
      _positions(positions)
    { prepare(); }

    void normalize()
    {
      double total_mass = 0.0;
      for (size_t i = 0; i < _masses.size(); ++i)
	total_mass += _masses[i];

      for (size_t i = 0; i < _masses.size(); ++i)
	_masses[i] /= total_mass;
    }

    void prepare()
    {
      assert(_masses.size() == _positions.size());
      normalize();

      _rpg = Discrete_random_point_generator_2 (_masses);

      OT_ASSERT(_positions.size() > 0);
      _bbox = CGAL::Bbox_2 (_positions[0].x(),
			    _positions[0].y(),
			    _positions[0].x(),
			    _positions[0].y());

      for (size_t i = 1; i < _positions.size(); ++i)
	{
	  _bbox = _bbox + CGAL::Bbox_2 (_positions[i].x(),
					_positions[i].y(),
					_positions[i].x(),
					_positions[i].y());
	}
    }

    void output(std::ostream &os) const
    {
      //assert(_masses.size() == _positions.size());
      os << _masses.size() << "\n";
      for(size_t i = 0; i < size(); ++i)
	os << _positions[i] << " " << _masses[i] << "\n";
    }

    double clamp_radius() const
    {
      return sqrt (pow(_bbox.xmax() - _bbox.xmin(), 2.0) + 
		   pow(_bbox.ymax() - _bbox.ymin(), 2.0));
    }

    template <class Engine>
    Point
    random_point_(Engine &eng) const
    {
      return _positions[_rpg(eng)];
    }

    size_t size() const
    {
      return _positions.size();
    }

    const std::vector<Point> &positions() const
    {
      return _positions;
    }

    const std::vector<double> &masses() const
    {
      return _masses;
    }

    const Point &position(size_t i) const
    {
      return _positions[i];
    }

    double mass(size_t i) const
    {
      return _masses[i];
    }

    void push_back (double m, const Point &p)
    {
      _masses.push_back(m);
      _positions.push_back(p);
    }
  };
}

#endif
