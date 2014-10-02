#ifndef OT_IMAGE_DENSITY_HPP
#define OT_IMAGE_DENSITY_HPP

#include <OT/Density_base_2.hpp>
#include <OT/intersection.hpp>
#include <OT/rasterization.hpp>
#include <OT/misc.hpp>

#include <QImage>

namespace OT
{

  template<class K, class Image>
  typename CGAL::Polygon_2<K> QImage_polygon(const Image &image)
  {
    double mx = image.width() - 1, my = image.height() - 1;

    typename CGAL::Polygon_2<K> p;
    typedef typename K::Point_2 Point;

    p.push_back(Point(0,0));
    p.push_back(Point(mx,0));
    p.push_back(Point(mx,my));
    p.push_back(Point(0,my));
    return p;
  }

  class QImage_mass_functor
  {
  protected:
    const QImage *_image;
    double _scale;
    double _mass;

  public:
    QImage_mass_functor (const QImage *img, double scale):
      _image(img), _scale(scale), _mass(0.0)
    {}

    inline void
    operator () (int i, int j, double cov = 1.0)
    {
      _mass += cov * _image->pixelIndex(i,j);
    }

    double mass() const
    {
      return _scale * _mass;
    }
  };

  template <class K>
  class QImage_centroid_functor : public QImage_mass_functor
  {
    typedef typename K::Point_2 Point;
    typedef typename K::Vector_2 Vector;

    Vector _centroid;
    size_t _N;

  public:
    QImage_centroid_functor (const QImage *img, double scale):
      QImage_mass_functor(img, scale), _centroid(0.0, 0.0), _N(0)
    {}
    
    inline void
    operator () (int i, int j, double cov = 1.0)
    {    
      int m = cov * _image->pixelIndex(i,j);
      _mass += m;
      _centroid = _centroid + m * Vector(i,j);
      _N++;
    }

    Point centroid() const
    {
      return CGAL::ORIGIN + (1.0/_mass) * _centroid;
    }
    
    bool is_valid()
    {
      return (_mass > 0);
    }
  };

  template <class Function>
  class QImage_integrate_functor
  {
  protected:
    const QImage *_image;
    double _scale;
    double _mass;
    Function _f;

  public:
    QImage_integrate_functor (const QImage *img, double scale, Function f):
      _image(img), _scale(scale), _mass(0.0), _f(f)
    {
      assert(_scale != 0.0);
    }

    inline void
    operator () (int i, int j, double cov = 1.0)
    {
      _mass += cov * _f(i,j) * _image->pixelIndex(i,j);
    }

    double mass() const
    {
      return _scale * _mass;
    }
  };

  typedef std::pair<int, int> Pixel_coordinate;
  typedef std::vector< std::pair<double, Pixel_coordinate> >
      Pixel_cumulative_distribution;

  class Compare_cum_pixel
  {
  public:
    bool operator () (const std::pair<double, Pixel_coordinate> &p, 
		      double r)
    {
      return p.first < r;
    }
  };
  
  class Cumulative_raster
  {
    const QImage &_image;
    Pixel_cumulative_distribution  &_cum;
    double _total;

  public:
    Cumulative_raster (const QImage &image, 
		       Pixel_cumulative_distribution &cum) : 
      _image(image), _cum(cum), _total(0)
    {
      _cum.clear();
    }
    
    void operator () (int x, int y, double cover)
    {
      _total += cover * _image.pixelIndex(x,y);
      _cum.push_back(std::make_pair(_total, std::make_pair(x,y)));
    }

    void normalize()
    {
      if (_total = 0.0)
	return;

      for (size_t i = 0; i < _cum.size(); ++i)
	_cum[i].first /= _total;
    }
  };

  template <class K>
  class QImage_random_point_generator
  {
    typedef typename K::Point_2 Point;
    
    Pixel_cumulative_distribution _cum;
    CGAL::Polygon_2<K> _poly;
    
  public:
    QImage_random_point_generator (const QImage &image, 
				   const CGAL::Polygon_2<K> &poly):
      _poly(poly)
    {
      OT::Cumulative_raster rast(image, _cum);
      rasterize_convex_polygon(poly, rast, true);
      rast.normalize();
    }
    
    QImage_random_point_generator () {}
    
    template <class Engine>
    Point operator () (Engine &engine) const
    {
      boost::uniform_real<> u01(0,1);
      OT_ASSERT(_cum.size() > 0);
      Point result;

      do
	{
	  double r = u01(engine);

	  Pixel_cumulative_distribution::const_iterator it = 
	    std::lower_bound (_cum.begin(), _cum.end(), u01(engine),
			      Compare_cum_pixel());
	  
	  OT_ASSERT(it != _cum.end());
	  size_t i = it - _cum.begin();
	  
	  double dx = u01(engine), dy = u01(engine);
	  const Pixel_coordinate &p = _cum[i].second;
	  result = Point (p.first + dx, p.second + dy);
	} while (_poly.bounded_side(result) != CGAL::ON_BOUNDED_SIDE);

      return result;
    }
  };

  template <class K>
  class QImage_density_2 : public Density_base_2<K>
  {  
  public:
    typedef typename OT::Density_base_2<K> Parent;
    typedef typename CGAL::Polygon_2<K> Polygon;
    typedef typename K::Point_2 Point;
    typedef typename K::Vector_2 Vector;
    typedef typename K::Segment_2 Segment;
    typedef typename K::FT FT;
    typedef typename K::Iso_rectangle_2 BBox;
    typedef QImage_random_point_generator<K> Random_point_generator;
    typedef OT::Random_point_in_polygon<K> Random_point_in_polygon;

  private:
    QImage _image;
    double _scale;
    bool _exact;
    Random_point_generator _rpg;
    Random_point_in_polygon _rpp;

    double compute_scale(bool exact)
    {
      assert (_image.depth() == 8);

      QImage_mass_functor f(&_image, 1.0);
      rasterize_convex_polygon(Parent::bounding_poly(), f, exact);
      double total_mass = f.mass();

      if (total_mass == 0.0)
	return 0.0;
      else
	return 1.0/total_mass;
    }

  public:
    QImage_density_2() :
      _scale(0.0)
    {}

    QImage_density_2(const QImage &image, bool exact = true) :
      Density_base_2<K>(OT::QImage_polygon<K>(image)),
      _image(image),
      _scale(compute_scale(exact)),
      _exact(exact),
      _rpg (image, Parent::bounding_poly()),
      _rpp (Parent::bounding_poly())
    {}

    // clip MUST be clipped against the original density _bounding_poly.
    QImage_density_2 (const QImage_density_2 &orig,
		      const Polygon &clip,
		      bool exact = true):
      Density_base_2<K>(clip),
      _image(orig._image),
      _scale(compute_scale(exact)),
      _exact(exact),
      _rpg (_image, clip),
      _rpp (clip)
    {}

    template <class Engine>
    Point
    random_point_(Engine &eng) const
    {
      if (_scale == 0.0)	
	return _rpg(eng);
      else
	return _rpp(eng);
    }
    
    bool
    centroid_ (const Polygon &poly,
	       std::pair <Point, FT> &p) const
    {
      QImage_centroid_functor<K> f(&_image, _scale);
      if (!OT::rasterize_convex_polygon(poly, f, _exact) ||
	  !f.is_valid())
	return false;

      p = std::make_pair(f.centroid(), f.mass());
      return true;
    }

    double
    mass_ (const Polygon &poly) const
    {
      QImage_mass_functor f(&_image, _scale);
      OT::rasterize_convex_polygon(poly, f, _exact);
      return f.mass();
    }

    double
    mass_ (const Segment &seg) const
    {
      QImage_mass_functor f(&_image, _scale);
      OT::rasterize_segment<K> (seg.source(), seg.target(), f);
      return f.mass();
    }

    template <class Function>
    double
    integrate_ (const Polygon &poly,
		const Function &f) const
    {
      QImage_integrate_functor<Function> ifunc(&_image, _scale, f);
      rasterize_convex_polygon(poly, ifunc, _exact);
      return ifunc.mass();
    }
  };


}

#endif
