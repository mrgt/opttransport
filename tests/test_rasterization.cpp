#include <fstream>
#include <OT/rasterization.hpp>
#include <OT/random.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef CGAL::Polygon_2<K> Polygon_2;

class Image
{
  size_t _width, _height;
  std::vector<size_t> _data;

public:
  Image (size_t width, size_t height):
    _width(width), _height(height),
    _data(width * height, 0)
  {}

  size_t width () const
  {
    return _width;
  }

  size_t height () const
  {
    return _height;
  }

  void fill_random()
  {
    boost::mt19937 rng;
    boost::uniform_int<> r(0,255);
    for (size_t i = 0; i < width(); ++i)
      for (size_t j = 0; j < height(); ++j)
      setPixel(i, j, r (rng));
  }

  size_t pixelIndex (size_t x, size_t y) const
  {
    return _data[x + y *_width];
  }

  void setPixel (size_t x, size_t y, size_t v)
  {
    _data[x + y *_width] = v;
  }
};

template<class K>
typename CGAL::Polygon_2<K> Image_polygon(const Image &image)
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


class Image_raster_functor
{
protected:
  const Image *_image;
  double _result;
  double _scale;

  public:
  Image_raster_functor (const Image *img):
    _image(img),
    _result(0),
    _scale(0)
  {}
  
  inline void
  operator () (int i, int j, double cov = 1.0)
  {
    _result += cov * _image->pixelIndex(i,j);
    _scale += cov;
  }

  inline void
  operator () (double x, double y)
  {
    _result += _image->pixelIndex((int) floor(x), (int) floor(y));
    _scale += 1;
  }
  
  double result() const
  {
    return _result / _scale;
  }
};

class Trivial_raster_functor
{
protected:
  double _result;

  public:
  Trivial_raster_functor ():
    _result(0)
  {}
  
  inline void
  operator () (int i, int j, double cov = 1.0)
  {
    _result += cov;
  }
  
  double result() const
  {
    return _result;
  }
};


double
raster_exact (const Image &image, 
	       const Polygon_2 &p)
{
  Image_raster_functor f(&image);
  OT::rasterize_convex_polygon(p, f, true);
  return f.result();
}

void show_polygon (const Polygon_2 &p)
{
  size_t N = p.size();
  std::cerr << "[ ";
  for (size_t i = 0; i < N; ++i)
    {
      std::cerr << "(" << p[i].x() << ", " << p[i].y() << ")";
      if (i < N-1)
	std::cerr << "; ";
    }
  std::cerr << " ] \n";
}

double
raster_exact_slow (const Image &image, 
		   const Polygon_2 &p)
{
  CGAL::Bbox_2 bb = p.bbox();
  int xmin = (int) bb.xmin();
  int xmax = (int) bb.xmax();
  int ymin = (int) bb.ymin();
  int ymax = (int) bb.ymax();

  double sum, scov;
  for (size_t x = xmin; x <= xmax; ++x)
    for (size_t y = ymin; y <= ymax; ++y)
      {
	double cov = OT::details::pixel_covering(p, x, y);
	sum += cov * image.pixelIndex(x,y);
	scov += cov;
      }
  return sum/scov;
}

bool
test_rasterization_triangle(const Image &image,
			    size_t ntests = 15)
{
  Polygon_2 rect (Image_polygon<K> (image));
  OT::Random_point_in_polygon<K> rr(rect);
  boost::mt19937 eng (time(NULL));

  size_t failure = 0;
  for (size_t k = 0; k < ntests; ++k)
    {
      Polygon_2 g_poly = Polygon_2();
      for (size_t i = 0; i < 3; ++i)
	{
	  Point_2 pt (rr(eng));
	  //g_poly.push_back ( Point_2 (floor(pt.x()), floor(pt.y())) );
	  g_poly.push_back ( Point_2 (pt.x(), pt.y()) );
	}
      if (g_poly.orientation () != CGAL::LEFT_TURN)
	g_poly.reverse_orientation();
      
      double got = raster_exact(image, g_poly),
	exp = raster_exact_slow(image, g_poly);
      std::cerr << "Test " << k << ": " 
		<<  got
		<< " (exp: "
		<< exp << ")\n";

      //OT_DEBUG_SHOW(raster_approx(image, g_poly));
      
      if (fabs(exp - got) >= 1e-4)
	++failure;
    }
  std::cerr << "Convex rasterization: "
	    << (100.0 * failure)/ntests << "% failure ("
	    << failure << "/" << ntests << ")\n";
  return true;
}

double
raster_segment_exact_slow (const Point_2 &a, 
			   const Point_2 &b)
{
  return sqrt(pow(a.x() - b.x(), 2.0) + pow(a.y() - b.y(), 2.0));
}

double
raster_segment_exact (const Point_2 &a, 
		      const Point_2 &b)
{
  Trivial_raster_functor f;
  OT::rasterize_segment(a, b, f);
  return f.result();
}

bool
test_rasterization_segment(const Image &image,
			   size_t ntests = 100)
{
  Polygon_2 rect (Image_polygon<K> (image));
  OT::Random_point_in_polygon<K> rr(rect);
  boost::mt19937 eng (time(NULL));

  for (size_t k = 0; k < ntests; ++k)
    {
      Point_2 a (rr(eng)), b(rr(eng));

      double got = raster_segment_exact(a,b),
	exp = raster_segment_exact_slow(a,b);
      std::cerr << "Segment test " << k << ": " 
		<<  got
		<< " (exp: "
		<< exp << ")\n";

      //OT_DEBUG_SHOW(raster_approx(image, g_poly));
      
      if (fabs(exp - got) >= 1e-4)
	return false;
    }
  return true;  
}

int main(int argc, const char **argv)
{
  Image image (512, 521);
  image.fill_random();

  test_rasterization_triangle (image);
  test_rasterization_segment (image);
  return 0;
}
