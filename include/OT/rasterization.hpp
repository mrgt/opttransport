#ifndef OT_RASTERIZATION_HPP
#define OT_RASTERIZATION_HPP

#include <OT/intersection.hpp>
#include <OT/misc.hpp>

#include <CGAL/Polygon_2.h>
#include <limits>

#define RAST_DEBUG 3

namespace OT
{
  namespace details
  {
    enum
      {
	RAST_NONE,
	RAST_FIRST,
	RAST_ONX,
	RAST_ONY,
	RAST_LAST
      };

    template <class K>
    double pixel_covering(const CGAL::Polygon_2<K> &p,
			  int x, int y)
    {
      typedef typename CGAL::Polygon_2<K> Polygon;
      typedef typename K::Point_2 Point;

      Polygon pixel;
      pixel.push_back(Point(x,y)) ;
      pixel.push_back(Point(x+1,y)) ;
      pixel.push_back(Point(x+1,y+1)) ;
      pixel.push_back(Point(x,y+1)) ;
    
      Polygon clipped;
      OT::convex_intersection (p, pixel, clipped);
      return fabs(clipped.area());
    }
  
    class Rasterization_data
    {
    private:
      size_t _xmin, _ymin;
      size_t _width, _height;
      std::vector<double> _data;
#if (RAST_DEBUG >= 3)
      //std::vector< std::pair<int, int> > _types;
#endif
      std::vector<int> _left;
      std::vector<int> _right;

    public:
      Rasterization_data(int xmin, int ymin, int xmax, int ymax) :
	_xmin(xmin), _ymin(ymin),
	_width(xmax - xmin + 1),
	_height(ymax - ymin + 1),
	_data(_width*_height, 1.0),
#if (RAST_DEBUG >= 4)
	_types(_width*_height,
	       std::make_pair((int) RAST_NONE, (int)RAST_NONE)),
#endif
	_left(_height, xmax),
	_right(_height, xmin)
      {}
    
      double &
      cover (int x, int y) 
      {
	OT_ASSERT(x >= _xmin);
	OT_ASSERT(x < _xmin + _width);
	OT_ASSERT(y >= _ymin);
	OT_ASSERT(y < _ymin + _height);

	return _data[ (y - _ymin) * _width + (x - _xmin) ];
      }
      
      int &
      left (int y)
      {
	OT_ASSERT(y >= _ymin);
	OT_ASSERT(y < _ymin + _height);
	return _left[y - _ymin];
      }


      int &
      right (int y)
      {
	OT_ASSERT(y >= _ymin);
	OT_ASSERT(y < _ymin + _height);
	return _right[y - _ymin];
      }

#if (RAST_DEBUG >= 4)
      std::pair<int, int> &
      type(int x, int y)
      {
	return _types[(y - _ymin) * _width + (x - _xmin)];
      }
#endif

      void 
      mark (int x, int y, int tx, int ty)
      {
	left(y) = std::min(left(y), x);
	right(y) = std::max(right(y), x);

#if (RAST_DEBUG >= 4)
	type(x,y) = std::make_pair(tx, ty);
#endif
      }

      size_t width() 
      {
	return _width;
      }

      size_t height() 
      {
	return _height;
      }
    };

    template <int SX>
    int next (double x)
    {
      if (SX > 0)
	{
	  double cx = ceil(x);
	  return (cx == x) ? (x + 1) : cx;
	}
      else
	{
	  double fx = floor(x);
	  return (fx == x) ? (x - 1) : fx;
	}
    }

    template <int SX>
    inline double prev (double x)
    {
      if (SX > 0)
	return floor(x);
      else
	return ceil(x);
    }

    template <int SX>
    inline bool before (double x, double nx)
    {
      return (SX > 0) ? (x < nx) : (nx < x);
    }

    const char * get_type_name(int type)
    {
      const char *T[] = {"none", "first", "on_x", "on_y", "last"};
      return T[type];
    }

    class Rasterization_functor
    {
    public:
      double _x0, _y0;
      int _type0;

      bool _has0;
      bool _exact;
      Rasterization_data &_data;

    public:
      Rasterization_functor (Rasterization_data &data,
			     bool exact = true) :
	_data(data),
	_has0(false),
	_exact(exact)
      {}

      template <int SX, int SY>
      void handle (int type, double x1, double y1)
      {
	if (!_has0)
	  {
	    _x0 = x1, _y0 = y1;
	    _type0 = type;
	    _has0 = true;
	    return;
	  }

	int px = (int) floor( (SX > 0) ? _x0 : x1);
	int py = (int) floor( (SY > 0) ? _y0 : y1);
	_data.mark(px, py, _type0, type);

#if (RAST_DEBUG >= 4)
	std::cerr << "Rasterization_functor::handle()\n";
	std::cerr << " a = (" << _x0 << ", " << _y0 << ")\n";
	std::cerr << " b = (" << x1 << ", " << y1 << ")\n";
	std::cerr << " p = (" << px << ", " << py << ")\n";
	std::cerr << " t = ("
		  << OT::details::get_type_name(_type0) << ", "
		  << OT::details::get_type_name(type) << ")\n";
#endif

	if (!_exact)
	  return;

	double area = 0.0;
	if (_type0 == RAST_ONX && type == RAST_ONX)
	  {
	    //   _______
	    //  |       |- yb
	    //  |    __/|          area =  (1/2) (ya+yb)
	    //  | __/...|          
	    //  |/......|- ya
	    //  |_______|

	    double ybottom = py; // floor(_y0); //prev<SY> 
	    double ya = fabs(_y0 - ybottom);
	    double yb = fabs(y1 - ybottom);
	    area = (ya + yb) / 2.0;

	    if (SX < 0)
	      area = 1.0 - area;
	  }
	else if (_type0 == RAST_ONY && type == RAST_ONY)
	  { 
	    //      xb
	    //   ___|___
	    //  |   /...| 
	    //  |   |...|  area = 1.0 - (1/2) (xa+xb)
	    //  |  /....|          
	    //  |  |....|
	    //  |_/_____|
	    //    xa
	    double xleft = px;
	    double xa = fabs(_x0 - xleft);
	    double xb = fabs(x1 - xleft);
	    area = (xb + xa) / 2.0;

	    if (SY > 0)
	      area = 1.0 - area;
	  }
	else if (_type0 == RAST_ONY && type == RAST_ONX)
	  {
	    //   _______
	    //  |       | 
	    //  |      _|- yb
	    //  |     /.|         area = [(1-xa) * yb]/2
	    //  |    /..|
	    //  |___/___|
	    //    xa
	    double xa = fabs(_x0 - prev<SX> (_x0));
	    double yb = fabs(y1 - prev<SY> (y1));
	    area = (1.0 - xa) * yb / 2.0;

	    if ( SX * SY < 0)
	      area = 1.0 - area;

	  }
	else if (_type0 == RAST_ONX && type == RAST_ONY)
	  {
	    //      xb
	    //     ___|___
	    //    |  /....| 
	    //    | /.....|
	    // ya |/......|         area = 1 - [(1-ya) * xb]/2
	    //    |.......|
	    //    |_______|
	    double ya = fabs(_y0 - prev<SY> (_y0));
	    double xb = fabs(x1 - prev<SX> (x1));
	    area = 1.0 - (1.0 - ya) * xb / 2.0;

	    if ( SX * SY < 0)
	      area = 1.0 - area;
	  }
	_data.cover(px, py) -= area;

	_x0 = x1, _y0 = y1;
	_type0 = type;
      }
    };


    template <int SX, int SY, class Functor>
    void draw_segment (double x, double y, double x1, double y1,
		       Functor &f)
    {
      double deltax (x1 - x), deltay (y1 - y);
      double length (sqrt(deltax * deltax + deltay * deltay));
      double sx = deltax/length, sy = deltay/length;
      double mx = fabs(1.0/sx), my = fabs(1.0/sy);
  
      int ny = next<SY>(y), nx = next<SX>(x);
      double tx = mx * fabs(nx - x);
      double ty = my * fabs(ny - y);

      f.template handle<SX,SY>(RAST_FIRST, x, y);
      while (1)
	{
	  while (ty > tx)
	    {
	      if (!before<SX>(nx, x1))	    
		break;

	      x = nx;
	      y += tx * sy;
	      nx += SX;
	      ty -= tx;
	      tx = mx;
	  
	      f.template handle<SX,SY>(RAST_ONX, x, y);
	    }
      
	  if (!before<SY>(ny, y1))
	    break;
	  
	  y = ny;
	  x += ty * sx;
	  ny += SY;
	  tx -= ty;
	  ty = my;
	  f.template handle<SX,SY>(RAST_ONY, x, y);
	}
      f.template handle<SX,SY>(RAST_LAST, x1, y1);
    }

    template <class Functor>
    void draw_segment (double x, double y, double x1, double y1,
		       Functor &f)
    {
      if (y1 > y)
	{
	  if (x1 > x)
	    draw_segment<1,1> (x, y, x1, y1, f);
	  else
	    draw_segment<-1,1> (x, y, x1, y1, f);
	}
      else
	{
	  if (x1 > x)
	    draw_segment<1,-1> (x, y, x1, y1, f);
	  else
	    draw_segment<-1,-1> (x, y, x1, y1, f);
	}
    }

    template <class DrawFunctor>
    class Rasterization_segment_functor
    {
      DrawFunctor &_draw;
      double _x0, _y0;
      bool _has0;
      
    public:
      Rasterization_segment_functor(DrawFunctor &draw): 
	_draw(draw), _x0(0), _y0(0), _has0(false)
      {}
      
      template <int SX, int SY>
      void handle (int type, double x1, double y1)    
      {
	if (!_has0)
	  {
	    _x0 = x1, _y0 = y1;
	    _has0 = true;
	    return;
	  }
	
	int px = (int) floor( (SX > 0) ? _x0 : x1);
	int py = (int) floor( (SY > 0) ? _y0 : y1);
	double dx = x1 - _x0, dy = y1 - _y0;
	double len = sqrt(dx*dx + dy * dy);
	
	_draw(px, py, len);
	_x0 = x1, _y0 = y1;
      }
    };    
  }

  template <class Functor>
  bool
  rasterize_segment(double xa, double ya, double xb, double yb,
		    Functor &drawSmallSeg)
  {
    typedef typename OT::details::Rasterization_segment_functor<Functor> RSF;

    RSF functor (drawSmallSeg);
    OT::details::draw_segment(xa, ya, xb, yb, functor);
    return true;
  }

  template <class K, class Functor>
  bool rasterize_segment (const typename CGAL::Point_2<K> &a,
			  const typename CGAL::Point_2<K> &b,
			  Functor &drawSmallSeg)
  {
    return rasterize_segment (a.x(), a.y(), b.x(), b.y(), drawSmallSeg);
  }

  template <class K, class Functor>
  bool rasterize_convex_polygon (const typename CGAL::Polygon_2<K> &upoly,
				 Functor &drawPixel,
				 bool exact = true)
  {
    typedef typename K::Point_2 Point;

    // FIXME: perturb polygon to avoid integer coordinates. ugly.
    size_t N = upoly.size();

    if (N == 0)
	return false;

    typename CGAL::Polygon_2<K> poly;
    for(size_t i = 0; i < N; i++)
      {
	double x = upoly[i].x(), y = upoly[i].y();
	if (fabs(round(x) - x) < 1e-8)
	  x += 1e-6;
	if (fabs(round(y) - y) < 1e-8)
	  y += 1e-6;
	poly.push_back(Point(x,y));
      }

    CGAL::Bbox_2 bb = poly.bbox();
    int xmin = bb.xmin(), xmax = bb.xmax();
    int ymin = bb.ymin(), ymax = bb.ymax();

    OT::details::Rasterization_data data(xmin, ymin, xmax, ymax);
    for(size_t i = 0; i < N; i++)
      {
	size_t j = (i+1)%N;
	OT::details::Rasterization_functor rf(data, exact);
	OT::details::draw_segment (poly[i].x(), poly[i].y(),
				   poly[j].x(), poly[j].y(), rf);
      }
    
    for(unsigned int i = 0; i < N; i++)
      {
	int x = (int) floor (poly[i].x());
	int y = (int) floor (poly[i].y());
	data.cover(x,y) = OT::details::pixel_covering(poly, x, y);
#if (RAST_DEBUG >= 4)
	std::cerr << "rasterize_convex_polygon: set pixel (" 
		  << x << ", " << y << ") to " << data.cover(x,y) << "\n";
#endif
      }
    
    for(int y = ymin; y <= ymax; y++)
      {
	const int left = data.left(y), right = data.right(y);
	double *p_cover = &data.cover(data.left(y), y);
	
	for(int x = left; x <= right; x++)
	  {
#if (RAST_DEBUG >= 4)
	    double real = details::pixel_covering(poly, x, y);
	    if ( fabs(*p_cover - real) > 1e-3)
	      {
		std::cerr
		  << __FUNCTION__ << ": "
		  << " pixel (" << x << ", " << y <<  ")"
#if (RAST_DEBUG >= 3)
		  << " with type ("
		  << details::get_type_name(data.type(x,y).first) << ", "
		  << details::get_type_name(data.type(x,y).second) << ")"
#endif
		  << ": \n";
		std::cerr 
		  << "\twanted " << real << " got " << (*p_cover) << "\n";
	      }
#endif
	    drawPixel (x, y, std::max(0.0, *(p_cover++)));
	  }
      }
    
    return true;
  }  
}


#endif
