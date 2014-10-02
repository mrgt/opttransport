#ifndef OT_TOOLS_MISC_HPP
#define OT_TOOLS_MISC_HPP

#include <fstream>
#include <OT/Image_density_2.hpp>
#include <OT/Uniform_density_2.hpp>
#include <OT/serialization.hpp>
#include <OT/Decomposition_2.hpp>

template <class K>
inline
typename OT::Uniform_density_2<K>
make_density (const typename CGAL::Polygon_2<K> &poly)
{
  return typename OT::Uniform_density_2<K>(poly);
}

template <class K>
inline
typename OT::QImage_density_2<K>
make_density (const QImage &img)
{
  return typename OT::QImage_density_2<K>(img);
}

template <class K>
inline
typename CGAL::Polygon_2<K> 
make_rectangle (double xmin, double ymin, double xmax, double ymax,
		K k)
{
  typedef typename K::Point_2 Point;
  typename CGAL::Polygon_2<K> P;
  P.push_back(Point(xmin,ymin));
  P.push_back(Point(xmax,ymin));
  P.push_back(Point(xmax,ymax));
  P.push_back(Point(xmin,ymax));
  return P;
}

  // OT::uvector _solution;
  // const OT::uvector &
  // solution ()
  // {
  //   return _solution;
  // }

  // void save_solution(const std::string &name)
  // {
  //   std::ofstream ofs(name.c_str());
  //   boost::archive::text_oarchive oa(ofs);
  //   oa << _solution;
  // }

  // void load_solution(const std::string &name)
  // {
  //   std::ifstream ifs(name.c_str());
  //   boost::archive::text_iarchive ia(ifs);
  //   ia >> _solution;
  // }

template <class K>
inline
void load (OT::Decomposition_2<K> &decomposition,
	   const std::string &fn)
{
  std::ifstream ifs(fn.c_str());
  boost::archive::text_iarchive ia(ifs);
  ia >> decomposition;
}

template <class K>
inline
void load (OT::Uniform_density_2<K> &density,
	   const std::string &fn)
{
  std::ifstream ifs(fn.c_str());
  double x, y;
  
  typename CGAL::Polygon_2<K> p;
  while (ifs >> x >> y)
    p.push_back (typename K::Point_2 (x,y));

  density = OT::Uniform_density_2<K> (p);
}

template <class K>
inline
void load (OT::QImage_density_2<K> &density,
	   const std::string &fn)
{
  QImage image (fn.c_str());
  density = OT::QImage_density_2<K> (image);
}

#endif

