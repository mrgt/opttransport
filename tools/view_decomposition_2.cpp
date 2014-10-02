#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Bbox_2.h>
#include <OT/Decomposition_2.hpp>
#include <OT/Image_density_2.hpp>
#include <OT/serialization.hpp>

#include <QImage>
#include <QPainter>
#include <QGraphicsItem>
#include <QGraphicsView>
#include <QGraphicsScene>
#include <QApplication>

template <class Polygon>
void
polygon_to_QPolygonF(const Polygon &p, QPolygonF &r)
{
  r.clear();
  for (size_t i = 0; i < p.size(); ++i)
    r << QPointF(p[i].x(), p[i].y());
}

inline QColor
double_to_QColor(double color)
{
  size_t c = (size_t) std::min(std::max(color, 0.0), 255.0);
  return QColor(c, c, c);
}

namespace OT
{
  template <typename K>
  class DiscreteMeasureGraphicsItem : public QGraphicsItem
  {
    typedef typename CGAL::Polygon_2<K> Polygon;
    typedef typename CGAL::Point_2<K> Point;
    typedef typename OT::Discrete_measure_2<K> Measure;

    std::vector<QColor> _colors;
    std::vector<QPolygonF> _polygons;
    std::vector<QPointF> _positions;
    CGAL::Bbox_2 _bb;

  public:
    DiscreteMeasureGraphicsItem()
    {}

    void 
    set_data (const Measure &mu,
	      const Polygon &clip,
	      double average_color)
    {
      typedef typename CGAL::Delaunay_triangulation_2<K> DT;
      typedef typename DT::Finite_vertices_iterator Iterator;
      
      std::map<Point, double> masses;      
      for (size_t i = 0; i < mu.size(); ++i)
	{
	  masses [mu.position(i)] = mu.mass(i);
	  _positions.push_back(QPointF(mu.position(i).x(),
				       mu.position(i).y()));
	}

      DT dt;
      dt.insert(mu.positions().begin(),
		mu.positions().end());

      double total_area = clip.area();

      _colors.clear();
      _polygons.clear();

      size_t N = 0;
      for (Iterator it = dt.finite_vertices_begin();
	   it != dt.finite_vertices_end(); ++it, ++N)
	{
	  Polygon upoly, poly;
	  OT::tessellate_voronoi_cell(dt, it, std::back_inserter(upoly),
				      2.0 * mu.clamp_radius());
	  OT::clip (clip, upoly, poly);

	  double mass = masses[it->point()];
	  double area = poly.area();
	  double gray = mass * average_color * (total_area/area);

	  QPolygonF qpoly;
	  polygon_to_QPolygonF(poly, qpoly);

	  _colors.push_back (double_to_QColor(gray));
	  _polygons.push_back (qpoly);
	}

      std::cerr << N << " Voronoi cells\n";

      _bb = clip.bbox();
    }

public:
  virtual QRectF boundingRect() const
  {
    return QRectF(_bb.xmin(),
		  _bb.ymin(),
		  _bb.xmax() - _bb.xmin(),
		  _bb.ymax() - _bb.ymin());
  }

  void paint (QPainter *painter,
	      const QStyleOptionGraphicsItem *option,
	      QWidget *widget)
  {
    size_t N = _polygons.size();

    for (size_t i = 0; i < N; ++i)
      {
	painter->setPen(QPen(_colors[i]));
	painter->setBrush(QBrush(_colors[i]));
	painter->drawConvexPolygon(_polygons[i]);
	painter->setPen(QPen(Qt::red));
	painter->drawEllipse(_positions[i], 1, 1);
      }
  } 
};

}

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Point_2<K> Point;

int main(int argc, char **argv)
{
  QApplication app (argc, argv);

  OT_ASSERT(argc == 3);

  OT::Decomposition_2<K> dec;

  std::ifstream ifs(argv[1]);
  boost::archive::text_iarchive ia(ifs);
  ia >> dec;

  //////////////////////////////////////////
  OT::DiscreteMeasureGraphicsItem<K> gi;

  Polygon clip;
  size_t W = 512, H = 512;
  clip.push_back(Point(0,0));
  clip.push_back(Point(W,0));
  clip.push_back(Point(W,H));
  clip.push_back(Point(0,H));
  OT_ASSERT(clip.is_simple());

  gi.set_data(dec.level(strtod(argv[2], NULL)),
	      clip, 128.0);
  
  QGraphicsScene scene;
  scene.addItem (&gi);

  QGraphicsView view(&scene);
  view.setFixedSize(W+10, H+10);
  view.show();

  return app.exec();
}

