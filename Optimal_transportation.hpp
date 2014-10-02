#ifndef OPTIMAL_TRANSPORTATION_HPP
#define OPTIMAL_TRANSPORTATION_HPP

void QImage_invert(QImage &image)
{
  image.invertPixels();
  // for (size_t i = 0; i < image.width(); ++i)
  //   for (size_t j = 0; j < image.height(); ++j)
  //     image.setPixel(i,j, 255 - image.pixelIndex(i,j));
}

void QImage_invert(QPixmap &pix)
{
  QImage img = pix.toImage();
  QImage_invert(img);
  pix = QPixmap::fromImage(img);
}


template <class K>
inline
double radius_of_convex_polygon (const typename CGAL::Polygon_2<K> &p,
				 const typename K::Point_2 &c)
{
  double sqr = 0.0;
  for (size_t i = 0; i < p.size(); ++i)
    sqr = std::max(sqr, CGAL::squared_distance(p[i],c));
  return sqrt(sqr);
}

template <class K>
inline typename K::Point_2
linear_interpolation_point (const typename K::Point_2 &a,
			    const typename K::Point_2 &b,
			    double t)
{
  return CGAL::ORIGIN + ((1.0 - t) * (a - CGAL::ORIGIN) +
			 t * (b - CGAL::ORIGIN));
}

inline double
linear_interpolation (double a,
		      double b,
		      double t)
{
  return (1.0 - t) * a + t * b;
}

inline QColor
double_to_QColor(double color)
{
  size_t c = (size_t) std::min(std::max(color, 0.0), 255.0);
  return QColor(c, c, c);
}

template <typename K>
class OptimalTransportationGraphicsItem : public QGraphicsItem
{
  typedef typename CGAL::Polygon_2<K> Polygon;
  typedef typename K::Point_2 Point;
  typedef typename K::Circle_2 Circle;

  std::vector<Point> _origin_centroids, _dest_centroids;
  std::vector<double> _origin_radii, _dest_radii;
  std::vector<double> _masses;

  double _origin_mass_pp, _dest_mass_pp;
  CGAL::Bbox_2 _bb;
  bool _bb_initialized;
  bool _display_gradients;

  std::vector<Point> _cur_centroids;
  std::vector<double> _cur_radii;
  std::vector<QColor> _cur_colors;
  size_t _N;

public:
  OptimalTransportationGraphicsItem() : _display_gradients(true)
  { clear(); }

  void clear()
  {
    _origin_mass_pp = 0;
    _dest_mass_pp = 0;
    _bb = CGAL::Bbox_2(0, 0, 0, 0);
    _bb_initialized = false;
    _N = 0;

    _origin_centroids.clear();
    _dest_centroids.clear();
    _origin_radii.clear();
    _dest_radii.clear();
    _masses.clear();

    clear_render_data();
  }
  
  void clear_render_data()
  {
    _cur_centroids.resize(_N);
    _cur_radii.resize(_N);
    _cur_colors.resize(_N);
  }

  void setTime (double t)
  {
    clear_render_data();
    double cur_mass_pp = linear_interpolation(_origin_mass_pp,
					      _dest_mass_pp, t);

    for (size_t i = 0; i < _N; ++i)
      {
	_cur_centroids[i] =
	  linear_interpolation_point<K> (_origin_centroids[i],
					 _dest_centroids[i], t);
	double radius = linear_interpolation(_origin_radii[i],
					     _dest_radii[i], t);
	_cur_radii[i] = radius;

	double gray = _masses[i] / (M_PI * radius * radius * cur_mass_pp);
	_cur_colors[i] = double_to_QColor(255.0 * gray);
      }
    update();
  }

  // defines the amount of mass one has to spend to fill one pixel to white.
  void setMassPerPixel (double origin_mass_pp,
			double dest_mass_pp)
  {
    _origin_mass_pp = origin_mass_pp;
    _dest_mass_pp = dest_mass_pp;
  }

  void addElement (Polygon &origin_cell,
		   Point &origin_centroid,
		   Polygon &dest_cell,
                   Point &dest_centroid,
                   double mass)
  {
    double origin_radius = radius_of_convex_polygon(origin_cell,
						    origin_centroid),
      dest_radius = radius_of_convex_polygon(dest_cell,
					     dest_centroid);
    _origin_centroids.push_back(origin_centroid);
    _origin_radii.push_back(origin_radius);
    _dest_centroids.push_back(dest_centroid);
    _dest_radii.push_back(dest_radius);
    _masses.push_back(mass);

    if (!_bb_initialized)
      {
	_bb = origin_centroid.bbox();
	_bb_initialized = true;
      }

    _bb = _bb + Circle(origin_centroid, origin_radius).bbox();
    _bb = _bb + Circle(dest_centroid, dest_radius).bbox();
    _N++;
  }

public:
  virtual QRectF boundingRect() const
  {
    return QRectF(_bb.xmin(),
		  _bb.ymin(),
		  _bb.xmax() - _bb.xmin(),
		  _bb.ymax() - _bb.ymin());
  }

  void setDisplayGradients(bool d)
  {
    _display_gradients = d;
    update();
  }

  void paint (QPainter *painter,
	      const QStyleOptionGraphicsItem *option,
	      QWidget *widget)
  {
    size_t N = _cur_centroids.size();

    painter->setPen(QPen(Qt::black));
    painter->setBrush(QBrush(Qt::black));
    painter->drawRect(boundingRect());

    painter->setCompositionMode(QPainter::CompositionMode_Plus);
    for (size_t i = 0; i < N; ++i)
      {
	painter->setPen(QPen(_cur_colors[i]));
	painter->setBrush(QBrush(_cur_colors[i]));
	painter->drawEllipse(_cur_centroids[i].x() - _cur_radii[i],
			     _cur_centroids[i].y() - _cur_radii[i],
			     2 * _cur_radii[i],
			     2 * _cur_radii[i]);
      }

    if (!_display_gradients)
      return;

    painter->setCompositionMode(QPainter::CompositionMode_SourceOver);
    for (size_t i = 0; i < N; ++i)
      {
	painter->setPen(QPen(Qt::red));
	painter->drawLine
	  (_cur_centroids[i].x(), _cur_centroids[i].y(),
	   _cur_centroids[i].x() + .1 * (_dest_centroids[i].x() - 
					  _origin_centroids[i].x()),
	   _cur_centroids[i].y() + .1 * (_dest_centroids[i].y() - 
					  _origin_centroids[i].y()));
      }
  } 
};

#endif
