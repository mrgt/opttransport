#include <fstream>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/point_generators_2.h>

#include <OT/Image_density_2.hpp>
#include <OT/Uniform_density_2.hpp>
#include <OT/Density_decomposition_2.hpp>
#include <OT/integration.hpp>
#include <OT/bfgs_descent.hpp>
#include <OT/lloyd.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_2 Circle_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef CGAL::Polygon_2<K> Polygon_2;

typedef double Weight;
typedef CGAL::Regular_triangulation_euclidean_traits_2<K,Weight>  Gt;
typedef CGAL::Regular_triangulation_2<Gt> Regular;

#include <QtGui>
#include <QString>
#include <QActionGroup>
#include <QFileDialog>
#include <QInputDialog>
#include <QGLWidget>

#include <CGAL/Qt/GraphicsViewCircleInput.h>
#include <CGAL/Qt/PowerdiagramGraphicsItem.h>
#include <CGAL/Qt/utility.h>
#include <CGAL/Qt/DemosMainWindow.h>

#include "ui_Optimal_transportation.h"
#include "Optimal_transportation.hpp"
#include <time.h>


template <class Polygon>
void
polygon_to_QPolygonF(const Polygon &p, QPolygonF &r)
{
  r.clear();
  for (size_t i = 0; i < p.size(); ++i)
    r << QPointF(p[i].x(), p[i].y());
}

Polygon_2
rectangle_to_polygon(const Iso_rectangle_2 &r)
{ 
  Polygon_2 p;
  for (size_t i = 0; i < 4; ++i)
    p.push_back(r.vertex(i));
  return p;
}

double QImage_total_mass(const QImage &image)
{
  double image_mass = 0.0;
  for (size_t i = 0; i < image.width(); ++i)
    for (size_t j = 0; j < image.height(); ++j)
      image_mass += (double)image.pixelIndex(i,j);
  return image_mass;
}

Iso_rectangle_2  QImage_rectangle(const QImage &image)
{
  return Iso_rectangle_2 (0,0, image.width() - 1, image.height() - 1);
}

double QImage_area(const QImage &image)
{
  return image.width() * image.height();
}

double QImage_average_color(const QImage &image)
{
  return QImage_total_mass(image) / QImage_area(image);;
}




class MainWindow :
  public CGAL::Qt::DemosMainWindow,
  public Ui::Optimal_transportation
{
  Q_OBJECT
  
private:
  typedef OT::QImage_density_2<K> Density_to;
  typedef OT::Density_decomposition<K, Density_to> Density_decomposition;
#define IMAGE_FROM
#ifdef IMAGE_FROM
  typedef OT::QImage_density_2<K> Density_from;
#else
  typedef OT::Uniform_density_2<K> Density_from;
#endif

  Regular dt; 
  QGraphicsScene scene;  
  CGAL::Qt::PowerdiagramGraphicsItem<Regular> * vgi;
  OptimalTransportationGraphicsItem<K> *otgi;

  QImage image, image_from;
  Density_to density_to;
  Density_from density_from;
  Density_decomposition decomposition;

  std::vector<double> _weights;
  size_t _level;

public:
  MainWindow();

public slots:
  void on_actionShowPowerdiagram_toggled(bool checked);  
  void on_actionDisplayGradients_toggled(bool checked)
  {
    otgi->setDisplayGradients(checked);
  }

  void on_actionInsertRandomPoints_triggered();
  void on_actionShowColoredPowerDiagram_triggered();
  void on_actionFPAurenhammerStep_triggered();
  void on_actionLoadPoints_triggered();
  void on_actionSavePoints_triggered();
  void on_actionClear_triggered();
  void on_actionRecenter_triggered();
  void on_time_changed(int t)
  {
    otgi->setTime(double(t)/100.0);  
  }

signals:
  void changed();
};

MainWindow::MainWindow()
  : DemosMainWindow()
{
  setupUi(this);

  // render using OpenGL
  graphicsView->setViewport(new QGLWidget);

  QStringList args = QCoreApplication::arguments();

  // Add a GraphicItem for the Powerdiagram diagram
  vgi = new CGAL::Qt::PowerdiagramGraphicsItem<Regular>(&dt);

  QObject::connect(this, SIGNAL(changed()),
		   vgi, SLOT(modelChanged()));

  vgi->setEdgesPen(QPen(Qt::blue, 0, Qt::SolidLine, Qt::RoundCap, Qt::RoundJoin));
  vgi->setVisible(true);
  scene.addItem(vgi);

  QObject::connect(this->actionQuit, SIGNAL(triggered()), 
		   this, SLOT(close()));

  // Check two actions 
  this->actionShowPowerdiagram->setChecked(true);

  OT_ASSERT(args.size() == 3);
  image.load(args[1]);
  QImage_invert(image);
  //image.load("nb.jpg");
  density_to = Density_to(image);
  decomposition = Density_decomposition(density_to);
  decomposition.compute_next_level();
  _level = 1;

  image_from.load(args[2]);
  QImage_invert(image_from);
#ifdef IMAGE_FROM
  density_from = Density_from(image_from);
#else
  density_from = Density_from(OT::QImage_polygon<K>(image));
#endif

  //
  // Setup the scene and the view
  //
  scene.setItemIndexMethod(QGraphicsScene::NoIndex);
  scene.setSceneRect(0, 0, image.width(), image.height());
  //scene.addPixmap(QPixmap::fromImage(image));
  this->graphicsView->setScene(&scene);
  this->graphicsView->setMouseTracking(true);

  // Turn the vertical axis upside down
  this->graphicsView->matrix().scale(1, -1);
                                                      
  // The navigation adds zooming and translation functionality to the
  // QGraphicsView
  this->addNavigation(this->graphicsView);

  this->setupStatusBar();
  this->setupOptionsMenu();
  this->addAboutCGAL();

  otgi = new OptimalTransportationGraphicsItem<K>;
  scene.addItem(otgi);
}

void
MainWindow::on_actionShowPowerdiagram_toggled(bool checked)
{
  vgi->setVisible(checked);
}


void
MainWindow::on_actionClear_triggered()
{
  dt.clear();
  emit(changed());
}

void
MainWindow::on_actionInsertRandomPoints_triggered()
{
  OT::Random_point_in_polygon<K> r (density_to.bounding_poly());
  boost::mt19937 eng;

  const int number_of_points = 
    QInputDialog::getInteger(this, 
                             tr("Number of random points"),
                             tr("Enter number of random points"), 100, 0);

  // wait cursor
  QApplication::setOverrideCursor(Qt::WaitCursor);

  std::vector<Point_2> points;

  points.reserve(number_of_points);
  for(int i = 0; i < number_of_points; ++i)
    points.push_back(r(eng));
  dt.insert(points.begin(), points.end());

  // default cursor
  QApplication::setOverrideCursor(Qt::ArrowCursor);
  emit(changed());
}

void
MainWindow::on_actionShowColoredPowerDiagram_triggered()
{
  const Density_decomposition::Level &lvl = decomposition.level(_level - 1);
  double avg_color_from = QImage_average_color(image_from);
  double avg_color_to = QImage_average_color(image);

  // find correct indices
  std::map<Point_2, size_t> indices;
  for (size_t i = 0; i < lvl.size(); ++i)
    indices[lvl._positions[i]] = i;

  std::vector<Regular::Vertex_handle> vertices;
  OT::rt_copy_vertices_in_order (dt, indices, vertices);

#if 1
  otgi->clear();
  for (size_t i = 0; i < vertices.size(); ++i)  
    {
      Regular::Vertex_handle v = vertices[i];

      double mass = lvl._masses[i];
      CGAL::Polygon_2<K> upoly_from, poly_from;
      std::pair<Point_2, double> c_from;

      OT::tessellate_voronoi_cell(dt, v, std::back_inserter(upoly_from),
				  density_from.clamp_radius());
      OT::clip(density_from.bounding_poly(), upoly_from, poly_from);
      if (!OT::centroid(density_from, poly_from, c_from, true))
	continue;

      CGAL::Polygon_2<K> poly_to = lvl._polygons[i];
      std::pair<Point_2, double> c_to;
      if (!OT::centroid(density_to, poly_to, c_to, true))
	continue;
      
      otgi->addElement(poly_from, c_from.first,
		       poly_to, c_to.first, mass);
    }

  otgi->setMassPerPixel(255.0 / (avg_color_from * QImage_area(image_from)),
			255.0 / (avg_color_to * QImage_area(image)));
  otgi->setTime(0);
  
  emit(changed());
#endif

#if 1
  size_t N = 50;
  for (size_t i = 0; i <= N; ++i)
  {
    QPixmap img (512, 512);
    QPainter painter;
    
    double t = double(i)/double(N);
    std::vector<double> weights = _weights;
    for (size_t j = 0; j < _weights.size(); ++j)
      {
	double w = t * _weights[j];
	weights[j] = w;
      }
    
    Regular rt;
    OT::rt_build (rt, lvl._positions.begin(), lvl._positions.end(), 
		  weights.begin(), weights.end());

    std::map<Point_2, size_t> indices;
    for (size_t j = 0; j < lvl.size(); ++j)
      indices[lvl._positions[j]] = j;
    
    std::vector<Regular::Vertex_handle> vertices;
    OT::rt_copy_vertices_in_order (rt, indices, vertices);
    
    painter.begin(&img);
    for (size_t j = 0; j < vertices.size(); ++j)  
      {
	Regular::Vertex_handle v = vertices[j];
	
	// we spread the mass contained in poly_from on the corresponding
	// polygon in lvl.
	std::pair<Point_2, double> c;
	double mass = lvl._masses[j];
	
	CGAL::Polygon_2<K> upoly_to, poly_to;
	OT::tessellate_voronoi_cell(rt, v, std::back_inserter(upoly_to),
				    std::max(density_from.clamp_radius(),
					     density_to.clamp_radius()));
	OT::clip(density_from.bounding_poly(), upoly_to, poly_to);
	double avg_color = t * avg_color_from + (1-t) * avg_color_to;
	double color =
	  mass * avg_color * (QImage_area(image) / poly_to.area());
	QPolygonF qpoly;
	polygon_to_QPolygonF(poly_to, qpoly);

#if 1
	painter.setPen(QPen(double_to_QColor(color)));
#else
	QPen pen(QColor(0,255,255));
	pen.setWidth(2);
	  
	painter.setPen(pen);
#endif
	painter.setBrush(QBrush(double_to_QColor(color)));
	painter.drawConvexPolygon(qpoly);
      }
    painter.end();

    std::ostringstream ss;
    ss << "images/results/test" << i << ".jpg";
    QImage_invert(img);
    img.save(ss.str().c_str());
  }

  return;
#endif

  for (size_t i = 0; i < vertices.size(); ++i)
    {
      Regular::Vertex_handle v = vertices[i];
      std::pair<Point_2, double> c;
      if (OT::centroid(density_from, dt, v, c))
	{
	  scene.addLine(lvl._positions[i].x(),
			lvl._positions[i].y(), 
			c.first.x(),
			c.first.y());
	}
    }
}

void 
MainWindow::on_actionFPAurenhammerStep_triggered()
{  
  std::vector<double> result;
  decomposition.compute_next_level(5);

  double target_errors[] = {0.005, 0.005, 0.01, 0.05, 0.15};

  //_weights.resize(0);
  decomposition.compute_weights(density_from, dt, _level, result, _weights,
				1e-6);
  _weights = result;
  _level ++;

  std::cerr << "Solved optimal transportation with " << _weights.size() 
	    << " points\n";
  emit(changed());
}

void
MainWindow::on_actionLoadPoints_triggered()
{
  QString fileName = QFileDialog::getOpenFileName(this,
						  tr("Open Points file"),
						  ".");
  if(! fileName.isEmpty()){
    std::ifstream ifs(qPrintable(fileName));

    Gt::Weighted_point_2 p;
    std::vector<Gt::Weighted_point_2> points;
    while(ifs >> p) {
      points.push_back(p);
    }
    dt.insert(points.begin(), points.end());

    actionRecenter->trigger();
    emit(changed());
  }
}


void
MainWindow::on_actionSavePoints_triggered()
{
  QString fileName = QFileDialog::getSaveFileName(this,
						  tr("Save points"),
						  ".");
  if(! fileName.isEmpty()){
    std::ofstream ofs(qPrintable(fileName));
    for(Regular::Finite_vertices_iterator 
          vit = dt.finite_vertices_begin(),
          end = dt.finite_vertices_end();
        vit!= end; ++vit)
    {
      ofs << vit->point() << std::endl;
    }
  }
}


void
MainWindow::on_actionRecenter_triggered()
{
}


#include "Optimal_transportation.moc"

int main(int argc, char **argv)
{
  QApplication app(argc, argv);
  app.setApplicationName("Optimal_transportation");

  // Import resources from libCGALQt4.
  // See http://doc.trolltech.com/4.4/qdir.html#Q_INIT_RESOURCE
  Q_INIT_RESOURCE(File);
  Q_INIT_RESOURCE(Optimal_transportation_2);
  Q_INIT_RESOURCE(Input);
  Q_INIT_RESOURCE(CGAL);

  MainWindow mainWindow;
  mainWindow.show();
  return app.exec();
}
