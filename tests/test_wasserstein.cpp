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
#include <OT/gradient_descent.hpp>
#include <OT/lloyd.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_2 Point_2;
typedef K::Point_2 Circle_2;
typedef K::Iso_rectangle_2 Iso_rectangle_2;
typedef CGAL::Polygon_2<K> Polygon_2;

typedef double Weight;
typedef CGAL::Regular_triangulation_euclidean_traits_2<K,Weight>  Gt;
typedef CGAL::Regular_triangulation_2<Gt> RT;

std::ofstream cresult ("result.txt");
size_t g_maxiter=300;

#define EXPERIMENT_NUM_ACTIVE

enum 
  {
    LOG_DISTANCE_TO_SOLUTION = 1,
    LOG_NUMBER_ACTIVE = 2
  };

int g_logtype = LOG_NUMBER_ACTIVE;

template <class RT, class Density>
class Objective_function_with_log:
  public OT::Objective_function<RT, Density>
{
public:
  typedef typename OT::Objective_function<RT, Density> Parent;
  typedef typename OT::Objective_function_workspace<Parent> Workspace;
  typedef typename Parent::Bare_point Bare_point;

  OT::uvector _solution;
  mutable size_t _neval;

public:
  Objective_function_with_log(const Density &density,		       
			      const std::vector<Bare_point> &positions,
			      const std::vector<double> &masses) :
    Parent(density, positions, masses)
  {}

  void set_solution(const OT::uvector solution)
  {
    _solution = solution;
    _neval = 0;
  }

  template <class InputIterator>
  bool eval (Workspace &wsp,
	     InputIterator begin,
	     InputIterator end) const
  {
    cresult << _neval;
    if (g_logtype == LOG_DISTANCE_TO_SOLUTION && 
	_solution.size() == end - begin && _neval < g_maxiter)
      {
	size_t N = _solution.size();
	OT::uvector r(N);
	double avg(0);

	for (size_t i = 0; i < N; ++i)
	  {
	    r[i] = *(begin+i) - _solution[i];
	    avg += r[i];
	  }
	avg /= N;

	for (size_t i = 0; i < N; ++i)
	  r[i] = (r[i] - avg)/N;

	cresult << " " << boost::numeric::ublas::norm_2(r) << "\n";
      }
    bool result = Parent::eval(wsp, begin, end);
    if (g_logtype & LOG_NUMBER_ACTIVE)
      {
	cresult << " " << wsp.num_active();
      }
    cresult << "\n";

    _neval++;
    return Parent::eval(wsp, begin, end);
  }
};

class Wasserstein_test
{
public:
  typedef OT::QImage_density_2<K> Density_to;
  typedef OT::Density_decomposition<K, Density_to> Density_decomposition;
  //#define IMAGE_FROM
#ifdef IMAGE_FROM
  typedef OT::QImage_density_2<K> Density_from;
#else
  typedef OT::Uniform_density_2<K> Density_from;
#endif
  RT dt; 

  QImage image, image_from;
  Density_to density_to;
  Density_from density_from;
  Density_decomposition decomposition;
  size_t _level;

public:
  Wasserstein_test() : _level(0) {}

  void init(int argc, char **argv)
  {
    OT_ASSERT(argc == 3);

    //std::cerr << "to = " << argv[1] << "\n";
    //std::cerr << "from = " << argv[2] << "\n";;
    
    image.load(argv[1]);
    density_to = Density_to(image);
    decomposition = Density_decomposition(density_to);
    decomposition.compute_next_level();
    _level = 1;
    
    image_from.load(argv[2]);
#ifdef IMAGE_FROM
    density_from = Density_from(image_from);
#else
    density_from = Density_from(OT::QImage_polygon<K>(image));
#endif
  }

  void decompose(size_t N, size_t k = 5)
  {
    for (size_t i = _level; i <= N; ++i)
      decomposition.compute_next_level(k);

    if (_level != N)
      {
	_level = N;
	std::cerr << "decomposed in " << N << " levels (maximum resolution: "
		  << pow(k, N) << ")\n";
      }
  }

  OT::uvector _solution;
  void wasserstein_2(size_t L, double eps=1e-5)
  {
    std::cerr.precision(4);

    std::vector<double> weights;
    boost::timer t;

    for (size_t i = 1; i <= L; ++i)
      {
	OT_DEBUG_SHOW(i);
	std::vector<double> result;
	decomposition.compute_weights(density_from, dt, i,
				      result, weights, eps);
	weights = result;
      }

    size_t N = weights.size();    
    _solution.resize(N);
    std::copy(weights.begin(), weights.end(), _solution.begin());

    std::cerr << "wasserstein_2 (" 
	      << weights.size() << " points): elapsed time ="
	      << t.elapsed() << "s\n";
  }

  void wasserstein_2_aurenhammer(size_t L, double eps=1e-5)
  {
    boost::timer t;
    size_t N = decomposition.level(L).size();
    OT::uvector weights (N, 0.0);

    Objective_function_with_log<RT, Density_from>
      objective (density_from,
		 decomposition.level(L)._positions,
		 decomposition.level(L)._masses);
    OT::rt_bfgs_descent (objective, weights, eps);

    std::cerr << "wasserstein_2_aurenhammer (" 
	      << N << " points): elapsed time = "
	      << t.elapsed() << "s\n";
  }

  void convex_optimization_comparison(size_t L,
				      size_t maxiter,
				      double eps=1e-8)
  {
    boost::timer t;
    size_t N = decomposition.level(L).size();

    Objective_function_with_log<RT, Density_from>
      objective (density_from,
		 decomposition.level(L)._positions,
		 decomposition.level(L)._masses);

    g_maxiter = maxiter;

#ifdef EXPERIMENT_COMP_CVXOPT
    // gradient descent, constant step
    {
      cresult.close();
      cresult.open("cst-step.d");
      OT::uvector weights (N, 0.0);
      objective.set_solution(_solution);
      OT::gradient_descent
	(objective, OT::Line_search_constant_step (1000.0), weights, eps,
	 maxiter);
    }
#endif

    // gradient descent, strong wolfe
    {
      cresult.close();
      cresult.open("strong-wolfe.d");
      OT::uvector weights (N, 0.0);
      objective.set_solution(_solution);
      OT::rt_aurenhammer_gradient_descent(objective, weights, eps,
					  maxiter);
    }

    // lbfgs, strong wolfe
    {
      cresult.close();
      cresult.open("lbfgs-strong-wolfe.d");
      OT::uvector weights (N, 0.0);
      objective.set_solution(_solution);
      OT::rt_bfgs_descent (objective, weights, eps, 
			   LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE,
			   maxiter);
    }

    // lbfgs, more-thuente
    {
      cresult.close();
      cresult.open("lbfgs-more-thuente.d");
      OT::uvector weights (N, 0.0);
      objective.set_solution(_solution);
      OT::rt_bfgs_descent (objective, weights, eps, 
			   LBFGS_LINESEARCH_MORETHUENTE,
			   maxiter);
    }
  }

};
int main(int argc, char **argv)
{
  size_t N = 6;

  Wasserstein_test test;
  test.init(argc, argv);

  test.decompose(N);
  //test.wasserstein_2(N, 1e-6);
  test.convex_optimization_comparison(N, 200);

  return 0;
}
