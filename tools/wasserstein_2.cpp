#include <fstream>
#include <boost/timer.hpp>

// CGAL headers
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/point_generators_2.h>

#include "misc.hpp"
#include <OT/integration.hpp>
#include <OT/bfgs_descent.hpp>
#include <OT/objective_function.hpp>
#include <OT/lloyd.hpp>

#include <boost/filesystem.hpp>
#include <boost/program_options.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K> Polygon;
typedef CGAL::Point_2<K> Point;
typedef CGAL::Regular_triangulation_euclidean_traits_2<K,double>  Gt;
typedef CGAL::Regular_triangulation_2<Gt> RT;

enum 
  {
    LOG_ITERATION =            (1<<1),
    LOG_TIME =                 (1<<2),
    LOG_DISTANCE_TO_SOLUTION = (1<<3),
    LOG_NUMBER_ACTIVE =        (1<<4),
    LOG_WASSERSTEIN_ERROR =    (1<<5),
    LOG_WASSERSTEIN =          (1<<6),
    LOG_WASSERSTEIN_MIN =      (1<<7),
    LOG_WASSERSTEIN_MAX =      (1<<8),
    LOG_FUNCTION_VALUE =       (1<<9),
    LOG_NORM_INF_GRADIENT =    (1<<10),
    LOG_NORM_1_GRADIENT =      (1<<11),
  };

struct Wasserstein_log : public OT::Display_value_gradient
{  
  typedef OT::Display_value_gradient Parent;

  int _log_type;
  size_t _max_iter;
  std::ostream *_log;
  boost::timer _timer;
  OT::uvector _solution;
  size_t _neval;
  double _quantization_error;
  double _wmin, _wmax;

public:
  Wasserstein_log(int log_type = 0) :
    _log(0),
    _log_type(log_type), 
    _max_iter(1e6),
    _neval(0),
    _wmin(0), _wmax(1e6)
  {}

  void set_log_file (std::ostream *os)
  {
    _log = os;
  }

  void set_max_iter(size_t m)
  {
    _max_iter = m;
  }

  void set_solution(const OT::uvector solution)
  {
    _solution = solution;
    _neval = 0;
  }

  bool need (int type)
  {
    return (_log_type & type) && (_neval < _max_iter);
  }

  void log (int type, double value)
  {
    if (need (type) && _log != NULL)
      {
	(*_log) << value << " ";
      }
  }

  void set_wasserstein_bounds (double wmin, 
			       double wmax)
  {
    _wmin = std::max (_wmin, wmin);
    _wmax = std::min (_wmax, wmax);
    log (LOG_WASSERSTEIN_MIN, _wmin);
    log (LOG_WASSERSTEIN_MAX, _wmax);
  }

  template <class Function, class InputIterator, class Workspace>
  void
  operator () (const Function &function,
	       InputIterator begin, InputIterator end,
	       const Workspace &wsp) 
  {
    Parent::operator () (function, begin, end, wsp);

    log(LOG_ITERATION,_neval);
    log(LOG_TIME, _timer.elapsed());
    log(LOG_FUNCTION_VALUE, wsp.value());
    
    if (need(LOG_NORM_INF_GRADIENT))
      {
	log (LOG_NORM_INF_GRADIENT, 
	     OT::norm_inf (wsp.gradient().begin(),
			   wsp.gradient().end()));
      }

    if (need(LOG_NORM_1_GRADIENT))
      {
	log (LOG_NORM_1_GRADIENT, 
	     OT::norm_1 (wsp.gradient().begin(),
			 wsp.gradient().end()));
      }

    if ( (need(LOG_DISTANCE_TO_SOLUTION) && 
	  _solution.size() == end - begin) )
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

	log(LOG_DISTANCE_TO_SOLUTION,
	    boost::numeric::ublas::norm_2(r));
      }
    else
      {
	log(LOG_DISTANCE_TO_SOLUTION, 1000);
      }
    double w = wsp.wasserstein();
    double werr = (function.quantization_error() + wsp.wasserstein_error());
    log(LOG_NUMBER_ACTIVE, wsp.num_active());
    log(LOG_WASSERSTEIN_ERROR, werr);
    log(LOG_WASSERSTEIN, w);
    set_wasserstein_bounds (w - werr, w + werr);

    if (_log) (*_log) << "\n";
    _neval++;
  }
};

class Stop_wasserstein_error: public OT::Stop_with_threshold
{
  mutable double _wmin, _wmax;
  bool _absolute;

 public:
  Stop_wasserstein_error(double eps, bool absolute) : 
    OT::Stop_with_threshold (eps),
    _wmin(0), _wmax(1e6), _absolute (absolute)
  {}

  template <class Function, class Workspace>
  bool
  operator () (const Function &f,
	       const Workspace &wsp) const
  {
    double werr = wsp.wasserstein_error() + f.quantization_error();
    double w = wsp.wasserstein();

    //OT_DEBUG_SHOW(wsp.wasserstein_error());
    //OT_DEBUG_SHOW(f.quantization_error());
    //OT_DEBUG_SHOW(wsp.wasserstein());

    _wmin = std::max (w - werr, _wmin);
    _wmax = std::min (w + werr, _wmax);
    double whalf = (_wmin + _wmax)/2.0;

    //OT_DEBUG_SHOW(OT::norm_inf(wsp.gradient().begin(), wsp.gradient().end()));
    double eps = _absolute ? threshold() : threshold() * whalf;

    if ( fabs (_wmax - whalf) < eps)
      {
	std::cerr << "W = " << whalf << " +/- " << fabs (_wmax - whalf) << "\n";
	return true;
      }

    if (f.quantization_error() >= eps &&
	OT::norm_inf(wsp.gradient().begin(), wsp.gradient().end()) <= 1e-5)
      {
	std::cerr << "change level!\n";
	return true;
      }
       
    return false;
  }
  
};

template <class Density, class K, class StopCriterion>
OT::uvector
wasserstein_2 (const Density &density,
	       const OT::Decomposition_2<K> &decomposition,
	       const StopCriterion &stop,
	       Wasserstein_log &g,
	       bool multiresolution = true,
	       const OT::uvector &sol = OT::uvector())
{
  g._timer.restart();
  
  OT::uvector weights; 
  size_t L =((multiresolution) 
	     ? decomposition.number_of_levels() - 1
	     : 0);

  for (int i = L; i >= 0; --i)
    {
      OT::uvector new_weights;
      decomposition.prepare_weights(i, new_weights, i+1, weights);
           
      OT::Objective_function<RT, Density>
	objective (density,
		   decomposition.level(i).positions(),
		   decomposition.level(i).masses());
      
      objective.set_quantization_error
	(decomposition.wasserstein_distance(0,i));
      
      OT::bfgs_descent (objective,
			stop,
			g, new_weights);
      weights = new_weights;
    }
  
  std::cerr << "wasserstein_2 (" 
	    << weights.size() << " points): elapsed time ="
	    << g._timer.elapsed() << "s\n";
  
  return weights;
}

template <class Density, class K>
OT::uvector
wasserstein_2 (const Density &density,
	       const OT::Decomposition_2<K> &decomposition,
	       std::string stop_criterion,
	       double eps,
	       Wasserstein_log &g,
	       bool multiresolution = true,
	       const OT::uvector &sol = OT::uvector())
{
  if (stop_criterion == "norm-inf-gradient")
    return ::wasserstein_2 (density, decomposition, OT::Stop_norm_inf_gradient(eps), g, multiresolution, sol);
  else if (stop_criterion == "wasserstein-error")
    return ::wasserstein_2 (density, decomposition, Stop_wasserstein_error(eps, false), g, multiresolution, sol);
  else if (stop_criterion == "wasserstein-abs-error")
    return ::wasserstein_2 (density, decomposition, Stop_wasserstein_error(eps, true), g, multiresolution, sol);
}	   

namespace po = boost::program_options;

int main(int argc, char **argv)
{
  double eps;
  std::string stop_criterion, source, target;

  po::options_description desc("Options");
  desc.add_options()
    ("help", "produce help message")
    ("stop", po::value<std::string>(&stop_criterion)->default_value("norm-inf-gradient"),
     "stop criterion used in descent")
    ("eps", po::value<double>(&eps)->default_value(1e-6), "threshold for stop criterion")
    ("source", po::value< std::string >(&source), "source measure with density")
    ("target", po::value< std::string >(&target), "target measure decomposition")
    ("single", "do not use multiscale algorithm")
    ;

  struct Log_option_desc {int value; const char *name, *desc;};
  Log_option_desc log_options [] = 
    {
      { LOG_ITERATION, "log-iteration", "log the number of iterations" },
      { LOG_TIME, "log-time", "log the elapsed time" },
      { LOG_DISTANCE_TO_SOLUTION, "log-distance-to-solution", "log the distance to the solution" },
      { LOG_NUMBER_ACTIVE, "log-number-actives", "log the number of active sites" },
      { LOG_WASSERSTEIN_ERROR, "log-wasserstein-error", "log the wasserstein error" },
      { LOG_WASSERSTEIN, "log-wasserstein", "log the estimated wasserstein distance" },
      { LOG_WASSERSTEIN_MIN, "log-wasserstein-min", "log the best lower bound on the wasserstein distance" },
      { LOG_WASSERSTEIN_MAX, "log-wasserstein-max", "log the best upper bound on the wasserstein distance" },
      { LOG_WASSERSTEIN_MAX|LOG_WASSERSTEIN_MIN, 
	"log-wasserstein-minmax", "log the best lower and upper bound on the wasserstein distance" },
      { LOG_FUNCTION_VALUE, "log-function-value", "log the value of the optimized function" },
      { LOG_NORM_INF_GRADIENT, "log-norm-inf-gradient", "log the INF-norm of the gradient" },
      { LOG_NORM_1_GRADIENT, "log-norm-inf-gradient", "log the 1-norm of the gradient" }
    };
  size_t num_log_options = sizeof(log_options)/sizeof(log_options[0]);
  for (size_t i = 0; i < num_log_options; ++i)
    desc.add_options() (log_options[i].name, log_options[i].desc);

  po::positional_options_description p;
  p.add("source", 1);
  p.add("target", 1);

  po::variables_map args;
  po::store(po::command_line_parser(argc, argv).options(desc)
	    .style (po::command_line_style::default_style |
		    po::command_line_style::allow_long_disguise)
	    .positional(p).run(), args);
  po::notify(args);

  if (args.count("help"))
    {
      std::cout << desc;
      return 1;
    }

  if (!args.count ("source") || !args.count ("target"))
    {
      std::cout << argv[0] << ": missing source or target measure\n";
      return 1;
    }

  bool multiresolution = true;
  if (args.count ("single"))
    multiresolution = false;

  std::cerr << "wasserstein_2 between " << source << " and " << target << "\n";

  int logopt = 0;
  for (size_t i = 0; i < num_log_options; ++i)
    {
      if (args.count (log_options[i].name))
	logopt |= log_options[i].value;
    }
  if (logopt == 0)
    logopt = LOG_ITERATION | LOG_FUNCTION_VALUE | LOG_NORM_INF_GRADIENT;
	

  OT::Decomposition_2<K> decomposition;

  load (decomposition, target);

  Wasserstein_log log (logopt);
  log.set_log_file(&std::cout);

  if (boost::filesystem::extension(source) == ".poly")
    {
      OT::Uniform_density_2<K> density;
      load (density, source);
      wasserstein_2(density, decomposition,
		    stop_criterion, eps, log, multiresolution);
    }
  else
    {
      OT::QImage_density_2<K> density;
      load (density, source);
      wasserstein_2(density, decomposition, stop_criterion, eps, log, multiresolution);
    }
 

  return 0;
}
