#ifndef AURENHAMMER_HPP
#define AURENHAMMER_HPP

#include <boost/timer.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <OT/line_search.hpp>

namespace OT
{
  template <class Function, class LineSearch, class Vector>
  void
  gradient_descent(const Function &objective,
		   const LineSearch &line_search,
		   Vector &position,
		   double epsilon,
		   size_t max_iter)
  {
   typename Function::Workspace wsp;

    Vector gradient;
    objective.eval(wsp, position.begin(), position.end());
    wsp.gradient(gradient);
    double value = wsp.value();

    double nr_gradient = 0.0;
    double step = 1.0;
    size_t niter = 0;

    do
      {
	const Vector direction = -gradient;
	int N = line_search(objective, wsp, position, value,
			    gradient, direction, step);	
	nr_gradient = boost::numeric::ublas::norm_inf(gradient);
	std::cerr << "[" << niter << "]:\t"
		  << "f(x) = " << value << "\t"
		  << "sup(Df(x)) = " << nr_gradient << "\n";
      }
    while (nr_gradient >= epsilon && niter++ < max_iter);
    std::cerr << __FUNCTION__ << ": " << niter << " iterations\n";
    std::cerr << "value: " << value << "\n";
  }

  template <class Function>
  void
  rt_aurenhammer_gradient_descent(const Function &objective,
				  uvector &weights,
				  double eps,
				  size_t max_iter = 1000)
  {
    gradient_descent
      (objective, OT::Line_search_backtracking<> (), weights, eps,
       max_iter);
  }

}

#endif
