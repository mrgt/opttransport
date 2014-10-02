#ifndef OT_BFGS_HPP
#define OT_BFGS_HPP

#include <stdio.h>
#include <lbfgs.h>

#include <OT/Function_with_workspace.hpp>
#include <OT/misc.hpp>

#include <boost/numeric/ublas/matrix_sparse.hpp>


namespace OT
{
  namespace details
  {
    const char *lbgfs_error_message(int nr)
    {
      switch (nr)
	{
	case LBFGS_SUCCESS:
	  return "LBFGS_SUCCESS";
	case LBFGS_STOP:
	  return "LBFGS_STOP";
	case LBFGS_ALREADY_MINIMIZED:
	  return "LBFGS_ALREADY_MINIMIZED";
	case LBFGSERR_UNKNOWNERROR:
	  return "LBFGSERR_UNKNOWNERROR";
	case LBFGSERR_LOGICERROR:
	  return "LBFGSERR_LOGICERROR";
	case LBFGSERR_OUTOFMEMORY:
	  return "LBFGSERR_OUTOFMEMORY";
	case LBFGSERR_CANCELED:
	  return "LBFGSERR_CANCELED";
	case LBFGSERR_INVALID_N:
	  return "LBFGSERR_INVALID_N";
	case LBFGSERR_INVALID_N_SSE:
	  return "LBFGSERR_INVALID_N_SSE";
	case LBFGSERR_INVALID_X_SSE:
	  return "LBFGSERR_INVALID_X_SSE";
	case LBFGSERR_INVALID_EPSILON:
	  return "LBFGSERR_INVALID_EPSILON";
	case LBFGSERR_INVALID_TESTPERIOD:
	  return "LBFGSERR_INVALID_TESTPERIOD";
	case LBFGSERR_INVALID_DELTA:
	  return "LBFGSERR_INVALID_DELTA";
	case LBFGSERR_INVALID_LINESEARCH:
	  return "LBFGSERR_INVALID_LINESEARCH";
	case LBFGSERR_INVALID_MINSTEP:
	  return "LBFGSERR_INVALID_MINSTEP";
	case LBFGSERR_INVALID_MAXSTEP:
	  return "LBFGSERR_INVALID_MAXSTEP";
	case LBFGSERR_INVALID_FTOL:
	  return "LBFGSERR_INVALID_FTOL";
	case LBFGSERR_INVALID_WOLFE:
	  return "LBFGSERR_INVALID_WOLFE";
	case LBFGSERR_INVALID_GTOL:
	  return "LBFGSERR_INVALID_GTOL";
	case LBFGSERR_INVALID_XTOL:
	  return "LBFGSERR_INVALID_XTOL";
	case LBFGSERR_INVALID_MAXLINESEARCH:
	  return "LBFGSERR_INVALID_MAXLINESEARCH";
	case LBFGSERR_INVALID_ORTHANTWISE:
	  return "LBFGSERR_INVALID_ORTHANTWISE";
	case LBFGSERR_INVALID_ORTHANTWISE_START:
	  return "LBFGSERR_INVALID_ORTHANTWISE_START";
	case LBFGSERR_INVALID_ORTHANTWISE_END:
	  return "LBFGSERR_INVALID_ORTHANTWISE_END";
	case LBFGSERR_OUTOFINTERVAL:
	  return "LBFGSERR_OUTOFINTERVAL";
	case LBFGSERR_INCORRECT_TMINMAX:
	  return "LBFGSERR_INCORRECT_TMINMAX";
	case LBFGSERR_ROUNDING_ERROR:
	  return "LBFGSERR_ROUNDING_ERROR";
	case LBFGSERR_MINIMUMSTEP:
	  return "LBFGSERR_MINIMUMSTEP";
	case LBFGSERR_MAXIMUMSTEP:
	  return "LBFGSERR_MAXIMUMSTEP";
	case LBFGSERR_MAXIMUMLINESEARCH:
	  return "LBFGSERR_MAXIMUMLINESEARCH";
	case LBFGSERR_MAXIMUMITERATION:
	  return "LBFGSERR_MAXIMUMITERATION";
	case LBFGSERR_WIDTHTOOSMALL:
	  return "LBFGSERR_WIDTHTOOSMALL";
	case LBFGSERR_INVALIDPARAMETERS:
	  return "LBFGSERR_INVALIDPARAMETERS";
	case LBFGSERR_INCREASEGRADIENT:
	  return "LBFGSERR_INCREASEGRADIENT";
	}

      return "LBFGSERR_UNKNOWN";
    }
  
    template <class Fww>
    lbfgsfloatval_t bfgs_evaluate(void *instance,
				  const lbfgsfloatval_t *x,
				  lbfgsfloatval_t *g,
				  const int n,
				  const lbfgsfloatval_t step)
    {
      Fww *p = reinterpret_cast<Fww*> (instance);
    
      p->eval (x, x+n);

      OT_ASSERT(p->gradient().size() == n);
      std::copy(p->gradient().begin(), p->gradient().end(), g);

      return p->value();
    }

    template <class Fww>
    int bfgs_progress(void *instance,
		      const lbfgsfloatval_t *x,
		      const lbfgsfloatval_t *g,
		      const lbfgsfloatval_t fx,
		      const lbfgsfloatval_t xnorm,
		      const lbfgsfloatval_t gnorm,
		      const lbfgsfloatval_t step,
		      int n,
		      int k,
		      int ls)
    {
      Fww *p = reinterpret_cast<Fww*> (instance);
      p->progress(x, x+n);
      return p->should_stop() ? 1 : 0;
    }
  }

  template <class Function, class StoppingCriterion, class Progress>
  void
  bfgs_descent(const Function &function,
	       const StoppingCriterion &stop,
	       Progress &prog,
	       uvector &weights,
	       int linesearch = LBFGS_LINESEARCH_MORETHUENTE,
	       int max_iterations = 0)
  {
    double value;
    lbfgs_parameter_t param;
    lbfgs_parameter_init(&param);

    param.m = 10;
    param.epsilon = 0;
    param.linesearch = linesearch;
    param.max_iterations = max_iterations;

    typedef typename OT::Function_with_workspace<Function, StoppingCriterion, Progress> Fww;
    Fww functionww (function, stop, prog);

    size_t linesearches [] = {LBFGS_LINESEARCH_MORETHUENTE,
			      LBFGS_LINESEARCH_BACKTRACKING_STRONG_WOLFE,
			      //LBFGS_LINESEARCH_BACKTRACKING_WOLFE,
			      //LBFGS_LINESEARCH_BACKTRACKING_ARMIJO
    };
    size_t num_linesearches = sizeof(linesearches)/sizeof(linesearches[0]);

    int current_line_search = 0;
    for (int i = 0; i < num_linesearches; ++i)
      {
	if (linesearches[i] == linesearch)
	  {
	    current_line_search = i;
	    break;
	  }
      }

    do 
      {
	param.linesearch = linesearches[current_line_search];

	int ret = lbfgs(weights.size(), &weights[0], &value,
			OT::details::bfgs_evaluate<Fww>,
			OT::details::bfgs_progress<Fww>,
			static_cast<void*>(&functionww), &param);

	if (ret == LBFGS_SUCCESS || ret == LBFGS_STOP)
	  return;

	if (ret != LBFGSERR_ROUNDING_ERROR && 
	    ret != LBFGSERR_MAXIMUMLINESEARCH)
	  {
	    std::cerr << "rt_bfgs_descent: "
		      << OT::details::lbgfs_error_message(ret) << "\n";
	    return;
	  }
	
	std::cerr << "rt_bfgs_descent: "
		  << OT::details::lbgfs_error_message(ret) << "\n"
		  << "trying again with weaker linesearch\n";	
	current_line_search++;
      } while (current_line_search < num_linesearches);
  }

  template <class Function, class StoppingCriterion, class Progress>
  void
  bfgs_descent_zm(const Function &function,
		  const StoppingCriterion &stop,
		  Progress &prog,
		  uvector &weights,
		  int linesearch = LBFGS_LINESEARCH_MORETHUENTE,
		  int max_iterations = 0)
  {
    size_t N = function.number_of_vertices();
    OT_ASSERT(weights.size() == N);

    for (size_t i = 0; i < weights.size(); ++i)
      weights[i] -= weights[N-1];
    weights.resize(N-1);

    bfgs_descent (function, stop, prog, weights, linesearch, max_iterations);

    weights.resize(N);
    weights[N-1] = 0;
  }

  template <class Function>
  void
  bfgs_descent(const Function &function,
	       uvector &weights,
	       double eps,
	       int linesearch = LBFGS_LINESEARCH_MORETHUENTE,
	       int max_iterations = 0)
  {
    Display_value_gradient dvg;
    bfgs_descent (function, Stop_norm_inf_gradient(eps), dvg, 
		  weights,
		  linesearch, max_iterations);
  }

  template <class Function>
  void
  bfgs_descent_zm(const Function &function,
		  uvector &weights,
		  double eps,
		  int linesearch = LBFGS_LINESEARCH_MORETHUENTE,
		  int max_iterations = 0)
  {
    Display_value_gradient dvg;
    bfgs_descent_zm (function, Stop_norm_inf_gradient(eps), dvg, 
		  weights,
		  linesearch, max_iterations);
  }
}

#endif

