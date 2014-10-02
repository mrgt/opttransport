#ifndef OT_LINESEARCH_HPP
#define OT_LINESEARCH_HPP

#include <OT/misc.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>

namespace OT
{
  enum
    {
      LINESEARCH_BACKTRACKING_ARMIJO,
      LINESEARCH_BACKTRACKING_WOLFE,
      LINESEARCH_BACKTRACKING_STRONG_WOLFE
    };

  enum
    {
      LINESEARCH_ERROR = -1
    };

  std::ostream &
  operator << (std::ostream &os, const uvector &v)
  {
    os << "[";
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(os, ", "));
    return os << "]";
  }

  template<int Type = LINESEARCH_BACKTRACKING_STRONG_WOLFE>
  class Line_search_backtracking
  {
    const double _ftol;
    const double _wolfe;
    const double _min_step, _max_step;
    size_t _max_num_steps;

  public:
    Line_search_backtracking(double ftol = 1e-4,
			     double wolfe = 0.9,
			     double min_step = 1e-20,
			     double max_step = 1e20,
			     size_t max_num_steps = 20) :
      _ftol(ftol),
      _wolfe(wolfe),
      _min_step(min_step),
      _max_step(max_step),
      _max_num_steps(max_num_steps)
    {}

    template <class Function, class Vector>
    int operator() (const Function &f, 
		    typename Function::Workspace &wsp,
		    Vector &position,
		    double &value,
		    Vector &gradient,
		    const Vector &direction,
		    double &step) const
    {
      typedef OT::uvector Vector;
      Vector start_position = position;
      double start_value = value;

      const double dec = 0.5, inc = 2.1;
      double dginit = ublas::inner_prod(gradient, direction);
      double dgtest = dginit * _ftol;

      if (step <= 0)
	return LINESEARCH_ERROR;

      if (dginit > 0)
	{
	  std::cerr << "Line_search_backtracking: direction is not a descent direction\n";
	  return LINESEARCH_ERROR;
	}

      size_t N = 0;
      for (;;)
	{
	  double width = 0.0;	 
	  position = start_position + step * direction;
	  f.eval(wsp, position.begin(), position.end());
	  value = wsp.value();
	  wsp.gradient(gradient);
	  ++N;

	  //OT_DEBUG_SHOW(start_position);
	  //OT_DEBUG_SHOW(value);
	  //OT_DEBUG_SHOW(start_value);
	  //OT_DEBUG_SHOW(step);
	  //OT_DEBUG_SHOW(direction);
	  //OT_DEBUG_SHOW(_ftol);
	
	  if (value > start_value + step * dgtest)
	    width = dec;
	  else
	    {
	      // The sufficient decrease condition (Armijo condition).
	      if (Type == LINESEARCH_BACKTRACKING_ARMIJO)
		return N;
	    
	      // Check the Wolfe condition. 
	      double dg = ublas::inner_prod(gradient, direction);
	      if (dg < _wolfe * dginit)
		  width = inc;
	      else
		{
		  // Exit with the regular Wolfe condition. 
		  if(Type == LINESEARCH_BACKTRACKING_WOLFE)
		    return N;
	      
		  // Check the strong Wolfe condition. 
		  if(dg > - _wolfe * dginit)
		    width = dec;
		  else
		    return N; 
		}
	    }
	  
	  if ( (step <= _min_step) ||
	       (step >= _max_step) ||
	       (N >= _max_num_steps))
	    return LINESEARCH_ERROR;
	  step *= width;
	}

      return 0;
    }
  };

  class Line_search_constant_step
  {
    const double _step;

  public:
    Line_search_constant_step(double step) :
      _step(step)
    {}

    template <class Function, class Vector>
    int operator() (const Function &f, 
		    typename Function::Workspace &wsp,
		    Vector &position,
		    double &value,
		    Vector &gradient,
		    const Vector &direction,
		    double &step) const
    {
      position = position + _step * direction;
      step = _step;

      // update value and gradient
      f.eval(wsp, position.begin(), position.end());
      value = wsp.value();
      wsp.gradient(gradient);
    }
  };

}  

#endif
