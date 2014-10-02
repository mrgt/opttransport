#ifndef OT_BUNDLE_DESCENT_HPP
#define OT_BUNDLE_DESCENT_HPP

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_ieee_utils.h>
#include "gsl_multimin_fsdf.h"
#include <OT/misc.hpp>

namespace OT
{
  namespace details
  {
    template<class X>
    void ot_fsdf (const gsl_vector *x, void *instance,
		  double *f, gsl_vector *g) 
    {
      *f = reinterpret_cast<X*>(instance)->evaluate
	(x->data, x->data + x->size, g->data);
    }

    template<class X>
    double ot_f (const gsl_vector *x, void *instance) 
    {
      return reinterpret_cast<X*>(instance)->evaluate_func
	(x->data, x->data + x->size);
    }

    template<class X>
    void ot_sdf (const gsl_vector *x, void *instance,
		 gsl_vector *g)
    {
      reinterpret_cast<X*>(instance)->evaluate_grad
	(x->data, x->data + x->size, g->data);
    }

    double
    norm_inf(const gsl_vector *v)
    {
      double res = 0.0;
      for (size_t i = 0; i < v->size; ++i)
	{
	  res = std::max(res, fabs(v->data[i]));
	}
      return res;
    }
  }


  template <class Function>
  void
  rt_bundle_descent(const Function &objective,
		    uvector &weights,
		    double eps)
  {
    const size_t max_iter = 100, N = weights.size();    
    size_t iter, bundle_size_max;
    int status;
	
    const gsl_multimin_fsdfminimizer_type *T =
      gsl_multimin_fsdfminimizer_bundle_method;
    gsl_multimin_fsdfminimizer *s =  gsl_multimin_fsdfminimizer_alloc(T,N);
	
    gsl_multimin_function_fsdf function;
    function.f = details::ot_f<Function>;
    function.sdf = details::ot_sdf<Function>;
    function.fsdf = details::ot_fsdf<Function>;
    function.n = N; 
    function.params = const_cast<void*> ((const void*) &objective);

    gsl_vector *start_point = gsl_vector_alloc(N);
    for (size_t i = 0; i < N; ++i)
      {
	OT_DEBUG_SHOW(weights[i]);
	gsl_vector_set(start_point, i, weights[i]);
      }
    
    bundle_size_max =  N+3;	
    status = gsl_multimin_fsdfminimizer_set
      (s, &function, start_point, bundle_size_max);
	
    printf("== k ===== f(x) ===== ||sgr_f(x)|| ======= eps ======  \n");
    printf("%4d  %14.7f  %13.8e  %13.8e\n",
	   iter, s->f,
	   details::norm_inf(gsl_multimin_fsdfminimizer_subgradient(s)),
	   s->eps);
    
    iter = 0;
    do
      {
	iter++;
	status = gsl_multimin_fsdfminimizer_iterate(s);
	status = gsl_multimin_test_convergence(s,eps);
	
	double nr = details::norm_inf(s->subgradient);
	printf("%4d  %14.7f  %13.8e  %13.8e\n",
	       iter, s->f, nr, s->eps);
	
	//if(nr <= eps)	 
	//    break;	 
      }
    while( /*status == GSL_CONTINUE &&*/ iter <= max_iter);

    gsl_vector *v = gsl_multimin_fsdfminimizer_x(s);

    // size_t j;
    // printf("\nMinimum is found at\n");
    // for(j=0; j < s->x->size; j++)
    //   {
    // 	printf("%9.6f ",
    // 	       gsl_vector_get(v, j));
    //   }
    // printf("\n\n");  

    //std::cerr << "with value " << s->f << "\n";

    std::copy(v->data, v->data + N, weights.begin());
    gsl_multimin_fsdfminimizer_free(s);
  }
}

#endif
