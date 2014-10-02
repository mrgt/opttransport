/***************************************************************************
 *   Copyright (C) 2004 by Ewgenij Huebner                                  *
 *   huebner@uni-trier.de                                                  *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef GSL_MULTIMIN_FSDF_H
#define GSL_MULTIMIN_FSDF_H

#include <gsl/gsl_multimin.h>

#ifdef __cplusplus
# define __BEGIN_DECLS extern "C" {
# define __END_DECLS }
#else
#define __BEGIN_DECLS 
#define __END_DECLS
#endif


__BEGIN_DECLS

/* Definition of an arbitrary subdifferentiable real-valued function */
/* with gsl_vector input and parameters */

struct gsl_multimin_function_fsdf_struct 
{
  double (* f) (const gsl_vector  * x, void * params);
  void (* sdf) (const gsl_vector * x, void * params, gsl_vector * sdf);
  void (* fsdf) (const gsl_vector * x, void * params, double *f, gsl_vector * sdf);
  size_t n;
  void * params;
};

typedef struct gsl_multimin_function_fsdf_struct gsl_multimin_function_fsdf;

#define GSL_MULTIMIN_FN_EVAL_F(F,x) (*((F)->f))(x,(F)->params)
#define GSL_MULTIMIN_FN_EVAL_SDF(F,x,sg) (*((F)->sdf))(x,(F)->params,(sg))
#define GSL_MULTIMIN_FN_EVAL_F_SDF(F,x,y,sg) (*((F)->fsdf))(x,(F)->params,(y),(sg))


/* minimisation of subdifferentiable functions */

typedef struct 
{
  const char *name;
  size_t size;
  int (*alloc) (void *state, size_t n);
  int (*set) (void *state, gsl_multimin_function_fsdf * fsdf,
              const gsl_vector * x, double * f, 
              gsl_vector * subgradient, double *accuracy, size_t bundle_size);
  int (*iterate) (void *state, gsl_multimin_function_fsdf * fsdf, 
                  gsl_vector * x, double * f, 
		  gsl_vector * subgradient, gsl_vector * dx, double *accuracy);
  int (*restart) (void *state, gsl_multimin_function_fsdf * fsdf,
                  const gsl_vector * x, double * f, 
                  gsl_vector * subgradient, double *accuracy);
  void (*free) (void *state);
}
gsl_multimin_fsdfminimizer_type;


typedef struct 
{
	/* multi dimensional part */
	const gsl_multimin_fsdfminimizer_type *type;
	gsl_multimin_function_fsdf *fsdf;

	double f;
	gsl_vector * x;
	gsl_vector * subgradient;
	gsl_vector * dx;
	double eps;
	
	void *state;
}
gsl_multimin_fsdfminimizer;

gsl_multimin_fsdfminimizer *
gsl_multimin_fsdfminimizer_alloc(const gsl_multimin_fsdfminimizer_type *T,
                                size_t n);

int 
gsl_multimin_fsdfminimizer_set (gsl_multimin_fsdfminimizer * s,
                               gsl_multimin_function_fsdf *fsdf,
                               const gsl_vector * x,
                               size_t bundle_size);

void
gsl_multimin_fsdfminimizer_free(gsl_multimin_fsdfminimizer *s);

const char * 
gsl_multimin_fsdfminimizer_name (const gsl_multimin_fsdfminimizer * s);

int
gsl_multimin_fsdfminimizer_iterate(gsl_multimin_fsdfminimizer *s);

int
gsl_multimin_fsdfminimizer_restart(gsl_multimin_fsdfminimizer *s);

gsl_vector * 
gsl_multimin_fsdfminimizer_x (gsl_multimin_fsdfminimizer * s);

gsl_vector * 
gsl_multimin_fsdfminimizer_dx (gsl_multimin_fsdfminimizer * s);

gsl_vector * 
gsl_multimin_fsdfminimizer_subgradient (gsl_multimin_fsdfminimizer * s);

int
gsl_multimin_test_convergence (gsl_multimin_fsdfminimizer * s, double epsabs);

double 
gsl_multimin_fsdfminimizer_minimum (gsl_multimin_fsdfminimizer * s);

GSL_VAR const gsl_multimin_fsdfminimizer_type *gsl_multimin_fsdfminimizer_bundle_method;


__END_DECLS

#endif
