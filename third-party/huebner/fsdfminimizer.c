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
#include <gsl/gsl_errno.h>

#include "gsl_multimin_fsdf.h"

gsl_multimin_fsdfminimizer *
gsl_multimin_fsdfminimizer_alloc (const gsl_multimin_fsdfminimizer_type * T,
                                  size_t n)
{
	int status;
	
	gsl_multimin_fsdfminimizer *s = (gsl_multimin_fsdfminimizer *) malloc (sizeof (gsl_multimin_fsdfminimizer));

	if (s == 0)
	{
		GSL_ERROR_VAL ("failed to allocate space for minimizer struct", GSL_ENOMEM, 0);
	}

	s->type = T;
	
	s->x = gsl_vector_calloc (n);
	if (s->x == 0) 
	{
		free (s);
		GSL_ERROR_VAL ("failed to allocate space for x", GSL_ENOMEM, 0);
	}

	s->subgradient = gsl_vector_calloc (n);
	if (s->subgradient == 0)
	{
		gsl_vector_free (s->x);
		free (s);
		GSL_ERROR_VAL ("failed to allocate space for subgradient", GSL_ENOMEM, 0);
	}
	
	s->dx = gsl_vector_calloc (n);
	if (s->dx == 0)
	{
		gsl_vector_free (s->x);
		gsl_vector_free (s->subgradient);
		free (s);
		GSL_ERROR_VAL ("failed to allocate space for dx", GSL_ENOMEM, 0);
	}
	
	s->state = malloc (T->size);
	if (s->state == 0)
	{
		gsl_vector_free (s->x);
		gsl_vector_free (s->subgradient);
		gsl_vector_free (s->dx);
		free (s);
		GSL_ERROR_VAL ("failed to allocate space for minimizer state",GSL_ENOMEM, 0);
	}
	
	status = (T->alloc) (s->state, n);
	  if (status != GSL_SUCCESS)
	  {
		free (s->state);
		gsl_vector_free (s->x);
		gsl_vector_free (s->subgradient);
		gsl_vector_free (s->dx);
		free (s);
		GSL_ERROR_VAL ("failed to initialize minimizer state", GSL_ENOMEM, 0);
	}
	
	
	return s;
}

int
gsl_multimin_fsdfminimizer_set (gsl_multimin_fsdfminimizer * s,
                                gsl_multimin_function_fsdf * fsdf,
                                const gsl_vector * x,
                                size_t bundle_size)
{
	if (s->x->size != fsdf->n)
	{
		GSL_ERROR ("function incompatible with solver size", GSL_EBADLEN);
	}
	if (x->size != fsdf->n) 
	{
		GSL_ERROR ("vector length not compatible with function", GSL_EBADLEN);
	}
	
	s->fsdf = fsdf;
	gsl_vector_memcpy (s->x,x);
	gsl_vector_set_zero (s->dx);
	
	return (s->type->set) (s->state, s->fsdf, s->x, &(s->f), s->subgradient, &(s->eps), bundle_size);
}

void
gsl_multimin_fsdfminimizer_free (gsl_multimin_fsdfminimizer * s)
{
	(s->type->free) (s->state);
	free (s->state);
	gsl_vector_free (s->dx);
	gsl_vector_free (s->subgradient);
	gsl_vector_free (s->x);
	free (s);
}

int
gsl_multimin_fsdfminimizer_iterate (gsl_multimin_fsdfminimizer * s)
{
	gsl_vector_set_zero (s->dx);
	return (s->type->iterate) (s->state, s->fsdf, s->x, &(s->f), s->subgradient, s->dx, &(s->eps));
}

int
gsl_multimin_fsdfminimizer_restart (gsl_multimin_fsdfminimizer * s)
{
	gsl_vector_set_zero (s->dx);
	return (s->type->restart) (s->state,  s->fsdf, s->x, &(s->f), s->subgradient, &(s->eps));
}

const char *
gsl_multimin_fsdfminimizer_name (const gsl_multimin_fsdfminimizer * s)
{
	return s->type->name;
}


gsl_vector *
gsl_multimin_fsdfminimizer_x (gsl_multimin_fsdfminimizer * s)
{
	return s->x;
}

gsl_vector *
gsl_multimin_fsdfminimizer_dx (gsl_multimin_fsdfminimizer * s)
{
	return s->dx;
}

gsl_vector *
gsl_multimin_fsdfminimizer_subgradient (gsl_multimin_fsdfminimizer * s)
{
	return s->subgradient;
}

double
gsl_multimin_fsdfminimizer_minimum (gsl_multimin_fsdfminimizer * s)
{
	return s->f;
}
