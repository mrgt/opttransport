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
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>


#include "gsl_cqp.h"
#include "qp_solver.h"


int
solve_qp_pdip(gsl_matrix *Q, gsl_vector *q, gsl_vector *solution)
{

	int status;
	
	size_t k=1, i;
	
	const gsl_cqpminimizer_type *T;
	gsl_cqpminimizer *s;
	
	
	gsl_cqp_data *cqp; 
	cqp = malloc(sizeof(gsl_cqp_data));
	
	/* problem data */
	/* Q */
	cqp->Q = Q;
	
	/* q */
	cqp->q = q;
	
	/* A */
	cqp->A = gsl_matrix_alloc(1, q->size);
	gsl_matrix_set_all(cqp->A, 1.0);
	
	
	/* b */
	cqp->b = gsl_vector_alloc(1);
	gsl_vector_set_all(cqp->b, 1.0);
	
	/* C */
	cqp->C = gsl_matrix_alloc(q->size,q->size);
	gsl_matrix_set_identity(cqp->C);
	
	/* d */
	cqp->d = gsl_vector_calloc(q->size);
	
	
	T = gsl_cqpminimizer_mg_pdip;
	s = gsl_cqpminimizer_alloc(T, q->size, 1, q->size);
	
	status = gsl_cqpminimizer_set(s, cqp);	
	
	
	status = gsl_cqpminimizer_test_convergence(s, 1.0e-10, 1.0e-10);
	
	
	while(status == GSL_CONTINUE && k < 1000)
	{
		status = gsl_cqpminimizer_iterate(s);	
		
		
		status = gsl_cqpminimizer_test_convergence(s, 1.0e-12, 1.0e-12);
		
		k++;
		
	}
	
	
	status = gsl_blas_dcopy(gsl_cqpminimizer_x(s),solution);
	
	
	gsl_cqpminimizer_free(s);
	
	free(cqp);
	
	return GSL_SUCCESS;
	
}


