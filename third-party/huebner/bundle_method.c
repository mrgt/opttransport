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
#include <string.h>
#include <stdio.h>


#include <math.h>
#include <gsl/gsl_blas.h>


#include "gsl_multimin_fsdf.h"
#include "qp_solver.h"


typedef struct bundle_element
{
	gsl_vector *sgr; 	/* subgradient */ 
	double lin_error;	/* linearization error */
	
	size_t state;		/* state = 0, if the subgradient belongs to the current iteration point
				state = 1, if the subgradient belong to an iteration point, which is calculated  in the iterations process
				state = 2, if the subgradient is the aggregate subgradient*/
	
	struct bundle_element *next;
	struct bundle_element *previous;
	
}bundle_element;

typedef struct
{
	bundle_element *head;
	bundle_element *tail;
	
	double v;
	
	/* t is the step length*/
	double t;
	double t_min;
	double t_max;
	
	size_t fixed_step_length;
	
	/* for t-update */
	int step_counter;
	
	/* descent tolerance */
	double m_ss; 	/*for serius step*/
	double m_t;	/*for update of the step length */
	double m_rel;	/*for the relaxation step */
	
	/* size of the bundle */
	size_t bundle_size_max;
	size_t bundle_size;
	
	/* tolerance of the computation */
	double lambda_min;
	double lm_accuracy;
	double rg;
	
	/* for the relaxation step */
	size_t relaxation;	/* relaxation =0 : no relaxation step will be executed */
	size_t fixed_relaxation;	/* relaxation step with constant rel. parameter */
	double rel_parameter;
	int rel_counter;
	int rel_counter_max;
	
	/* counter for the function and subgradient evaluation */
	size_t f_eval;
	size_t sgr_eval;
	
	size_t serious_step;
	
} bundle_method_state_t;

/* add the new subgradient und lin. error (a new bundle element) to the bundle */
static int
update_bundle(bundle_method_state_t *state, const gsl_vector *new_sgr, const double lin_error_sgr, const gsl_vector *lambda,
	      const gsl_vector *aggr_sgr, const double lin_error_aggr_sgr, const size_t serious_step);

/* remove an element 'item' from bundle and free the space */
static int
remove_element(bundle_element *item, bundle_element **head, bundle_element **tail);

/* insert a new element before the tail */
static int
insert_element(bundle_element *item, bundle_element **head, bundle_element **tail);

/* build the data for the dual problem */
static int
build_cqp_data(const bundle_method_state_t *state, gsl_matrix *Q, gsl_vector *q);

/* get out the bundle elements */
static void 
bundle_out_liste(bundle_method_state_t * state);


static int bundle_method_alloc (void *vstate, size_t n)
{
	bundle_method_state_t * state = (bundle_method_state_t *) vstate;
  
	bundle_element *item;

	/* allocate the space for the initial bundle */ 
	item = malloc(sizeof(bundle_element));
	if (item == 0)
	{
		GSL_ERROR ("failed to allocate the space for the initial element of the bundle", GSL_ENOMEM);
	}
	
	item->sgr = gsl_vector_alloc(n);
	if (item->sgr == 0)
	{
		free(item);
		GSL_ERROR ("failed to allocate the space for a subgradient in the initial element of the bundle", GSL_ENOMEM);
	}
	
	state->head = item;
	state->tail = item;
	
	state->head->next = NULL;
	state->head->previous = NULL;
	
	return GSL_SUCCESS;
}

/********************************************************************************/

static int
bundle_method_set(void *vstate, gsl_multimin_function_fsdf * fsdf,
                  const gsl_vector * x, double *f, gsl_vector * subgradient,
                  double *eps, size_t bundle_size_max)
{
	bundle_method_state_t *state = (bundle_method_state_t *) vstate;
	
	int status, debug=0;
	
	state->rg = 1.0e-20;
	
	/* initialize the first bundle element with a subgradient at the start point,
	the lin. error at this point (it is of course zero),
	the state of the bundle elemet is seted to 0 (current iteration point) */
	GSL_MULTIMIN_FN_EVAL_F_SDF (fsdf, x, f, subgradient);
	
	status = gsl_blas_dcopy(subgradient, state->head->sgr);
	state->head->lin_error = 0.0;
	state->head->state = 0;
	
	*eps = 0.0;
	
	state->f_eval = 1;
	state->sgr_eval = 1;
	
	state->serious_step = 0;
	
	if(debug)
	{
		size_t i;
		printf("\n the start point:\n");
		for(i=0; i<x->size; i++)
			printf("%12.8f ", gsl_vector_get(x,i));
		printf("\n");
		printf("the function value at the start point:   f=%f\n",*f);
		printf("\n the subgradient at the start point:\n");
		for(i=0; i<subgradient->size; i++)
			printf("%12.8f ", gsl_vector_get(subgradient,i));
		printf("\n");
	}
	
	
	/* the initial length step, and lower and upper bound for t */
	state->fixed_step_length = 0;
	
	if(state->fixed_step_length)
	{
		state->t_max = 1;
		state->t = 1;
		state->t_min = 1;
	}
	else
	{
		state->t_max = 1.0e+10;
	
		if(gsl_blas_dnrm2(subgradient) > state->rg)
			state->t = 1.0/gsl_blas_dnrm2(subgradient);
		else
			state->t = state->t_max;
	
		state->t_min = (state->t)*(1.0e-10);
	}
	
	
	state->step_counter = 0;
	
	state->m_ss  = 1.0e-1;
	state->m_t   = 5.0e-1;
	state->m_rel = 3.0e-1;
	
	
	if(bundle_size_max < 3)
	{
		GSL_ERROR ("the maximal number of bundle elements must be greater or equal to 3", GSL_ENOMEM);
	}
	
	state->bundle_size_max = bundle_size_max; 
	state->bundle_size = 1;
	
	state->lm_accuracy = 1.0e-12;
	state->lambda_min = (1.0/(state->t))*state->lm_accuracy;
	
	state->relaxation = 1;
	state->fixed_relaxation = 0;
	state->rel_parameter = 1.2;
	state->rel_counter = 0;
	state->rel_counter_max = 3;
	
	
	return GSL_SUCCESS;
}

/********************************************************************************/

static int
update_bundle(bundle_method_state_t *state, const gsl_vector *new_sgr, const double lin_error_sgr, const gsl_vector *lambda,
	      const gsl_vector *aggr_sgr, const double lin_error_aggr_sgr, const size_t serious_step)
{
	bundle_element *current;
	bundle_element *item;
	bundle_element *item_aggr_sgr;
	bundle_element *item_largest_lin_error;
	
	size_t i;
	
	int status;

	
	/* at first: drop all inactive bundle elements, i.e. elements with lambda(j)=0, because they do not contribute to the trial point y */
	/* second: drop the bundle elemen (if necessary ) with the largest linearization error */
	current = state->head;
	item_largest_lin_error = NULL;
	item_aggr_sgr = NULL;
	
	for(i=0; i<lambda->size; i++)
	{
		if(!serious_step && current->state == 2)
		{
			item_aggr_sgr = current;
			current = current->next;
			
		}
		else if((fabs(gsl_vector_get(lambda,i))<state->lambda_min && current->state != 0) || (serious_step && current->state == 2))
		{
			item = current;
			current = current->next;
			
			status = remove_element(item, &(state->head), &(state->tail));
			
			(state->bundle_size)--;
			
		}
		else 
		{
			if(item_largest_lin_error == NULL || fabs(current->lin_error) > fabs(item_largest_lin_error->lin_error))
				item_largest_lin_error = current;
			
			if (current->state ==0 && serious_step)
				current->state = 1;
			
			current = current->next;
			
			
		}
	}
	
	if(state->bundle_size >= state->bundle_size_max)
	{
		status = remove_element(item_largest_lin_error, &(state->head), &(state->tail));
		(state->bundle_size)--;
		
		if(state->bundle_size >= state->bundle_size_max-1 && !serious_step && item_aggr_sgr == NULL)
		{
			if(state->head->state != 0)
				item_largest_lin_error = state->head;
			else
				item_largest_lin_error = state->head->next;
			
			for(item = item_largest_lin_error->next; item != NULL; item=item->next)
			{
				if(fabs(item->lin_error) > fabs(item_largest_lin_error->lin_error) && item->state != 0)
					item_largest_lin_error = item;
			}
			
			status = remove_element(item_largest_lin_error, &(state->head), &(state->tail));
			(state->bundle_size)--;
		}
	}
	
	
	/* add the new element to the bundle */
	item = malloc(sizeof(bundle_element));
	if (item == 0)
	{
		GSL_ERROR ("failed to allocate space for a new bundle element", GSL_ENOMEM);
	}
	
	item->sgr = gsl_vector_alloc(new_sgr->size);
	if (item->sgr == 0)
	{
		free(item);
		GSL_ERROR ("failed to allocate space for a subgradient in the new bundle element", GSL_ENOMEM);
	}
	
	status = gsl_blas_dcopy(new_sgr,item->sgr);
	item->lin_error = lin_error_sgr;
	
	if(serious_step)
		item->state = 0;
	else
		item->state = 1;
	
	status = insert_element(item, &(state->head), &(state->tail));
	
	(state->bundle_size)++;
	
	
	if(!serious_step)
	{
		if(item_aggr_sgr == NULL)
		{
			item = malloc(sizeof(bundle_element));
			if (item == 0)
			{
				GSL_ERROR ("failed to allocate space for a new bundle element", GSL_ENOMEM);
			}
	
			item->sgr = gsl_vector_alloc(new_sgr->size);
			if (item->sgr == 0)
			{
				free(item);
				GSL_ERROR ("failed to allocate space for a subgradient in the new bundle element", GSL_ENOMEM);
			}
	
			status = gsl_blas_dcopy(aggr_sgr,item->sgr);
			item->lin_error = lin_error_aggr_sgr;
			item->state = 2;
	
			status = insert_element(item, &(state->head), &(state->tail));
	
			(state->bundle_size)++;
		}
		else
		{
			status = gsl_blas_dcopy(aggr_sgr,item_aggr_sgr->sgr);
			item_aggr_sgr->lin_error = lin_error_aggr_sgr;
		}
	}
	
	return GSL_SUCCESS;

} 

/********************************************************************************/

static int
bundle_method_iterate (void *vstate, gsl_multimin_function_fsdf * fsdf, gsl_vector * x, double * f, 
                       gsl_vector * subgradient, gsl_vector * dx, double * eps)
{
	bundle_method_state_t *state = (bundle_method_state_t *) vstate;
	
	bundle_element *item;
	
	size_t i, debug=0;
	
	int status;
	double tmp_d, t_old, t_int_l; /* local variables */
	
	gsl_vector *y;		/* a trial point (the next iteration point by the serios step) */
	gsl_vector *sgr_y;	/* subgradient at y */
	double f_y;		/* the function value at y */
	
	gsl_vector *p;			/* the aggregate subgradient */
	double p_norm, lin_error_p;	/* norm of p, the aggregate linear. error */ 
	gsl_vector *tmp_v;
	
	/* data for the convex quadratic problem (for the dual problem) */
	gsl_vector *q;		/* elements of the array are the linearization errors */
	gsl_matrix *Q;		/* Q=G^T*G (G is matrix which collumns are subgradients) */
	gsl_vector *lambda;	/*  the convex combination coefficients of the subgradients (solution of the dual problem) */
	
	
	lambda = gsl_vector_alloc(state->bundle_size);
	if(lambda == 0)
	{
		GSL_ERROR_VAL ("failed to allocate workspace", GSL_ENOMEM, 0);
	}
	
	q = gsl_vector_alloc(lambda->size);
	if(q == 0)
	{
		gsl_vector_free(lambda);
		GSL_ERROR_VAL ("failed to allocate workspace", GSL_ENOMEM, 0);
	}
	
	y = gsl_vector_calloc(x->size);
	if(y == 0)
	{
		gsl_vector_free(q);
		gsl_vector_free(lambda);
		GSL_ERROR_VAL ("failed to allocate workspace", GSL_ENOMEM, 0);
	}
	
	sgr_y = gsl_vector_calloc(x->size);
	if(sgr_y == 0)
	{
		gsl_vector_free(y);
		gsl_vector_free(q);
		gsl_vector_free(lambda);
		GSL_ERROR_VAL ("failed to allocate workspace", GSL_ENOMEM, 0);
	}
	
	Q = gsl_matrix_alloc(state->bundle_size, state->bundle_size);
	if(Q == 0)
	{
		gsl_vector_free(sgr_y);
		gsl_vector_free(y);
		gsl_vector_free(q);
		gsl_vector_free(lambda);
		GSL_ERROR_VAL ("failed to allocate workspace", GSL_ENOMEM, 0);
	}
	
	p = gsl_vector_calloc(x->size);
	if(p == 0)
	{
		gsl_matrix_free(Q);
		gsl_vector_free(sgr_y);
		gsl_vector_free(y);
		gsl_vector_free(q);
		gsl_vector_free(lambda);
		GSL_ERROR_VAL ("failed to allocate workspace", GSL_ENOMEM, 0);
	}
	
	tmp_v = gsl_vector_calloc(x->size);
	if(tmp_v == 0)
	{
		gsl_vector_free(p);
		gsl_matrix_free(Q);
		gsl_vector_free(sgr_y);
		gsl_vector_free(y);
		gsl_vector_free(q);
		gsl_vector_free(lambda);
		GSL_ERROR_VAL ("failed to allocate workspace", GSL_ENOMEM, 0);
	}
	
	/* solve the dual problem */
	status = build_cqp_data(state, Q, q);
	
	status = solve_qp_pdip(Q, q, lambda);	
	
	gsl_matrix_free(Q);
	gsl_vector_free(q);
	
	
	/* compute the aggregate subgradient (it is called p in the documantation)*/
	/* and the appropriated linearization error */
	
	lin_error_p = 0.0;
	item = state->head;
	for(i=0; i<lambda->size; i++)
	{
		status = gsl_blas_daxpy(gsl_vector_get(lambda,i), item->sgr, p);
		lin_error_p += gsl_vector_get(lambda,i)*(item->lin_error);
		
		item = item->next;
	}
	
	
	if(debug)
	{
		printf("the dual problem solution:\n");
		for(i=0;i<lambda->size;i++)
			printf("%7.6e ",gsl_vector_get(lambda,i));
		printf("\n\n");
		
		printf("the aggregate subgradient: \n");
		for(i=0;i<p->size;i++)
			printf("%.6e ",gsl_vector_get(p,i));
		printf("\n");
		
		printf("lin. error for aggr subgradient = %e\n",lin_error_p);
	}
	
	/* the norm of the aggr subgradient */
	p_norm = gsl_blas_dnrm2(p);
		
	/* search direction dx=-t*p (t is the length of step) */
	status = gsl_vector_memcpy(dx,p);
	status = gsl_vector_scale(dx,-1.0*state->t);
	
	
	/* v =-t*norm(p)^2-alpha_p */
	state->v = -gsl_pow_2(p_norm)*(state->t)-lin_error_p;
	
	/* the subgradient is the aggegate sungradient */
	status = gsl_blas_dcopy(p,subgradient);
		
	/* iteration step */	
	/* y=x+dx */
	status = gsl_blas_dcopy(dx,y);
	status = gsl_blas_daxpy(1.0,x,y);
	
	/* function value at y */
	f_y = GSL_MULTIMIN_FN_EVAL_F(fsdf, y);
	
	state->f_eval++;
	
	/* for t-update */
	if(!state->fixed_step_length)
	{
		t_old = state->t;
		if(fabs(state->v-(f_y-*f)) < state->rg || state->v-(f_y-*f) > state->rg)
			t_int_l = state->t_max;
		else
			t_int_l = 0.5*t_old*(state->v)/(state->v-(f_y-*f));
	}
	else
	{
		t_old = state->t;
		t_int_l = state->t;
	}
	
	
	if( f_y-*f <= state->m_ss*state->v ) /* Serious-Step */
	{
		
		if(debug)
			printf("\nSerious-Step\n");
		
		/* the relaxation step */
		if(state->relaxation)
		{
			if(f_y-*f <= state->v*state->m_rel)
			{
				double f_z;
			
				gsl_vector * z = gsl_vector_alloc(y->size);
			
				/* z = y+dx = x+2*dx */
				status = gsl_blas_dcopy(x,z);
				status = gsl_blas_daxpy(2.0,dx,z);
			
				f_z = GSL_MULTIMIN_FN_EVAL_F(fsdf, z);
				state->f_eval++;
				
				if(0.5*f_z-f_y+0.5*(*f) > state->rg)
					state->rel_parameter = GSL_MIN_DBL(-0.5*(-0.5*f_z+2.0*f_y-1.5*(*f))/(0.5*f_z-f_y+0.5*(*f)),1.999);
				else if (fabs(0.5*f_z-f_y+0.5*(*f)) > state->rg)
					state->rel_parameter = 1.999;
				else
					/* something is wrong */
					state->rel_parameter = 1.0;
								
				
				/* save the old iteration point */
				status = gsl_blas_dcopy(y,z);
				
				/* y = (1-rel_parameter)*x+rel_parameter*y */
				gsl_blas_dscal(state->rel_parameter,y);
				status = gsl_blas_daxpy(1.0-state->rel_parameter,x,y);
				
				/* f(y) und sgr_f(y) */
				tmp_d = GSL_MULTIMIN_FN_EVAL_F(fsdf, y);
				state->f_eval++;
				if(tmp_d > f_y)
				{
					/* keep y as the current point */
					status = gsl_blas_dcopy(z,y);
					
					state->rel_counter++;	
					
				}				
				else
				{
					f_y = tmp_d;
					/* dx = y-x */
					status = gsl_blas_dcopy(y,dx);
					status = gsl_blas_daxpy(-1.0,x,dx);
					
					/* if iteration points bevor and after the rel. step are closly,
					the rel_step counte will be increased */
					/* |1-rel_parameter| <= 0.1*/
					if( fabs(1.0-state->rel_parameter) < 0.1)
						state->rel_counter++;	
				}
				
				
				GSL_MULTIMIN_FN_EVAL_SDF(fsdf, y, sgr_y);
				state->sgr_eval++;
				
				if(state->rel_counter > state->rel_counter_max)
					state->relaxation = 0;
				
				/* */
				status = gsl_blas_daxpy(-1.0,y,z);
				status = gsl_blas_ddot(p, z, &tmp_d);
				*eps = f_y-*f-(state->v)+tmp_d;
				
				gsl_vector_free(z);
			}
			else
			{
				*eps = f_y-(state->v)-*f;
				GSL_MULTIMIN_FN_EVAL_SDF(fsdf, y, sgr_y);
				state->sgr_eval++;
			}
		}
		else
		{
			*eps = f_y-(state->v)-*f;
			
			GSL_MULTIMIN_FN_EVAL_SDF(fsdf, y, sgr_y);
			state->sgr_eval++;
		}
		
		/* calculate linearization errors at new iteration point  */
		item = state->head;
		for(i=0; i<state->bundle_size; i++)
		{
			status = gsl_blas_ddot(item->sgr, dx, &tmp_d);
			item->lin_error += f_y-*f-tmp_d;
			
			item = item->next;
		}
		
		/*  linearization error at new iteration point  */
		status = gsl_blas_ddot(p, dx, &tmp_d);
		lin_error_p += f_y-*f-tmp_d;
		
		/* update the bundle  */
		status = update_bundle(state, sgr_y, 0.0, lambda, p, lin_error_p, 1);
		
		/* adapt the step length */
		if(!state->fixed_step_length)
		{
			if(f_y-*f <= state->v*state->m_t && state->step_counter > 0)
				state->t = t_int_l;
			else if(state->step_counter>3)
				state->t=2.0*t_old;
		
			state->t = GSL_MIN_DBL(GSL_MIN_DBL(state->t,10.0*t_old),state->t_max);
			/*state->eps_v = GSL_MAX_DBL(state->eps_v,-2.0*state->v);*/
		
			state->step_counter = GSL_MAX_INT(state->step_counter+1,1);
				
			if(fabs(state->t-t_old) > state->rg) 
				state->step_counter=1;
		}
		
		
		/* x=y, f=f(y) */
		status = gsl_blas_dcopy(y,x);
		*f = f_y;
	 
		
	}
	else /* Null-Step */
	{	
		
		if(debug)
		  printf("\nNull-Step\n");
		
		GSL_MULTIMIN_FN_EVAL_SDF(fsdf, y, sgr_y);
		state->sgr_eval++;
		
		/* eps for the eps_subdifferential */
		*eps = lin_error_p;
		
		/*calculate the liniarization error at y */
		status = gsl_blas_ddot(sgr_y,dx,&tmp_d);
		tmp_d += *f-f_y;
		
		/* Bundle update */
		status = update_bundle(state, sgr_y, tmp_d, lambda, p, lin_error_p, 0);
		
		/* adapt the step length */
		if(!state->fixed_step_length)
		{
			/*state->eps_v = GSL_MIN_DBL(state->eps_v,lin_error_p);*/
		
			if(tmp_d > GSL_MAX_DBL(p_norm,lin_error_p) && state->step_counter < -1)
				state->t = t_int_l;
			else if(state->step_counter < -3)
				state->t = 0.5*t_old;
		
			state->t = GSL_MAX_DBL(GSL_MAX_DBL(0.1*t_old,state->t),state->t_min);
		
			state->step_counter = GSL_MIN_INT(state->step_counter-1,-1);
				
			if(fabs(state->t-t_old) > state->rg) 
				state->step_counter = -1;
		}

		
	}
	
	
	state->lambda_min = p_norm * state->lm_accuracy;

	if(debug)
	{  
	  
	  printf("\nthe new bundle:\n");
	  bundle_out_liste(state);
  
	  printf("\n\n");
	
	  printf("the curent itarationspoint (1 x %d)\n",x->size);
	  for(i=0;i<x->size;i++)
		  printf("%12.6f ",gsl_vector_get(x,i)); 
	  printf("\n\n");	
	
	  printf("functions value at current point: f=%.8f\n",*f);
	
	  printf("\nstep length t=%.5e\n",state->t);
	  
	  printf("\nstep_counter sc=%d\n",state->step_counter);
	
	  printf("\naccuracy: v=%.5e\n",state->v);
	
	  printf("\nlambda_min=%e\n",state->lambda_min);
  
	  printf("\n");
	}
	
	gsl_vector_free(lambda);
	gsl_vector_free(y);
	gsl_vector_free(sgr_y);
	gsl_vector_free(p);
	
	return GSL_SUCCESS;
}

/* for the debug */
static void 
bundle_out_liste(bundle_method_state_t * state)
{

	size_t i;
	bundle_element *item;
	
	printf("bundle: (%d x %d)\n",state->head->sgr->size, state->bundle_size);
	
	for(item = state->head; item != NULL; item=item->next)
	{
		printf("%14d ",item->state);
	}
		
	printf("\n");
	
	i=0;
	for(i=0; i<state->head->sgr->size; i++)
	{
		for(item = state->head; item != NULL; item=item->next)
		{
			printf("%14.8f ",gsl_vector_get(item->sgr,i));
		}
		printf("\n");
	}
	
	printf("\n");
	
	printf("lin. error: (1 x %d)\n",state->bundle_size);
	for(item = state->head; item != NULL; item=item->next)
	{
		printf("%10.8e ",item->lin_error);
	}
		
	printf("\n");
	
}

static int
remove_element(bundle_element *item, bundle_element **head, bundle_element **tail)
{
	if( ((*head) == NULL) || (item == NULL))
		return -1;
	
	if(item->next == NULL) 
		(*tail) = item->previous;
	else
		item->next->previous = item->previous;
	
	if(item->previous == NULL) 
		(*head) = item->next;
	else
		item->previous->next = item->next;
	
	gsl_vector_free(item->sgr);
	free(item);
	
	return GSL_SUCCESS; 
}

static int
insert_element(bundle_element *item, bundle_element **head, bundle_element **tail)
{
	if(item == NULL)
		return -1;
	
	if(tail == NULL)
		(*head) = item;
	else
		(*tail)->next = item;
	
	item->previous = (*tail);
	item->next = NULL;
	
	(*tail) = item;
	
	return GSL_SUCCESS;
}


static int
build_cqp_data(const bundle_method_state_t *state, gsl_matrix *Q, gsl_vector *q)
{
	size_t i,j;
	
	bundle_element *item_i;
	bundle_element *item_j;
	
	double r;
	
	int status;
	
	item_i = state->head;
	
	for(i=0; i<state->bundle_size; i++)
	{
		/* matrix of the object function (the matrix is symmetric)*/
		item_j = item_i;
		for(j=i; j<state->bundle_size; j++)
		{
			status = gsl_blas_ddot(item_i->sgr, item_j->sgr, &r);
			gsl_matrix_set(Q, i, j, r);
			
			if(j!=i)
				gsl_matrix_set(Q, j, i, r);
			
			item_j = item_j->next;
		}
		
		/* vector of the object function */
		gsl_vector_set(q, i, (item_i->lin_error)/(state->t));
		
		item_i = item_i->next;
	}
	
	return GSL_SUCCESS;
}


static void
bundle_method_free (void *vstate)
{
	bundle_method_state_t *state = (bundle_method_state_t *) vstate;
	
	bundle_element *item;
	
	item=state->head;
	
	while(item != NULL)
	{
		if(item->next != NULL)
		{
			item = item->next;
			gsl_vector_free(item->previous->sgr);
			free(item->previous);
		}
		else
		{
			gsl_vector_free(item->sgr);
			free(item);
			item = NULL;
		}
	}
}

static int
bundle_method_restart (void *vstate, gsl_multimin_function_fsdf *fsdf, const gsl_vector *x, double *f, gsl_vector *subgradient,
                       double *eps)
{
	bundle_method_state_t *state = (bundle_method_state_t *) vstate;
		
	int status;
		
	bundle_method_free(vstate);
	status = bundle_method_alloc(vstate, x->size);
	
	
	status = bundle_method_set(vstate, fsdf, x, f, subgradient, eps, state->bundle_size_max);
	
	return GSL_CONTINUE;
}

static const gsl_multimin_fsdfminimizer_type bundle_method_type = {
  "bundle_method",               /* name */
  sizeof (bundle_method_state_t),
  &bundle_method_alloc,
  &bundle_method_set,
  &bundle_method_iterate,
  &bundle_method_restart,
  &bundle_method_free
};

const gsl_multimin_fsdfminimizer_type
  * gsl_multimin_fsdfminimizer_bundle_method = &bundle_method_type;
