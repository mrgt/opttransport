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
#include <math.h>
#include <gsl/gsl_blas.h>

#include "gsl_multimin_fsdf.h"

int
gsl_multimin_test_convergence (gsl_multimin_fsdfminimizer * s, double epsabs)
{
	
	if (epsabs < 0.0)
	{
		GSL_ERROR ("absolute tolerance is negative", GSL_EBADTOL);
	}
	
	if(gsl_blas_dnrm2(s->subgradient) <= epsabs && s->eps <= epsabs)
	{
		return GSL_SUCCESS;
	}


  return GSL_CONTINUE;
}
