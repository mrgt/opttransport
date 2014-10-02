#ifndef OT_NEWTON_DESCENT_HPP
#define OT_NEWTON_DESCENT_HPP

#include <boost/timer.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#if O
#include <boost/numeric/bindings/ublas/matrix_sparse.hpp>
#include <boost/numeric/bindings/umfpack/umfpack.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
//#include <boost/numeric/bindings/mumps/mumps_driver.hpp>
#endif

#include "pcg.hpp"

#include <OT/line_search.hpp>

namespace OT
{
  namespace ublas = ::boost::numeric::ublas;
  //typedef ublas::mapped_matrix<double> sumatrix;
  typedef ublas::compressed_matrix<double, ublas::column_major> sumatrix;


  template <class Matrix>
  void
  clamp_matrix(const Matrix &in, Matrix &out, size_t n)
  {
    out = Matrix(n,n);
    
    typedef typename Matrix::const_iterator1 it1_t;
    typedef typename Matrix::const_iterator2 it2_t;
    
    for (it1_t it1 = in.begin1(); it1 != in.end1(); it1++)
      {
	for (it2_t it2 = it1.begin(); it2 != it1.end(); it2++)
	  {
	    if (it2.index1() < n && it2.index2() < n)
	      out.push_back(it2.index2(), it2.index1(), *it2);
	  }
      }
  }

  template <class Matrix>
  void extract_submatrix(const Matrix &in, Matrix &out,
			 const std::vector<size_t> &v)
  {
    std::map<size_t, size_t> s;
    size_t n = v.size() - 1;

    for (size_t i = 0; i < n; ++i)
      {
	s[v[i]] = i;
	std::cerr << v[i] << " ";
      }
    std::cerr << "\n";
    
    out = Matrix(n,n);

    typedef typename Matrix::const_iterator1 it1_t;
    typedef typename Matrix::const_iterator2 it2_t;

    for (it1_t it1 = in.begin1(); it1 != in.end1(); it1++)
      {
	for (it2_t it2 = it1.begin(); it2 != it1.end(); it2++)
	  {
	    std::map<size_t, size_t>::iterator s1 = s.find(it2.index1());
	    std::map<size_t, size_t>::iterator s2 = s.find(it2.index2());

	    if (s1 != s.end() && s2 != s.end())
	      out.push_back(s2->second, s1->second, *it2);
	  }
      }
  }

  template <class Vector>
  void extract_subvector(const Vector &in, Vector &out,
			 const std::vector<size_t> &v)
  {
    size_t n = v.size() - 1;

    out = Vector(n);
    for (size_t i = 0; i < n; ++i)
      out[i] = in[v[i]];
  }

  template <class Function, class LineSearch, class Vector>
  void
  newton_descent(const Function &objective,
		 const LineSearch &line_search,
		 Vector &position,
		 double epsilon,
		 size_t max_iter = 100)
  {
    typename Function::Workspace wsp;

    Vector gradient;
    objective.eval(wsp, position.begin(), position.end());
    wsp.gradient(gradient);
    double value = wsp.value();

    double nr_gradient = 0.0;
    double step = 1.0;
    size_t niter = 0;
    size_t N = position.size();
    
    do
      {
	Vector direction = -gradient;
	sumatrix hm (N,N);
	std::vector<size_t> active;

	if (wsp.hessian(hm, active))  // && active.size() == N)
	  {
	    active.resize(N);
	    for (size_t i = 0; i < N; ++i)
	      active[i] = i;

	    size_t n = active.size() - 1;

	    std::cerr << "newton_descent: " << n << " active sites\n";

	    sumatrix subm (n,n);
	    extract_submatrix(hm, subm, active);

#if 0
	    double mean = 0.0;
	    for (size_t i = 0; i < n; ++i)
	      mean += subm(i,i)/n;

	    OT_DEBUG_SHOW(mean);
	    for (size_t i = 0; i < n; ++i)
	        subm(i,i) += mean/10;
#endif

	    Vector subd;
	    extract_subvector(direction, subd, active);
	    Vector resd = subd;

	    IncompleteCholeskyPreconditioner<sumatrix> P(subm);
	     ::pcg_solve(subm, resd, subd, P);

	    // OT_DEBUG_SHOW(direction);
	    // OT_DEBUG_SHOW(subd);
	     OT_DEBUG_SHOW(resd);
	     OT_DEBUG_SHOW(gradient);

	     std::cerr << "copying back\n";
	     direction = Vector(N);
	    for (size_t i = 0; i < n; ++i)
	      {
		direction[active[i]] = resd[i];
		//std::cerr << "direction [" << active[i] << "] = " << resd[i]
		//	  << "\n";
	      }
	  }

#if 1
	int N = line_search(objective, wsp, position, value,
			    gradient, direction, step);	
#else
	position += 0.9 * direction;
	objective.eval(wsp, position.begin(), position.end());
	value = wsp.value();
	wsp.gradient(gradient);
#endif
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
  newton_descent(const Function &objective,
		 uvector &weights,
		 double eps,
		 size_t max_iter = 1000)
  {
    newton_descent
      (objective, OT::Line_search_backtracking<> (), weights, eps,
       max_iter);
  }

}

#endif
