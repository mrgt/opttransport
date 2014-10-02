#include <uQuadProg++.hh>

namespace OT
{
  template <class Container>
  bool project_on_convex_hull(const Container &v,
			      const uvector &origin,
			      uvector &result)
  {
    OT_ASSERT(v.size() != 0);

    // min F(x) = 0.5 * x G x + g0 x
    // s.t.
    // CE^T x + ce0 = 0
    // CI^T x + ci0 >= 0

    // Projection on convex hull of (v_1,...,v_N):
    // minimize F(x) = 0.5 * ||\sum x_i v_i - origin||^2
    //               = 0.5 * ||\sum x_i v_i||^2 + <\sum x_i v_i, origin>
    // under contraints x_1 >= 0, ..., x_N >= 0 and \sum x_i = 1.0
    //
    // <\sum x_i v_i, origin> = \sum x_i <v_i, origin>

    typedef typename  Container::const_iterator Iterator;
    size_t N = v.size(), D = origin.size();

    uvector ce0 = ublas::scalar_vector<double>(1, -1.0);
    umatrix CE = ublas::scalar_matrix<double>(N, 1, 1.0);

    uvector ci0 = ublas::zero_vector<double>(N);
    umatrix CI = ublas::identity_matrix<double>(N);

    uvector g0(N);
    size_t i = 0;
    for (Iterator it = v.begin(); it != v.end(); ++it)
      g0[i] = ublas::inner_prod(*it, origin);

    ublas::matrix<double> G(N,N, 0.0);
    i = 0;
    for (Iterator it = v.begin(); it != v.end(); ++it, ++i)
      {
	size_t j = 0;
	for (Iterator jt = v.begin(); jt != v.end(); ++jt, ++j)
	  G(i,j) = ublas::inner_prod(*it, *jt);
      }

    uvector x;
    double value = uQuadProgPP::solve_quadprog(G, g0, CE, ce0, CI, ci0, x);

    OT_DEBUG_SHOW(value);
    if (value == std::numeric_limits<double>::infinity())
      return false;

    OT_DEBUG_SHOW(x);
    result = ublas::zero_vector<double>(origin.size());
    i = 0;
    for (Iterator it = v.begin(); it != v.end(); ++it, ++i)
      {
	OT_DEBUG_SHOW(*it);
	result += x[i] * (*it);
      }
    OT_DEBUG_SHOW(result);

    return true;
  }

  template <class Function, class LineSearch, class Vector>
  void
  simple_subgradient_descent(const Function &objective,
			     const LineSearch &line_search,
			     Vector &position,
			     size_t num_gradients,
			     double epsilon)
  {
    Vector gradient;
    typename std::list<Vector> gradients;
    double value = objective.evaluate(position, gradient);
    uvector origin = ublas::zero_vector<double>(position.size());

    do
      {
	double step = 1.0;

	gradients.push_back(gradient);
	if (gradients.size() > num_gradients)
	  gradients.pop_front();

	// find direction
	uvector projection;
	if (!OT::project_on_convex_hull(gradients, origin, projection))
	  projection = gradient;
	OT_DEBUG_SHOW(gradient);
	OT_DEBUG_SHOW(projection);
	int N = line_search(objective, position, value,
			    gradient, uvector(-projection), step);	
	OT_DEBUG_SHOW(gradient);
	OT_DEBUG_SHOW(N);
      }
    while (boost::numeric::ublas::norm_inf(gradient) >= epsilon);    
  }


}
