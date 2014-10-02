// This code is described in "Computational Geometry in C" (Second Edition),
// Chapter 7.  It is not written to be comprehensible without the
// explanation in that book.

// Written by Joseph O'Rourke.
// Last modified: December 1997
// Questions to orourke@cs.smith.edu.
// --------------------------------------------------------------------
// This code is Copyright 1997 by Joseph O'Rourke.  It may be freely
// redistributed in its entirety provided that this copyright notice is
// not removed.
// --------------------------------------------------------------------


#ifndef INTERSECTION_OROURKE_HPP
#define INTERSECTION_OROURKE_HPP

#include <CGAL/Polygon_2.h>
#include <CGAL/intersection_2.h>
#include <CGAL/Kernel/global_functions_2.h>

namespace OT
{
  namespace detail
  {
    enum tInFlag {Unknown, Qin, Pin};
    
    // Advances and prints out an inside vertex if appropriate.
    template<class K>
    size_t advance (size_t a, size_t *aa, size_t n,
		    bool inside,
		    typename CGAL::Polygon_2<K> &poly,
		    const typename K::Point_2 &v)
    {
      if (inside)
	poly.push_back(v);
      
      (*aa)++;
      return  (a+1) % n;
    }
  }

  
  // P and Q _have_ to be in ccw order.
  template<class K>
  bool convex_intersection (const typename CGAL::Polygon_2<K> &P,
			    const typename CGAL::Polygon_2<K> &Q,
			    typename CGAL::Polygon_2<K> &r)
  {
    using namespace detail;
    typedef typename K::Segment_2 Segment;
    typedef typename K::Point_2 Point;

    size_t n = P.size(), m = Q.size();
    tInFlag inflag = Unknown; // {Pin, Qin, Unknown}: which inside 
    bool FirstPoint = true;   
    size_t a = 0, b = 0;      // indices on P and Q (resp.)
    size_t aa = 0, ba = 0;    // # advances on a & b indices (after 1st inter.) 

    if (n == 0 || m == 0)
      {
	r = CGAL::Polygon_2<K>();
	return false;
      }
  
    do
      {
	const size_t a1 = (a + n - 1) % n;
	const size_t b1 = (b + m - 1) % m;
      
	CGAL::Orientation
	  cross = CGAL::orientation (P[a] - P[a1], Q[b] - Q[b1]),
	  aHB = CGAL::orientation (Q[b1], Q[b], P[a]),
	  bHA = CGAL::orientation (P[a1], P[a], Q[b]);
	
	// A = [P[a1],P[a]], B = [Q[b1], Q[b]]
	// If A & B intersect, update inflag. 
	CGAL::Object isect_object = CGAL::intersection(Segment(P[a1], P[a]),
						       Segment(Q[b1], Q[b]));
	Point p;
	if (CGAL::assign(p, isect_object))
	  {
	    if (inflag == Unknown && FirstPoint)
	      {
		FirstPoint = false;
		aa = ba = 0;
	      }
	    r.push_back(p);

	    // Update inflag. 
	    if (aHB == CGAL::POSITIVE)
	      inflag = Pin;
	    else if (bHA == CGAL::POSITIVE)
	      inflag = Qin;
	  }

	//-----Advance rules-----
	// Special case: A & B overlap and oppositely oriented. */
	Segment s;
	if (CGAL::assign(s, isect_object))
	  {
	    // the interior of the intersection is empty
	    if ( (P[a1] - P[a]) * (Q[b1] - Q[b]) < 0)
	      return false;

	    // FIXME: not sure about the other cases...
	  }

	if ( (cross == CGAL::COLLINEAR) &&
	     (aHB == CGAL::NEGATIVE) &&
	     (bHA == CGAL::NEGATIVE) )
	  {
	    // Special case: A & B parallel and separated -> P&Q disjoint. 
	    return false;
	  }
	else if ( (cross == CGAL::COLLINEAR) &&
		  (aHB == CGAL::COLLINEAR) &&
		  (bHA == CGAL::COLLINEAR) )
	  {
	    // Special case: A & B collinear. 
	    // Advance but do not output point.
	    if (inflag == Pin)
	      b = advance<K>(b, &ba, m, inflag == Qin, r, Q[b]);
	    else
	      a = advance<K>(a, &aa, n, inflag == Pin, r, P[a]);
	  }      
	else if ( (cross == CGAL::COLLINEAR) ||
		  (cross == CGAL::POSITIVE))
	  {
	    // Generic cases. 
	    if (bHA == CGAL::POSITIVE)
	      a = advance<K> (a, &aa, n, inflag == Pin, r, P[a]);
	    else
	      b = advance<K> (b, &ba, m, inflag == Qin, r, Q[b]);
	  }
	else // if ( cross < 0 ) 
	  {
	    if ( aHB == CGAL::POSITIVE)
	      b = advance<K> (b, &ba, m, inflag == Qin, r, Q[b]);
	    else
	      a = advance<K> (a, &aa, n, inflag == Pin, r, P[a]);
	  }
      }
    while ( ((aa < n) || (ba < m)) && (aa < 2*n) && (ba < 2*m) );

    if (inflag == Unknown) 
      {
	if (P.has_on_bounded_side(Q[0]))
	  r = Q;
	else if (Q.has_on_bounded_side(P[0]))
	  r = P;
      }

    return true;
  }
}

#endif
