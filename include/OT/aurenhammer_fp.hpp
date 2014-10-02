#ifndef OT_AURENHAMMER_FP_HPP
#define OT_AURENHAMMER_FP_HPP

#define DEBUG_AURENHAMMER 0

#if (DEBUG_AURENHAMMER >= 4)
#define OT_AUR_DEBUG_SHOW(x) OT_DEBUG_SHOW(x)
#else
#define OT_AUR_DEBUG_SHOW(x) (void) 0
#endif

namespace OT
{
  namespace details
  {
    template <class RT>
    typename RT::Vertex_handle
    rt_set_weight (RT &dt, typename RT::Vertex_handle v, double weight)
    {
      typename RT::Point wp = v->point();
      typename RT::Point new_wp (wp.point(), weight);
      
      dt.remove (v);
      typename RT::Vertex_handle it =  dt.insert (new_wp);

      return it;
    }

    template <class RT, class OutputIterator>
    void rt_copy_vertices(const RT &rt, 
			  OutputIterator out)
    {
      typedef typename RT::Finite_vertices_iterator Finite_vertices_iterator;
      typedef typename RT::Hidden_vertices_iterator Hidden_vertices_iterator;
      typedef typename RT::Vertex_handle Vertex_handle;

      for (Finite_vertices_iterator it = rt.finite_vertices_begin();
	   it != rt.finite_vertices_end(); ++it)
	*out++ = it;
      for (Hidden_vertices_iterator it = rt.hidden_vertices_begin();
	   it != rt.hidden_vertices_end(); ++it)
	*out++ = it;
    }

    template <class RT, class Density, class MassMap>
    void
    rt_display_result(RT &dt,
		      Density density,
		      MassMap mass)
    {
#if (DEBUG_AURENHAMMER >= 2)
      typedef typename RT::Finite_vertices_iterator Finite_vertices_iterator;

      size_t i = 0;
      for (Finite_vertices_iterator it = dt.finite_vertices_begin();
	   it != dt.finite_vertices_end(); ++it, ++i)
	{
	  std::cerr << "vertex " << (i+1) << ": mass = "
		    << OT::mass(dt, it, density)
		    << " (target_mass = " << mass[it->point().point()] << "),"
		    << " weight = " << it->point().weight() << "\n";
	}
#endif
    }


    template <class RT, class Density>
    typename RT::Vertex_handle
    rt_aurenhammer_vertex_step(RT &dt, typename RT::Vertex_handle v,
			       Density density,
			       double target_mass,
			       double target_error,
			       double *new_weight,
			       double *new_mass = 0.0)
    {
      double min_weight, max_weight, weight = v->point().weight();

      //v = rt_set_weight(dt, v, weight);
      double orig_mass = OT::mass(dt, v, density);

      OT_AUR_DEBUG_SHOW(orig_mass);
      OT_AUR_DEBUG_SHOW(target_mass);

      if (orig_mass <= target_mass)
	{
	  min_weight = weight;

	  double add_weight = 1.0;
	  double current_mass;
	  do 
	    {
	      max_weight = weight + add_weight;
	      add_weight *= 2.0;
	      v = rt_set_weight(dt, v, max_weight);
	      OT_AUR_DEBUG_SHOW(max_weight);
	      current_mass = OT::mass(dt, v, density);
	      OT_AUR_DEBUG_SHOW(current_mass);
	    } while (current_mass <= target_mass);
	}
      else
	{
	  max_weight = weight;

	  double sub_weight = 1.0;
	  double current_mass;
	  do 
	    {
	      min_weight = weight - sub_weight;
	      sub_weight *= 2.0;
	      v = rt_set_weight(dt, v, min_weight);
	      OT_AUR_DEBUG_SHOW(min_weight);
	      current_mass = OT::mass(dt, v, density);
	      OT_AUR_DEBUG_SHOW(current_mass);
	    } while (OT::mass(dt, v, density) >= target_mass);
	}

      v = rt_set_weight(dt, v, weight);

      OT_AUR_DEBUG_SHOW(min_weight);
      OT_AUR_DEBUG_SHOW(max_weight);

      size_t N = 0;
      do 
	{
	  double mid_weight = (min_weight + max_weight)/2.0;

	  OT_AUR_DEBUG_SHOW(N++);
	  OT_AUR_DEBUG_SHOW(target_mass);
	  OT_AUR_DEBUG_SHOW(min_weight);
	  OT_AUR_DEBUG_SHOW(mid_weight);
	  OT_AUR_DEBUG_SHOW(max_weight);      

	  v = rt_set_weight(dt, v, mid_weight);
	  double mass = OT::mass(dt, v, density);
	  OT_AUR_DEBUG_SHOW(mass);

	  if (fabs (mass - target_mass) <= target_error)
	    {
	      *new_weight = mid_weight;
	      if (new_mass) *new_mass = mass;
	  
	      // set back original weight.
	      return rt_set_weight(dt, v, weight);
	    }

	  if (mass >= target_mass)
	    {
	      min_weight = min_weight;
	      max_weight = mid_weight;
	    }
	  else
	    {
	      min_weight = mid_weight;
	      max_weight = max_weight;
	    }
	} while (N < 20);

      return v;
    }

  }

  template <class RT, class Density, class MassMap>
  void
  rt_aurenhammer_step(RT &dt,
		      Density density,
		      MassMap mass,
		      double target_error)
  {
    boost::timer t;

    typedef typename RT::Vertex_handle Vertex_handle;
    typedef typename RT::Point Point;

    std::vector<Vertex_handle> vertices;
    std::vector<Point> updated_positions;

    OT::details::rt_copy_vertices (dt, std::back_inserter(vertices));
  
    typename std::vector<Vertex_handle>::iterator jt = vertices.begin();
    for (size_t i = 0; jt != vertices.end(); ++jt, ++i)
      {
	Vertex_handle v = *jt;

	double target_mass = mass[v->point().point()];
	double new_weight, new_mass;

	v = details::rt_aurenhammer_vertex_step(dt, v, density,
						target_mass, 
						target_error,
						&new_weight,
						&new_mass);
#if (DEBUG_AURENHAMMER >= 3)
	std::cerr << "vertex " << (i+1) << ":\n"
		  << "\told weight = " << v->point().weight() << "\n"
		  << "\tnew weight = " << new_weight << "\n"
		  << "\tnew mass = " << new_mass << "\n";
#endif
	updated_positions.push_back
	  (Point(v->point().point(), new_weight));
      }

    dt.clear();
    dt.insert(updated_positions.begin(), updated_positions.end());  
    OT::details::rt_display_result(dt, density, mass);

    std::cerr << __FUNCTION__ << ": " << t.elapsed() << "\n";
  }


  template <class T>
  class Constant_map
  {
    T _value;

  public:
    Constant_map (const T &value): _value(value) {}

    template<class V>
    T operator [] (V v) const
    {return _value;}

    template<class Pair>
    void insert (Pair p) {}

    template<class V>
    void erase (V v) {}
  };

  template <class RT, class Density>
  void
  rt_aurenhammer_step(RT &dt,
		      Density density,
		      double target_mass,
		      double target_error)
  {
    return rt_aurenhammer_step (dt, density, 
				Constant_map<double>(target_mass),
				target_error);
  }
}

#endif
