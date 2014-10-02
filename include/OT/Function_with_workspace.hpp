#ifndef FUNCTION_WITH_WORKSPACE
#define FUNCTION_WITH_WORKSPACE

namespace OT
{
  template <class Function, class StoppingCriterion, class Progress>
  class Function_with_workspace
  { 
  public:
    typedef typename Function::Workspace Workspace;
    typedef typename Workspace::Gradient_type Gradient_type;

  private:
    Workspace _workspace;
    const Function &_function;
    const StoppingCriterion &_criterion;
    Progress &_progress;
    
  public:
    Function_with_workspace(const Function &function,
			    const StoppingCriterion &criterion,
			    Progress &progress) :
      _function(function),
      _criterion(criterion),
      _progress(progress)
    {}

    template <class InputIterator>
    bool eval(InputIterator begin, InputIterator end)
    {
      bool result = function().eval(workspace(), begin, end);
      return result;
    }

    double value ()
    {
      return workspace().value();
    }

    const Gradient_type &gradient () const
    {
      return workspace().gradient();
    }

    bool should_stop ()
    {
      return _criterion (_function, _workspace);
    }

    template <class InputIterator>
    void progress (InputIterator begin, InputIterator end)
    {
      _progress (_function, begin, end, _workspace);
    }
    
    const Function &function() const
    {
      return _function;
    }
    
    const Workspace &workspace () const
    {
      return _workspace;
    }

    Workspace &workspace()
    {
      return _workspace;
    }
  };
  
  template <class InputIterator>
  double norm_p(InputIterator begin, InputIterator end, double p)
  {
    OT_ASSERT (p >= 1.0);

    double psum (0);
    for (; begin != end; ++begin)      
      psum += pow(fabs(*begin), p);

    return pow(psum, 1/p);    
  }

  template <class InputIterator>
  double norm_2(InputIterator begin, InputIterator end)
  {
    return norm_p(begin, end, 2.0);
  }

  template <class InputIterator>
  double norm_1(InputIterator begin, InputIterator end)
  {
    return norm_p(begin, end, 1.0);
  }

  template <class InputIterator>
  double norm_inf(InputIterator begin, InputIterator end)
  {
    double m (0);    
    for (; begin != end; ++begin)
	  m = std::max(fabs(*begin), m);
    return m;
  }

  class Display_value_gradient
  {
    size_t _k;
  public:
    Display_value_gradient () : 
      _k(0)
    {}

    template <class Function, class InputIterator, class Workspace>
    void
    operator () (const Function &function,
		 InputIterator begin, InputIterator end,
		 const Workspace &workspace) 
    {
      double ninf = OT::norm_inf(workspace.gradient().begin(),
				 workspace.gradient().end());
      double n2 = OT::norm_2(workspace.gradient().begin(),
			     workspace.gradient().end());
      std::cerr << "[" << _k << "]: "
		<< "\tf = " << workspace.value() << ", "
		<<  "\tsup(Df) = " << ninf 
		<< "\teucl(Df) = " << n2 << "\n";
      _k++;
    }

    void reset()
    {
      _k = 0;
    }
  };

  class Stop_with_threshold
  {
  protected:
    double _threshold;

  public:
    Stop_with_threshold (double eps) : _threshold(eps) 
    {}

    double threshold() const
    {
      return _threshold;
    }
  };

  class Stop_norm_p_gradient : public Stop_with_threshold
  {
    double _p;
  public:
    Stop_norm_p_gradient(double threshold, double p) : 
      Stop_with_threshold(threshold),
      _p(p)
    {}

    template <class Function, class Workspace>
    bool
    operator () (const Function &function,
		 const Workspace &workspace) const
    {
      return OT::norm_p(workspace.gradient().begin(), 
			workspace.gradient().end(), _p) <= threshold();
    }
  };

  class Stop_norm_inf_gradient : public Stop_with_threshold
  {
  public:
    Stop_norm_inf_gradient(double threshold) : 
      Stop_with_threshold(threshold)
    {}

    template <class Function, class Workspace>
    bool
    operator () (const Function &function,
		 const Workspace &workspace) const
    {
      return OT::norm_inf(workspace.gradient().begin(), 
			  workspace.gradient().end()) <= threshold();
    }
  };

  // template <class Function, class StoppingCriterion>
  // typename Function_with_workspace<Function, StoppingCriterion, class Progress &progress>
  // make_function_with_workspace(const Function &function,
  // 			       const StoppingCondition &condition,
  // 			       Progress &progress) 
  // {
  //   return typename Function_with_workspace<Function, StoppingCondition> (function, condition);
  // }

  // template <class Function, class StoppingCriterion>
  // typename Function_with_workspace<Function, StoppingCriterion, Display_value_and_gradient>
  // make_function_with_workspace(const Function &function,
  // 			       const StoppingCondition &condition) 
  // {
  //   Display_value_and_gradient dvv;
  //   return make_function_with_workspace(function, condition);
  // }
}

#endif
