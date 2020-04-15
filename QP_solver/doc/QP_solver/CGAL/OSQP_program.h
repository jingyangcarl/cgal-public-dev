
namespace CGAL {

  /*!
    \ingroup PkgQPSolverClasses

  \f$
  \newcommand{\qprel}{\gtreqless}
  \newcommand{\qpx}{\mathbf{x}}
  \newcommand{\qpl}{\mathbf{l}}
  \newcommand{\qpu}{\mathbf{u}}
  \newcommand{\qpc}{\mathbf{c}}
  \newcommand{\qpb}{\mathbf{b}}
  \newcommand{\qpy}{\mathbf{y}}
  \newcommand{\qpw}{\mathbf{w}}
  \newcommand{\qplambda}{\mathbf{\lambda}}
  \f$

    An object of class `OSQP_program` describes a convex quadratic program of the form 
    \f{eqnarray*}{
    \mbox{(QP)}& \mbox{minimize} 
    & \qpx^{T}D\qpx+\qpc^{T}\qpx+c_0 \\ 
    &\mbox{subject to} & A\qpx\qprel \qpb, \\ 
    & & \qpl \leq \qpx \leq \qpu 
    \f}
    in \f$ n\f$ real variables \f$ \qpx=(x_0,\ldots,x_{n-1})\f$. 

    Here, 
    <UL>

    <LI>\f$ A\f$ is an \f$ m\times n\f$ matrix (the constraint matrix), 
    <LI>\f$ \qpb\f$ is an \f$ m\f$-dimensional vector (the right-hand side), 
    <LI>\f$ \qprel\f$ is an \f$ m\f$-dimensional vector of relations 
    from \f$ \{\leq, =, \geq\}\f$, 

    <LI>\f$ \qpl\f$ is an \f$ n\f$-dimensional vector of lower 
    bounds for \f$ \qpx\f$, where \f$ l_j\in\mathbb{R}\cup\{-\infty\}\f$ for all \f$ j\f$ 
    <LI>\f$ \qpu\f$ is an \f$ n\f$-dimensional vector of upper bounds for 
    \f$ \qpx\f$, where \f$ u_j\in\mathbb{R}\cup\{\infty\}\f$ for all \f$ j\f$ 

    <LI>\f$ D\f$ is a symmetric positive-semidefinite \f$ n\times n\f$ matrix (the 
    quadratic objective function), 

    <LI>\f$ \qpc\f$ is an \f$ n\f$-dimensional vector (the linear objective 
    function), and 
    <LI>\f$ c_0\f$ is a constant. 

    </UL> 

    This particular implementation is a wrapper of the external \ref thirdpartyOSQP "OSQP solver".

    \cgalModels `QuadraticProgram`
  */
  template<typename FT>
  class OSQP_program { 

  public:
    
    /// \name Types 
    /// @{

    /*!
      The number type of the program entries. 
    */ 
    typedef unspecified_type FT; 

    /// @} 

    /// \name Creation 
    /// @{

    /*!
      constructs a quadratic program with no variables and no constraints, ready 
      for data to be added. Unless relations are explicitly set, they will 
      be of type `default_r`. Unless bounds are explicitly set, they 
      will be as specified by `default_fl` (finite lower bound?), 
      `default_l` (lower bound value if lower bound is finite), 
      `default_fu` (finite upper bound?), and 
      `default_l` (upper bound value if upper bound is finite). If all 
      parameters take their default values, we thus get equality constraints 
      and bounds \f$ x\geq0\f$ by default. Numerical entries that are not 
      explicitly set will default to \f$ 0\f$.\pre if `default_fl == default_fu == true`, then `default_l <= default_u`. 
    */ 
    OSQP_program(
      CGAL::Comparison_result relation = CGAL::SMALLER,
      bool finite_lower = true,
      FT lower = 0,
      bool finite_upper = false,
      FT upper = 0) : 
    m_default_r(relation), 
    m_default_fl(finite_lower),
    m_default_l(lower), 
    m_default_fu(finite_upper), 
    m_default_u(upper), 
    m_is_valid(true) {
    
    CGAL_qpe_assertion(
      !finite_lower || !finite_upper || lower <= upper);
    }

    /// @} 

    /// \name Operations 
    /// @{

    /*!
      sets the entry \f$ A_{ij}\f$ 
      in column \f$ j\f$ and row \f$ i\f$ of the constraint matrix \f$ A\f$ of `qp` to 
      `val`. An existing entry is overwritten. `qp` is enlarged if 
      necessary to accommodate this entry. 
    */ 
    void set_a (int j, int i, const NT& val); 

    /*!
      sets the entry \f$ b_i\f$ 
      of `qp` to `val`. An existing entry is overwritten. 
      `qp` is enlarged if necessary to accommodate this entry. 
    */ 
    void set_b (int i, const NT& val); 

    /*!
      sets the entry \f$ \qprel_i\f$ of `qp` to `rel`. `CGAL::SMALLER` 
      means that the \f$ i\f$-th constraint is of type "\f$ \leq\f$", `CGAL::EQUAL` 
      means "\f$ =\f$", and `CGAL::LARGER` encodes "\f$ \geq\f$". An existing entry 
      is overwritten. `qp` is enlarged if necessary to accommodate this entry. 
    */ 
    void set_r (int i, CGAL::Comparison_result rel); 

    /*!
      if `is_finite`, this sets the entry \f$ l_j\f$ of `qp` to `val`, 
      otherwise it sets \f$ l_j\f$ to \f$ -\infty\f$. An existing entry is overwritten. 
      `qp` is enlarged if necessary to accommodate this entry. 
    */ 
    void set_l (int j, bool is_finite, const NT& val = NT(0)); 

    /*!
      if `is_finite`, this sets the entry \f$ u_j\f$ of `qp` to `val`, 
      otherwise it sets \f$ u_j\f$ to \f$ \infty\f$. An existing entry is overwritten. 
      `qp` is enlarged if necessary to accommodate this entry. 
    */ 
    void set_u (int j, bool is_finite, const NT& val = NT(0)); 

    /*!
      sets the entry \f$ c_j\f$ 
      of `qp` to `val`. An existing entry is overwritten. 
      `qp` is enlarged if necessary to accommodate this entry. 
    */ 
    void set_c (int j, const NT& val); 

    /*!
      sets the entry \f$ c_0\f$ 
      of `qp` to `val`. An existing entry is overwritten. 
    */ 
    void set_c0 (const NT& val); 

    /*!
      sets the entries 
      \f$ 2D_{ij}\f$ and \f$ 2D_{ji}\f$ of `qp` to `val`. Existing entries are 
      overwritten. `qp` is enlarged if necessary to accommodate these entries. 
      \pre `j <= i` 
    */ 
    void set_d (int i, int j, const NT& val); 

    /// @}
  };

} // namespace CGAL
