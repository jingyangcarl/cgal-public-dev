/*!
\ingroup PkgShapeRegularizationRefConcepts
\cgalConcept

A concept that describes the set of methods used by the class
`CGAL::Shape_regularization::QP_regularization`
to solve the corresponding quadratic programming (QP) problem.

This concept is an extension of the concept `QuadraticProgram` from the
\ref PkgQPSolver package. It includes all the methods from that concept and
all the methods below.

\cgalHasModel
- `CGAL::Shape_regularization::CGAL_quadratic_program`,
- `CGAL::Shape_regularization::OSQP_quadratic_program`
*/
class QPSolver {

public:

  /*!
    allocates memory for `n` non-zero items from the matrix `D`.
  */
  void reserve_d(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero items from the vector `c`.
  */
  void reserve_c(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero items from the matrix `A`.
  */
  void reserve_a(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero items from the vector `b`.
  */
  void reserve_b(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero items from the vector `l`.
  */
  void reserve_l(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero items from the vector `u`.
  */
  void reserve_u(const std::size_t n) { }

  /*!
    \brief solves the quadratic program.

    \param solution
    a vector with the solution

    \returns a status of the computation `success == true`
  */
  bool solve(
    std::vector<FT>& solution) { }
};
