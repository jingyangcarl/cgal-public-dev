/*!
\ingroup PkgSolverInterfaceConcepts
\cgalConcept

A concept that describes the set of methods used by the class
`CGAL::Shape_regularization::QP_regularization`
to solve the corresponding quadratic programming (QP) problem.

This concept is an extension of the concept `QuadraticProgram` from the
\ref PkgQPSolver package. It includes all the methods from that concept and
all the methods below.
*/
class LinearProgramTraits {

public:

  /*!
    allocates memory for `n` non-zero values from the matrix `D`.
  */
  void reserve_d(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero values from the vector `c`.
  */
  void reserve_c(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero values from the matrix `A`.
  */
  void reserve_a(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero values from the vector `b`.
  */
  void reserve_b(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero values from the vector `l`.
  */
  void reserve_l(const std::size_t n) { }

  /*!
    allocates memory for `n` non-zero values from the vector `u`.
  */
  void reserve_u(const std::size_t n) { }

  /*!
    \brief solves the quadratic program.

    Number of values in `solution` equals to the number n of geometric objects being
    regularized + the number m of neighbor pairs between these objects.

    \param solution
    a vector with the solution

    \returns a status of the computation `success == true`
  */
  bool solve(
    std::vector<FT>& solution) { }
};
