/*!
\ingroup PkgSolverInterfaceConcepts
\cgalConcept

A concept that describes the set of methods used to define and solve a
quadratic programming (QP) problem of the general form:
<center>
\f{eqnarray*}{
& \mbox{minimize} & \mathbf{x}^{T}D\mathbf{x} + \mathbf{c}^{T}\mathbf{x} + c_0 \\
& \mbox{subject to} & A\mathbf{x} \gtreqless \mathbf{b} \\
& \mbox{and} & \mathbf{l} \leq \mathbf{x} \leq \mathbf{u}
\f}
</center>
in \f$ n \f$ real variables \f$ \mathbf{x} = (x_0, \ldots, x_{n-1}) \f$.

Here,
<UL>
<LI>\f$ D \f$ is a symmetric positive-semidefinite \f$ n \times n\f$ matrix (the quadratic objective function),
<LI>\f$ \mathbf{c} \f$ is an \f$ n \f$-dimensional vector (the linear objective function),
<LI>\f$ c_0 \f$ is a constant,
<LI>\f$ A \f$ is an \f$ m\times n\f$ matrix (the constraint matrix),
<LI>\f$ \mathbf{b} \f$ is an \f$ m \f$-dimensional vector (the right-hand side),
<LI>\f$ \gtreqless \f$ is an \f$ m \f$-dimensional vector of relations
from \f$ \{\leq, =, \geq\} \f$,
<LI>\f$ \mathbf{l} \f$ is an \f$ n \f$-dimensional vector of lower bounds for
\f$ \mathbf{x} \f$, where \f$ l_j \in \mathbb{R} \cup \{-\infty\} \f$ for all \f$ j \f$,
<LI>\f$ \mathbf{u} \f$ is an \f$ n \f$-dimensional vector of upper bounds for
\f$ \mathbf{x} \f$, where \f$ u_j \in \mathbb{R} \cup \{+\infty\} \f$ for all \f$ j \f$.
</UL>
*/
class QuadraticProgramTraits {

public:

  /*!
    Allocates memory for `k` non-zero values from the matrix `D`.
  */
  void reserve_d(const std::size_t k) { }

  /*!
    Allocates memory for `k` non-zero values from the vector `c`.
  */
  void reserve_c(const std::size_t k) { }

  /*!
    Allocates memory for `k` non-zero values from the matrix `A`.
  */
  void reserve_a(const std::size_t k) { }

  /*!
    Allocates memory for `k` non-zero values from the vector `b`.
  */
  void reserve_b(const std::size_t k) { }

  /*!
    Allocates memory for `k` non-zero values from the vector `l`.
  */
  void reserve_l(const std::size_t k) { }

  /*!
    Allocates memory for `k` non-zero values from the vector `u`.
  */
  void reserve_u(const std::size_t k) { }

  /*!
    Sets the entries `2Dij` and `2Dji` of `lp` to `val`.
  */
  void set_d(
    const std::size_t i,
    const std::size_t j,
    const FT val) { }

  /*!
    Sets the entry `cj` of `lp` to `val`.
  */
  void set_c(
    const std::size_t j,
    const FT val) { }

  /*!
    Sets the entry `c0` of `lp` to `val`.
  */
  void set_c0(
    const FT val) { }

  /*!
    Sets the entry `Aij` in the column `j` and row `i` of
    the constraint matrix `A` of `lp` to `val`.
  */
  void set_a(
    const std::size_t j,
    const std::size_t i,
    const FT val) { }

  /*!
    Sets the entry `bi` of `lp` to `val`.
  */
  void set_b(
    const std::size_t i,
    const FT val) { }

  /*!
    If `is_finite`, this sets the entry `lj` of `lp` to `val`,
    otherwise it sets `lj` to \f$-\infty\f$.
  */
  void set_l(
    const std::size_t j,
    const bool is_finite,
    const FT val) { }

  /*!
    If `is_finite`, this sets the entry `uj` of `lp` to `val`,
    otherwise it sets `uj` to \f$+\infty\f$.
  */
  void set_u(
    const std::size_t j,
    const bool is_finite,
    const FT val) { }

  /*!
    \brief Solves the quadratic program.

    Number of values in `solution` equals to the number `n` of values in the
    vector `x`.

    \param solution
    a vector with the solution

    \returns a status of the computation `success == true`
  */
  bool solve(
    std::vector<FT>& solution) { }
};
