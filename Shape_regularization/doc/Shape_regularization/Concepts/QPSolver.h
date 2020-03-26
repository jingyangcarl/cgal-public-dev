/*!
\ingroup PkgShapeRegularizationRef_Concepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_regularization::QP_regularization` 
to solve the quadratic programming problem.

\cgalHasModel 
- `CGAL::Shape_regularization::CGAL_solver`
- `CGAL::Shape_regularization::OSQP_solver`
*/
class QPSolver {

public:

  /*!  
    computes the quadratic programming problem:
    
    1/2x^{T}Px + q^{T}x subject to the constraints

    l <= Ax <= u, where x in R^{n} is the optimization variable.

    `CGAL::Shape_regularization::QP_regularization` calls this function once per
    each data set.

    The user has to provide two extra variables:
    - `number_of_items` - number of all input items,
    - `number_of_edges` - number of edges in the connectivity graph.

    The matrix and vector types are:
    - `Sparse_matrix = typename Eigen::SparseMatrix<FT, Eigen::ColMajor>`
    - `Dense_vector = typename Eigen::Matrix<FT, Eigen::Dynamic, 1>`

    where `FT` is the field number type.
  */
  void solve(
    const std::size_t number_of_items,
    const std::size_t number_of_edges, 
    const Sparse_matrix& P, 
    const Sparse_matrix& A,
    const Dense_vector& q,
    const Dense_vector& l,
    const Dense_vector& u,
    std::vector<FT>& x) {

  }
};
