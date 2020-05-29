/*!
\ingroup PkgShapeRegularizationRefConcepts
\cgalConcept

A concept that describes the set of methods used by the class
`CGAL::Shape_regularization::QP_regularization` to access various data
required for setting up the the global regularization problem.

\cgalHasModel
- `CGAL::Shape_regularization::Segments::Angle_regularization_2`,
- `CGAL::Shape_regularization::Segments::Offset_regularization_2`.
*/
class RegularizationType {

public:

  /*!
    returns the max bound on a regularization characteristic (angle-orientation/
    distance-offset/etc.) with respect to which a geometric object with the index
    `query_index` is being regularized.

    `CGAL::Shape_regularization::QP_regularization` calls this method
    once for each object from the input range.
  */
  typename GeomTraits::FT bound(
    const std::size_t query_index) const {

  }

  /*!
    returns an objective function value between two geometric objects, which are
    direct neighbors, that is they form a neighbor pair `i <-> j`. Neighbors are
    provided by the concept `NeighborQuery`.

    `CGAL::Shape_regularization::QP_regularization` calls this method
    once for each neighbor pair being regularized.
  */
  typename GeomTraits::FT target(
    const std::size_t i,
    const std::size_t j) {

  }

  /*!
    updates regularization characteristics (angle-orientation/distance-offset/etc.)
    of the geometric objects being regularized using values from `solution`, one
    value per one regularized object. These values depend on what is being regularized,
    they could be angle or offset differences for example. The solution vector is
    computed by the `QPSolver`.

    Number of values in `solution` equals to the number n of geometric objects being
    regularized + the number m of neighbor pairs between these objects. The first
    n values are the values that should be used.

    `CGAL::Shape_regularization::QP_regularization` calls this method
    once after the global regularization QP problem has been solved.
  */
  void update(
    const std::vector<GeomTraits::FT>& solution) {

  }
};
