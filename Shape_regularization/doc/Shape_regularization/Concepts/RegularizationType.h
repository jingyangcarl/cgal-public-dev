/*!
\ingroup PkgShapeRegularizationRefConcepts
\cgalConcept

A concept that describes the set of methods used by the class 
`CGAL::Shape_regularization::QP_regularization` to access various data 
required for setting up the the global regularization problem.

\cgalHasModel 
- `CGAL::Shape_regularization::Segments::Angle_regularization_2`, 
- `CGAL::Shape_regularization::Segments::Offset_regularization_2`
*/
class RegularizationType {

public:

  /*!  
  returns the max bound on an item value that is regularized.

  `CGAL::Shape_regularization::QP_regularization` calls this function for each item 
  with the index `query_index` that participates in the regularization process.
  */
  typename GeomTraits::FT bound(
    const std::size_t query_index) const { 

  }

  /*!  
    returns an objective function value between two items, which are direct neighbors.

    `CGAL::Shape_regularization::QP_regularization` calls this function for each neighbor pair 
    `query_index_i <-> query_index_j` that participates in the regularization process.
  */
  typename GeomTraits::FT target(
    const std::size_t query_index_i, 
    const std::size_t query_index_j) {

  }

  /*!
    applies the `result` from the QP solver to the initial items.

    `CGAL::Shape_regularization::QP_regularization` calls this function once after
    the global regularization QP problem has been solved.

    The `result` vector contains values, one value per one regularized item. These 
    values depend on what is being regularized, they could be angles or offsets for example.
  */
  void update(
    const std::vector<GeomTraits::FT>& result) {
    
  }
};
