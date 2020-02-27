/*!
\ingroup PkgShapeRegularization_Concepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_regularization::QP_regularization` 
to apply regularization.

\cgalHasModel 
- `CGAL::Shape_regularization::Angle_regularization_2`, 
- `CGAL::Shape_regularization::Ordinate_regularization_2`
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
  typename GeomTraits::FT target_value(
    const std::size_t query_index_i, 
    const std::size_t query_index_j) {

  }

  /*!
    applies the results from the QP solver to the initial items.

    `CGAL::Shape_regularization::QP_regularization` calls this function once, after
    the QP problem has been solved during the regularization process.
  */
  void update(
    const std::vector<GeomTraits::FT>& result) {
    
  }
};
