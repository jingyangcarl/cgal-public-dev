/*!
\ingroup PkgShapeRegularizationRefConcepts
\cgalConcept

A concept that describes the set of methods used by the class 
`CGAL::Shape_regularization::Contour_regularization_2` 
to estimate principal directions of the contour.

\cgalHasModel 
- `CGAL::Shape_regularization::Contours::Longest_direction_2`,
- `CGAL::Shape_regularization::Contours::Multiple_directions_2`,
- `CGAL::Shape_regularization::Contours::User_defined_directions_2`
*/
class ContourDirections {

public:

  /*!  
    \brief orients a given `segment` with the index `query_index` with respect 
    to the best estimated contour direction.

    `CGAL::Shape_regularization::Contour_regularization_2` calls this function once  
    for each contour edge.
  */
  void orient(
    const std::size_t query_index, 
    Segment_2& segment) { 

  }
};
