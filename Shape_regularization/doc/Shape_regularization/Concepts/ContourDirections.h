/*!
\ingroup PkgShapeRegularizationRef_Concepts
\cgalConcept

A concept that describes the set of methods used by the `CGAL::Shape_regularization::Contour_regularization_2` 
to estimate principal directions of the contour.

\cgalHasModel 
- `CGAL::Shape_regularization::Longest_principal_direction_2`,
- `CGAL::Shape_regularization::Multiple_principal_directions_2`,
- `CGAL::Shape_regularization::User_defined_principal_directions_2`
*/
class ContourDirections {

public:

  /*!  
    estimates principal `directions` of the open/closed contour, sets angle 
    `bounds` on each found direction, and fills `assigned` that assigns an index of 
    the direction from `directions` to each contour edge.

    `CGAL::Shape_regularization::Contour_regularization_2` calls this function once per contour.
  */
  void estimate(
    const bool is_closed, 
    std::vector< std::pair<FT, FT> >& bounds,
    std::vector<Direction_2>& directions,
    std::vector<std::size_t>& assigned) {
        
  }
};
