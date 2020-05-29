/*!
\ingroup PkgShapeRegularizationRefConcepts
\cgalConcept

We assume that each contour has one or several principal directions. By implementing
a model of this cocept, the user sets such directions and provides a way to orient
contour edges towards these directions. All contour regularization functions
in this package are parameterized by this concept.

\cgalHasModel
- `CGAL::Shape_regularization::Contours::Longest_direction_2`,
- `CGAL::Shape_regularization::Contours::Multiple_directions_2`,
- `CGAL::Shape_regularization::Contours::User_defined_directions_2`.
*/
class ContourDirections {

public:

  /*!
    \brief orients a given `segment` with the index `query_index` towards the
    best-fit direction of the contour.
  */
  void orient(
    const std::size_t query_index,
    Segment_2& segment) {

  }
};
