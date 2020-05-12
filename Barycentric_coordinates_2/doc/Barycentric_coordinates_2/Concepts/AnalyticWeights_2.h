/*!
\ingroup PkgBarycentricCoordinates2RefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all
barycentric weight functions, which can be computed or evaluated analytically.

\cgalHasModel
- `CGAL::Barycentric_coordinates::Wachspress_weights_2`
- `CGAL::Barycentric_coordinates::Mean_value_weights_2`
- `CGAL::Barycentric_coordinates::Discrete_harmonic_weights_2`
- `CGAL::Barycentric_coordinates::Harmonic_coordinates_2`
*/
class AnalyticWeights_2 {

public:

  /*!
    fills `weights` with barycentric weights computed at a given `query` point.
  */
  template<
  typename Point_2,
  typename OutputIterator>
  OutputIterator operator()(
    const Point_2& query,
    OutputIterator weights) {

  }
};
