/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all 2D
generalized weight functions, which can be computed analytically.

\cgalHasModel
- `CGAL::Generalized_weights::Wachspress_weights_2`
- `CGAL::Generalized_weights::Mean_value_weights_2`
- `CGAL::Generalized_weights::Discrete_harmonic_weights_2`
*/
class AnalyticWeights_2 {

public:

  /*!
    fills `weights` with 2D generalized weights computed at a
    given `query` point.
  */
  template<
  typename Point_2,
  typename OutputIterator>
  OutputIterator operator()(
    const Point_2& query,
    OutputIterator weights) {

  }
};
