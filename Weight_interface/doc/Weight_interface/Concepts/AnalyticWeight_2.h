/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all 2D
generalized weight functions, which can be computed analytically.

\cgalHasModel
- `CGAL::Generalized_weights::Uniform_weight_2`
- `CGAL::Generalized_weights::Wachspress_weight_2`
- `CGAL::Generalized_weights::Mean_value_weight_2`
- `CGAL::Generalized_weights::Discrete_harmonic_weight_2`
*/
class AnalyticWeight_2 {

public:

  /*!
    description
  */
  void operator()();
};
