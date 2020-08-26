namespace CGAL {
namespace Weights {

/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that defines the signature of all free functions for computing
2D and 3D weights with a simple analytic expression based on a query point
and three other points.

Given a query point `q` and three ordered points `p0`, `p1`, and `p2`, this
function returns the computed one-value weight at `q` with respect to these points.
The parameter `traits` here is the traits class with all necessary geometric
objects, predicates, and constructions.

The type `GeomTraits::Point` must be either
`GeomTraits::Point_2` or `GeomTraits::Point_3`.

\cgalHasModel
- `CGAL::Weights::uniform_weight`
- `CGAL::Weights::shepard_weight`
- `CGAL::Weights::inverse_distance_weight`
- `CGAL::Weights::three_point_family_weight`
- `CGAL::Weights::wachspress_weight`
- `CGAL::Weights::authalic_weight`
- `CGAL::Weights::mean_value_weight`
- `CGAL::Weights::tangent_weight`
- `CGAL::Weights::discrete_harmonic_weight`
- `CGAL::Weights::cotangent_weight`
*/
template<typename GeomTraits>
const typename GeomTraits::FT analytic_weight(
  const typename GeomTraits::Point& p0,
  const typename GeomTraits::Point& p1,
  const typename GeomTraits::Point& p2,
  const typename GeomTraits::Point& q,
  const GeomTraits& traits) { }

} // namespace Weights
} // namespace CGAL
