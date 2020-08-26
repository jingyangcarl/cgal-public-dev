namespace CGAL {
namespace Weights {

/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that defines the signature of all free functions for computing
2D and 3D three-point weights that is weights, which require only three input points.

Given three points `p`, `q`, and `r` this functions returns the computed one-value weight
with respect to these points. The parameter `traits` here is the traits class with
all necessary geometric objects, predicates, and constructions.

The type `GeomTraits::Point` must be either
`GeomTraits::Point_2` or `GeomTraits::Point_3`.

\cgalHasModel
- `CGAL::Weights::tangent`
- `CGAL::Weights::cotangent`
- `CGAL::Weights::uniform_area`
- `CGAL::Weights::triangular_area`
- `CGAL::Weights::barycentric_area`
- `CGAL::Weights::voronoi_area`
- `CGAL::Weights::mixed_voronoi_area`
*/
template<typename GeomTraits>
const typename GeomTraits::FT three_point_weight(
  const typename GeomTraits::Point& p,
  const typename GeomTraits::Point& q,
  const typename GeomTraits::Point& r,
  const GeomTraits& traits) { }

} // namespace Weights
} // namespace CGAL
