/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all 2D
generalized weights, which can be computed analytically.

\cgalHasModel
- `CGAL::Generalized_weights::Uniform_weight_2`
- `CGAL::Generalized_weights::Three_point_family_weight_2`
- `CGAL::Generalized_weights::Inverse_distance_weight_2`
- `CGAL::Generalized_weights::Shepard_weight_2`
- `CGAL::Generalized_weights::Wachspress_weight_2`
- `CGAL::Generalized_weights::Authalic_weight_2`
- `CGAL::Generalized_weights::Mean_value_weight_2`
- `CGAL::Generalized_weights::Tangent_weight_2`
- `CGAL::Generalized_weights::Discrete_harmonic_weight_2`
- `CGAL::Generalized_weights::Cotangent_weight_2`
*/
class AnalyticWeight_2 {

public:

  /*!
    computes a chosen weight at the `query` point given its three 2D neighbors
    `vm` - previous neighbor (m stands for minus), `vj` - jth neighbor, and
    `vp` - next neighbor (p stands for plus).

    Here, `query` is usually a query point inside a polygon, `vj` is the jth vertex
    of a polygon, `vm` is the previous polygon vertex, and `vp` is the next polygon
    vertex in the counterclockwise order.

    This configuration is arising when the weight is computed on a 2D plane.

    \return weight
  */
  const FT operator()(
    const Point_2& query,
    const Point_2& vm,
    const Point_2& vj,
    const Point_2& vp) const
  { }

  /*!
    computes a chosen weight at the `query` point given its three 3D neighbors
    `vm` - previous neighbor (m stands for minus), `vj` - jth neighbor, and
    `vp` - next neighbor (p stands for plus).

    Here, `query` is usually the center vertex of the one-ring neighborhood in a
    polygon mesh, `vj` is the jth neighbor of the center vertex, `vm` is the previous
    neighbor, and `vp` is the next neighbor in the counterclockwise order.

    This configuration is arising when the weight is computed on a 2D surface
    of a polygon mesh.

    \return weight
  */
  const FT operator()(
    const Point_3& query,
    const Point_3& vm,
    const Point_3& vj,
    const Point_3& vp) const
  { }
};
