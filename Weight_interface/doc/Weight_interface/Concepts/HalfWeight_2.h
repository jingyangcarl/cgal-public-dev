/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all 2D
generalized weights, whose half value can be computed analytically and per triangle
while the full value is obtained as a sum of half values computed per two adjacent triangles.

\cgalHasModel
- `CGAL::Generalized_weights::Authalic_weight_2`
- `CGAL::Generalized_weights::Tangent_weight_2`
- `CGAL::Generalized_weights::Cotangent_weight_2`
*/
class HalfWeight_2 {

public:

  /*!
    computes a half weight at the `query` point given its two 2D neighbors
    `vj` - j neighbor and `vp` - next neighbor.

    This configuration is arising when the half weight is computed on a 2D plane.

    \return half weight
  */
  const FT operator()(
    const Point_2& query,
    const Point_2& vj,
    const Point_2& vp) const
  { }

  /*!
    computes a half weight at the `query` point given its two 3D neighbors
    `vj` - j neighbor and `vp` - next neighbor.

    This configuration is arising when the half weight is computed on a 2D surface
    of a polygon mesh.

    \return half weight
  */
  const FT operator()(
    const Point_3& query,
    const Point_3& vj,
    const Point_3& vp) const
  { }
};
