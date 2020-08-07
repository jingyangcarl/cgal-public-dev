/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that describes the set of methods that should be defined for all 2D
generalized weights, which can be computed analytically.

\cgalFigureBegin{neighborhood, neighborhood.svg}
  Notation used for the local neighborhood.
\cgalFigureEnd

\cgalHasModel
- `CGAL::Generalized_weights::Uniform_weight`
- `CGAL::Generalized_weights::Three_point_family_weight`
- `CGAL::Generalized_weights::Inverse_distance_weight`
- `CGAL::Generalized_weights::Shepard_weight`
- `CGAL::Generalized_weights::Wachspress_weight`
- `CGAL::Generalized_weights::Authalic_weight`
- `CGAL::Generalized_weights::Mean_value_weight`
- `CGAL::Generalized_weights::Tangent_weight`
- `CGAL::Generalized_weights::Discrete_harmonic_weight`
- `CGAL::Generalized_weights::Cotangent_weight`
*/
class AnalyticWeight_2 {

public:

  /*!
    computes a weight at the `q` point given its three 2D neighbors
    `t` - previous neighbor, `r` - jth neighbor, and `p` - next neighbor.

    Here, `q` is usually a query point inside a polygon, `r` is the jth vertex
    of a polygon, `t` is the previous polygon vertex, and `p` is the next polygon
    vertex in the counterclockwise order (see the figure above).

    This configuration is arising when the weight is computed on a 2D plane.

    \return weight
  */
  const FT operator()(
    const Point_2& q,
    const Point_2& t,
    const Point_2& r,
    const Point_2& p) const
  { }

  /*!
    computes a weight at the `q` point given its three 3D neighbors
    `t` - previous neighbor, `r` - jth neighbor, and `p` - next neighbor.

    Here, `q` is usually the center vertex of the one-ring neighborhood in a
    polygon mesh, `r` is the jth neighbor of the center vertex, `t` is the previous
    neighbor, and `p` is the next neighbor in the counterclockwise order (see the figure above).

    This configuration is arising when the weight is computed on a 2D surface
    embedded in 3D.

    \return weight
  */
  const FT operator()(
    const Point_3& q,
    const Point_3& t,
    const Point_3& r,
    const Point_3& p) const
  { }
};
