/*!
\ingroup PkgBarycentricCoordinates2RefConcepts
\cgalConcept

A concept that describes the set of requirements of the template parameter
`GeomTraits` used to parameterize all classes and functions with 2D barycentric
coordinates from the namespace `CGAL::Barycentric_coordinates`.

\cgalHasModel All models of `Kernel`.
*/
class BarycentricTraits_2 {

public:

/// \name Types
/// @{

/*!
	A model of `FieldNumberType`.
*/
typedef unspecified_type FT;

/*!
	`CGAL::Comparison_result` or `Uncertain<CGAL::Comparison_result>`.
*/
typedef unspecified_type Comparison_result;

/*!
	`CGAL::Orientation` or `Uncertain<CGAL::Orientation>`.
*/
typedef unspecified_type Orientation;

/// @}

/// \name 2D Geometric Objects
/// @{

/*!
	A model of `Kernel::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
	A model of `Kernel::Vector_2`.
*/
typedef unspecified_type Vector_2;

/// @}

/// \name 2D Generalized Constructions
/// @{

/*!
	A model of this concept must provide an `FT operator(const Point_2& p, const Point_2& q, const Point_2& r)`
	that returns the signed area of the triangle defined by the points `p`, `q`, and `r`.
*/
typedef unspecified_type Compute_area_2;

/*!
  A model of this concept must provide an `FT operator(const Point_2& p, const Point_2& q)`
  that returns the squared Euclidean distance between the points `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_2;

/*!
  A model of this concept must provide an `FT operator(const Vector_2& p)`
  that returns the squared length of the vector `p`.
*/
typedef unspecified_type Compute_squared_length_2;

/*!
  A model of this concept must provide an `FT operator(const Vector_2& p, const Vector_2& q)`
  that returns the scalar product of the vectors `p` and `q`.
*/
typedef unspecified_type Compute_scalar_product_2;

/*!
  A model of this concept must provide an `FT operator(const Vector_2& p, const Vector_2& q)`
  that returns the determinant of the vectors `p` and `q`.
*/
typedef unspecified_type Compute_determinant_2;

/// @}

/// \name 2D Generalized Predicates
/// @{

/*!
	A model of this concept must provide a `bool operator(const Point_2& p, const Point_2& q)`
	that returns `true` if `p = q` and `false` otherwise.
*/
typedef unspecified_type Equal_2;

/*!
	A model of this concept must provide a `bool operator(const Point_2& p, const Point_2& q, const Point_2& r)`
	that returns `true` if the points `p`, `q`, and `r` are collinear and `false` otherwise.
*/
typedef unspecified_type Collinear_2;

/*!
	A model of this concept must provide a `bool operator(const Point_2& p, const Point_2& q, const Point_2& r)`
	that returns `true` if the point `q` lies between the points `p` and `r` and all three points are collinear.
*/
typedef unspecified_type Collinear_are_ordered_along_line_2;

/*!
	A model of this concept must provide a `bool operator(const Point_2& p, const Point_2& q)`
	that returns `true` iff the x-coordinate of `p` is smaller than the x-coordinate of `q` or
	if they are the same and the y-coordinate of `p` is smaller than the y-coordinate of `q`.
*/
typedef unspecified_type Less_xy_2;

/*!
	A model of this concept must provide a `Comparison_result operator(const Point_2& p, const Point_2& q)`
	that compares the Cartesian x-coordinates of the points `p` and `q`.
*/
typedef unspecified_type Compare_x_2;

/*!
	A model of this concept must provide a `Comparison_result operator(const Point_2& p, const Point_2& q)`
	that compares the Cartesian y-coordinates of the points `p` and `q`.
*/
typedef unspecified_type Compare_y_2;

/*!
	A model of this concept must provide an `Orientation operator(const Point_2& p, const Point_2& q, const Point_2& r)`
	that returns `CGAL::LEFT_TURN` if `r` lies to the left of the oriented line `l` defined by `p` and `q`,
	returns `CGAL::RIGHT_TURN` if `r` lies to the right of `l`,
	and returns `CGAL::COLLINEAR` if `r` lies on `l`.
*/
typedef unspecified_type Orientation_2;

/// @}

};
