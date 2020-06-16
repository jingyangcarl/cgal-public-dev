/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that describes the set of requirements of the template parameter
`GeomTraits` used to parameterize all classes and functions with 2D generalized
weights from the namespace `CGAL::Generalized_weights`.

\cgalHasModel All models of `Kernel`.
*/
class AnalyticTraits_2 {

public:

/// \name Types
/// @{

/*!
	A model of `FieldNumberType`.
*/
typedef unspecified_type FT;

/// @}

/// \name Geometric Objects
/// @{

/*!
	A model of `Kernel::Point_2`.
*/
typedef unspecified_type Point_2;

/*!
	A model of `Kernel::Point_3`.
*/
typedef unspecified_type Point_3;

/*!
	A model of `Kernel::Vector_2`.
*/
typedef unspecified_type Vector_2;

/*!
	A model of `Kernel::Vector_3`.
*/
typedef unspecified_type Vector_3;

/// @}

/// \name Generalized Constructions
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
  A model of this concept must provide an `FT operator(const Vector_2& v)`
  that returns the squared length of the vector `v`.
*/
typedef unspecified_type Compute_squared_length_2;

/*!
  A model of this concept must provide an `FT operator(const Vector_2& v, const Vector_2& w)`
  that returns the scalar product of the vectors `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_2;

/*!
  A model of this concept must provide an `FT operator(const Vector_2& v, const Vector_2& w)`
  that returns the determinant of the vectors `v` and `w`.
*/
typedef unspecified_type Compute_determinant_2;

/*!
  A model of this concept must provide an `FT operator(const Point_3& p, const Point_3& q)`
  that returns the squared Euclidean distance between the points `p` and `q`.
*/
typedef unspecified_type Compute_squared_distance_3;

/*!
  A model of this concept must provide an `FT operator(const Vector_3& v)`
  that returns the squared length of the vector `v`.
*/
typedef unspecified_type Compute_squared_length_3;

/*!
  A model of this concept must provide an `FT operator(const Vector_3& v, const Vector_3& w)`
  that returns the scalar product of the vectors `v` and `w`.
*/
typedef unspecified_type Compute_scalar_product_3;

/*!
  A model of this concept must provide an `Vector_3 operator(const Vector_3& v, const Vector_3& w)`
  that returns the cross product of the vectors `v` and `w`.
*/
typedef unspecified_type Construct_cross_product_vector_3;

/// @}

/// \name Generalized Predicates
/// @{

/*!
	A model of this concept must provide a `bool operator(const Point_2& p, const Point_2& q)`
	that returns `true` if `p = q` and `false` otherwise.
*/
typedef unspecified_type Equal_2;

/*!
	A model of this concept must provide a `bool operator(const Point_3& p, const Point_3& q)`
	that returns `true` if `p = q` and `false` otherwise.
*/
typedef unspecified_type Equal_3;

/// @}

};
