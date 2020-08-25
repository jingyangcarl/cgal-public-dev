/*!
\ingroup PkgWeightInterfaceRefConcepts
\cgalConcept

A concept that describes the set of requirements of the template parameter
`GeomTraits` used to parameterize several classes and functions with 3D generalized
weights from the namespace `CGAL::Weights`.

\cgalHasModel All models of `Kernel`.
*/
class AnalyticWeightTraits_3 {

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
	A model of `Kernel::Point_3`.
*/
typedef unspecified_type Point_3;

/*!
	A model of `Kernel::Vector_3`.
*/
typedef unspecified_type Vector_3;

/*!
	A model of `Kernel::Angle_3`.
*/
typedef unspecified_type Angle_3;

/// @}

/// \name Generalized Constructions
/// @{

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
  A model of this concept must provide a `Vector_3 operator(const Vector_3& v, const Vector_3& w)`
  that returns the cross product of the vectors `v` and `w`.
*/
typedef unspecified_type Construct_cross_product_vector_3;

/*!
  A model of this concept must provide a `Point_3 operator(const Point_3& p, const Point_3& q, const Point_3& r)`
  that returns the center of the circle passing through the points `p`, `q`, and `r`.
*/
typedef unspecified_type Construct_circumcenter_3;

/*!
  A model of this concept must provide a `Vector_3 operator(const Point_3& p, const Point_3& q)`
  that returns the vector through the points `p` and `q`.
*/
typedef unspecified_type Construct_vector_3;

/// @}

};
