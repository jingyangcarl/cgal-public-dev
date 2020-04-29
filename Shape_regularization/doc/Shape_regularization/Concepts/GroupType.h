/*!
\ingroup PkgShapeRegularizationRefConcepts
\cgalConcept

A concept that describes the set of methods used by the classes, which 
group items with respect to a certain regularity type.

\cgalHasModel 
- `CGAL::Shape_regularization::Segments::Parallel_groups_2`,
- `CGAL::Shape_regularization::Segments::Collinear_groups_2`,
- `CGAL::Shape_regularization::Segments::Orthogonal_groups_2`.
*/
class GroupType {

public:

  /*!  
    fills `groups` with `std::vector<std::size_t>`, where each vector contains 
    indices of all items, which belong to the same group.

    \tparam OutputIterator 
    must be a model of `OutputIterator`

    \param groups
    an instance of `OutputIterator`, 
    whose value type is `std::vector<std::size_t>`

    \return an output iterator
  */
  template<typename OutputIterator>
  OutputIterator groups( 
    OutputIterator& groups) const {
  
  }
};
