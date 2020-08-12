* When computing weights one by one, we recompute certain areas/distances multiple times and so using them e.g. in barycentric coordinates is not efficient. What should we do about that?
* I should try to combine mvc and dhc in the orbifold parameterization.
* Why do we have in the orbifold_Tutte_parameterizer "+tangent_weight" while in the MVC_post_processor it is "-tangent_weight"? I have fixed that.
* Do I need to cut-zero the dh weights in the orbifold_Tutte_parameterizer? (line 837) I do cut it for the moment. Otherwise, the results are not equal to the original version.
* Skeletonization uses the weird secure version for the cotangent weights.
* In skeletonization, the final example results are not determenistic.
* Should I move polygon weights from the internal namespace in BC to this package?
* I do not want to make Surface_mesh_shortest_path depend on the Barycentric_coordinates_2 package due to the license issues. In addition, they also use a 3D version for 3D triangles.
* Some packages require traits, which do not have a Point_2. What should we do in that case? See e.g. heat_method_3_concept.
* I should add more construct_objects from the kernel instead of constructors.
* Change to
  using GeomTraits = typename CGAL::Kernel_traits<
      typename boost::property_traits<VertexPointMap>::value_type>::type;
* Remove all classes, change them to free functions with dimensions. Add overloads with no traits.
* Move secure cotangent to the polygon_mesh_tools, rename polygon_mesh_tools to tools.
* Rename Weights.h to pmp_depr_weights.h and clean it up.
* Add a concept test as in the heat_method.
* Add more overloads to the uniform weight.
* Remove my current examples from the user manual and use them only in the reference manual. While in the user manual add more complete examples:
  one that chooses a weight from the list, one that shows how to set up a Laplacian, one that shows how to compute barycentric coordinates, one that shows how to define your own traits class. Make the user manual short.
* Clean up the reference manual.