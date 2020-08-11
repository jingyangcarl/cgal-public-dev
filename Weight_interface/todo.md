* When computing weights one by one, we recompute certain areas/distances multiple times and so
  using them e.g. in barycentric coordinates is not efficient. What should we do about that?
* I should try to combine mvc and dhc in the orbifold parameterization.
* Why do we have in the orbifold_Tutte_parameterizer "+tangent_weight" while in the MVC_post_processor it is "-tangent_weight"? I have fixed that.
* Do I need to cut-zero the dh weights in the orbifold_Tutte_parameterizer? (line 837) I do cut it for the moment. Otherwise, the results are not equal to
  the original version.
* Skeletonization uses the weird secure version fot the cotangent weights.
* In skeletonization, the final example results are not determenistic.
* Should I move polygon weights from the internal namespace in BC to this package?
* I do not want to make Surface_mesh_shortest_path depend on the Barycentric_coordinates_2 package due to the license issues. In
  addition, they also use a 3D version for 3D triangles.