* Split weights into dimension-based categories: Generalized_weights_1/2/3/d.
* All analytic weights should follow the same concept.
* I should first submit the Barycentric_coordinates_2 package taking into account that
  it will be later modified by using the Weight_interface package.
* Then I should submit the Weight_interface package that will also modify the internals
  of the Barycentric_coordinates_2 package.
* Finally, I will submit a set of features that will modify internals of such packages
  as Interpolation, Polygon_mesh_processing, Surface_mesh_parameterization, Surface_mesh_skeletonization,
  Surface_mesh_deformation, Surface_mesh_shortest_path, Heat_method_3, and possibly Kernel_d.
* When computing weights one by one, we recompute certain areas/distances multiple times and so
  using them e.g. in barycentric coordinates is not efficient. What should we do with that?
* Use the code:
  ```
    /// \cond SKIP_IN_MANUAL
    template<
    typename PolygonMesh,
    typename VertexDescriptor,
    typename VertextAroundTargetCirculator>
    const FT operator()(
      const PolygonMesh& polygon_mesh,
      const VertexDescriptor vdi,
      const VertextAroundTargetCirculator vcj) const {

      const auto point_map = get(vertex_point, polygon_mesh);
      const Point_3& query = get(point_map, vdi);

      auto vcm = vcj; vcm--;
      auto vcp = vcj; vcp++;

      const Point_3& vm = get(point_map, vcm);
      const Point_3& vj = get(point_map, vcj);
      const Point_3& vp = get(point_map, vcp);

      return weight_3(query, vm, vj, vp);
    }
    /// \endcond
  ```
* I should change interface of the half weights and make it such that query is always in the middle +
  rename query, vj, and vp to p, q, r.
* I should try to combine mvc and dhc in the orbifold parameterization.
* I should use free functions in the parameterization.
* In some weights we use area_3, which is always positive. Is it ok?