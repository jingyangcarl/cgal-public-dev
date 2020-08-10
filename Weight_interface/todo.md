* When computing weights one by one, we recompute certain areas/distances multiple times and so
  using them e.g. in barycentric coordinates is not efficient. What should we do about that?
* I should try to combine mvc and dhc in the orbifold parameterization.
* Why do we have in the orbifold_Tutte_parameterizer "+tangent_weight" while in the MVC_post_processor it is "-tangent_weight"?
* Do I need to cut-zero the dh weights in the orbifold_Tutte_parameterizer? (line 837)