* Comment the code.
* Spellcheck all files.
* Remove all unnecessary functions or put them in experimental.

* Install OSQP on all testing platforms.
* Add a 3D version of the QP regularization algorithm.
* Add an opencv tutorial with the global segments regularization over the Chicago map.

* Should we carry out a rigorous performance analysis?
* Add a figure with the inserted orthogonal edges. Or maybe three figures for all steps
  for creating a contour.
* Move more parameters to the named parameters. As a part of this work, undocument all
  unnecessary overloads and mention instead that GeomTraits e.g. can be ommitted if the kernel
  can be deduced from the input value type.
* In case the neighbor query class and the regularization function MUST get the same set of segments, the neighbor query class could be generated internally.