* Remove all unnecessary functions.
* Improve the contour base algorithm.
* Install OSQP on all testing platforms.
* When adding the QP_solver concept, I should put it in the Solver_interface and refactor the MixedIntegerConcept, too.
* Check the closed/open contour generating examples. They should avoid an almost collinear vertex at the left top corner.
* test_open_contour does not perform well.
* Add a few tests: does the code change the order of the input segments; if the input segments, an input contour (closed/open)
  are already regularized, what happens when running the code.
* Add a free function (and also in offsets) that merges all segments, which form a collinear group; I should first
  find collinear groups and then compute the central segment of each group. Let me call this function unique_segments().
  Add a figure for that function. Unique segments can later be converted into lines.
* Add in docs: the contours must be counterclockwise oriented.
* Add preserve_order parameter to np and also check if I correctly pass all the parameters to this class when I instatiate it
  and maybe the same for parallel/orthogonal/collinear groups.