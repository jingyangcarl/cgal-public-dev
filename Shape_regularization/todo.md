* Regenerate contour figures.
* Remove all unnecessary functions.
* Improve the contour base algorithm.
* Add benchmark graphs to the user manual.
* Add unique segments section in the user manual with the figure. Mention: Unique segments can later be converted into lines.
* Add a few tests: does the code change the order of the input segments; if the input segments, an input contour (closed/open)
are already regularized, what happens when running the code, test unique segments; improve current tests; contours with 0/1/2 edges.
* Why the first point in the closed contour is slightly different from the first point in the open contour?

* Install OSQP on all testing platforms.
* When adding the QP_solver concept in the next package release, I should put it in the Solver_interface and refactor the MixedIntegerConcept, too.
* Add a 3D version of the QP regularization algorithm.