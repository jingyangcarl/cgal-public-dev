* Recompile docs.
* Spellcheck all files.
* Regenerate contour figures.
* Remove all unnecessary functions.
* I should check how the code compiles where there is not OSQP on the system. I should probably turn off some functions in that case. For example, in the    regularize_segments.h header, I use OSQP as default.
* Add in the user manual in the performance section which solver has been used for testing.

* Install OSQP on all testing platforms.
* Add a 3D version of the QP regularization algorithm.
* Add an opencv tutorial with the global segments regularization over the Chicago map.