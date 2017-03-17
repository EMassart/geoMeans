# geoMeans
Means and consensus on manifolds - Estimation of the Karcher mean of positive definite matrices

Codes used for my master's thesis.

The folder Algorithms_geomeans contains the implementations of the different algorithms considered, enabling to compute or estimate geometric means.

The folder Basis contains inplementations of the different internal functions required by the algorithms, for example the function evaluating the intrinsic distance between two matrices or the one computing the geometric mean of two matrices. It contains also the implementation used to compute the Karcher mean by a steepest descent algorithm.
Most of the functions located in this subfolder are extracted from existing toolboxes, and the references are mentionned as the first lines of the functions.

The folder Input contains two implementations generating symmetric positive definite matrices.

Finally, the folder Tests contains a script comparing the performances of all the algorithms. The execution of this script is quite long (it may take several hours).
