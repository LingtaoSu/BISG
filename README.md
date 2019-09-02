# BISG
A method for cancer survival related gene sets detection.

BISG based on the adaptation of the RFN based biclustering methods, which is an unsupervised technique that learns a non-linear, high-dimensional representation of its input. The underlying algorithm has been submitted to AAAI20.


librfn is implemented in C++ and can be easily integrated in existing code bases. 
Installation:
Adjust the Makefile to your needs
Type make to start the building process

Requirements:
To run the code, you require a CUDA 7.5 (or higher) compatible GPU. While in theory CUDA 7.0 is also supported.

Note that BISG makes heavy use of BLAS and LAPACK, so make sure to link it to a high-quality implementation to get optimal speed (e.g. OpenBLAS or MKL) by modifying the Makefile.

Implementation Note
The BISG algorithm is based on the EM algorithm. Within the E-step, the algorithm includes a projection procedure that can be implemented in several ways. To make sure no optimzation constraints are violated during this projection, the original publication tries the simplest method first, but backs out to more and more complicated updates if easier method fail. 
