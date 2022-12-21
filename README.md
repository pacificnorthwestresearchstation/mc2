This is the source code for MC2 dynamic global vegetation model, revision 124.

The source code files comprise all the .h and .cpp files in this directory and
in the century subdirectory. To compile, run make in the century subdirectory
first, which produces a library file name "libcentmc2.a". Then run make in 
this directory, which produces an executable named "mc2".

MC2 is MC1 dynamic global vegetation model (DGVM) re-written in C++ for improved 
code organization and execution efficiency. The century portion of the model
remains in Fortran. MC1 (and therefore MC2) model design are most directly described 
in PNW GTR 508 and PNW GTR 926. Many publications describe its features and 
applications to various research projects. Please see docs/BIBLIOGRAPHY for a 
partial list of publications.

MC1 and MC2 were created, used and maintained by many scientists, with support
from the USDA Forest Service Pacific Northwest Research Station and other
organizations.
