Compiling on Linux:
In the file EnvSpec.h: You'll need to switch the comment from Defining Windows to Defining Linux.
The included Makefile should generate a run file that will launch the program.  

Windows: Not tested.

Design: Attempted implementation of this paper: https://pages.cpsc.ucalgary.ca/~brosz/research/terrain/tsbyexample.pdf. Terrain is designed to maintain subdivision information inline with geometry information. Allowing Multiresolution to be achieved. Then, a sub area of that information is designed to be selected and blended with another piece of geometry to combine features. Mainly, the design is to take the subdivision information of a very detailed high resolution terrain and blend that information with a low resolution terrain to have that detail over a piece of it's larger terrain. This allows subdivision of this lower resolution terrain to have a detailed section inherited from the higher resolution subject.

Problem: This implementation isn't working. The issue we ran into is adequate storage of the subdivision information. We were unable to find an algorithm that maintained an onto blending in all cases. After Multiresolution, the subdivsion information is stored in a subset D with the current geometry of the terrain stored in C. If the original terrain before multiresolution was stored in, say, A, then |D| = |C|+1. Basically, these sizes are uneven in the algorithm. This makes it strange for odd and even sizes of A as elements need to be added or removed to maintain |D| + |C| = |A|. We tried multiple implementations of storing this information but were unable to find an algorithm or understand how to accomplish this task elegantly or smoothly. It's something I may want to explore in the future.

Controls:
Mouse used to rotate (with Right Click) and Select a section of a Terrain (with Left Click). Multiresolution achieved with the up and down arrows. Generally, multiresolution works; there may be some odd cases where some graphical issues start to occur. Selecting an area and blending it with the other terrain works in a base case of an even square selection. Tab cycles between the two terrains. 

Code Re-use:
Majority of Code was written by me with some small functions provided by previous instructors for Image Reading.



