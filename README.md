Compiling on Windows:
Compiled using MSVS 2015.  Will add one for 2013 as well in the future.  
You'll need to set up OpenGL, GLEW, and trimesh.lib in the settings of the solution as required libraries.

Compiling on Linux:
In the file EnvSpec.h: You'll need to switch the comment from Defining Windows to Defining Linux.
The included Makefile should generate a run file that will launch the program.  

HypoCycloid Commands:
'1' & '2' - Increment and Decrement the U-Increment
'Left' - Move NurbsMan backwards along the curve
'Right' - Move NurbsMan forward along the curve
'Up' - Increment Order
'Down' - Decrement Order
'B' - Toggles option to draw the non-affine Nurbs Curve
'CTRL+Z' - Removes the last Control Point that was Added.

Mouse Wheel - Zooms Camera in and Out, or adjusts the weight of a control point
      	      if hovering over control point.
Left Click - Adds a Control Point or hold and drag over existing Control Point
     	     to move point around.
Right Click - Rotate around the origin in 3d space.

Code Re-use:
Majority of Code was written by me with some small functions provided by previous instructors for Image Reading (Not used in HypoCycloid implementation). Code was repurposed from my Animation/Rendering Framework to use a 2D/Orothogonal View Camera to view 3D images on a 2D plane. Same framework will be used for later modelling assignments.



