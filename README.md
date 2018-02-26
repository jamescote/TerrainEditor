Compiling on Windows:
Compiled using MSVS 2015.  Will add one for 2013 as well in the future.  
You'll need to set up OpenGL, GLEW, and trimesh.lib in the settings of the solution as required libraries.

Compiling on Linux:
In the file EnvSpec.h: You'll need to switch the comment from Defining Windows to Defining Linux.
The included Makefile should generate a run file that will launch the program.  

HypoCycloid Commands:
'A' - increases the radius of the Small circle by 0.1 units
'Z' - Decreases the radius of the Small circle by 0.1 units
'S' - Increases the radius of the Big Circle by 0.1 units
'X' - Decreases the radius of the Small Circle by 0.1 units
'D' - Increases the Number of Cycles by 1
'C' - Decreases the Number of Cycles by 1. Min = 1
'SPACE' - Pauses/Resumes Animation of HypoCycloid
'R' - Toggles the Animation On/Off. Will also reset animation.

Mouse Wheel - Zoom in and out.
Right Mouse Button + Drag = Pan image around.

Code Re-use:
Majority of Code was written by me with some small functions provided by previous instructors for Image Reading (Not used in HypoCycloid implementation). Code was repurposed from my Animation/Rendering Framework to use a 2D/Orothogonal View Camera to view 3D images on a 2D plane. Same framework will be used for later modelling assignments.

Additional Note:
There is a Boolean modifier in the Camera class that can be flipped to use a 3D camera instead of a 2D camera. The 3D camera is spherical and can be rotated to view the Hypocycloid in 3D.

