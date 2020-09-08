# cgal-public-dev
This is the Git repository used to store some on-going work by CGAL developers. The repository that hosts the `master` branch of CGAL is [`CGAL/cgal`](http://github.com/CGAL/cgal).

# Start working
To start working with this repository, you must also set the a remote repository that contains the `master` branch of CGAL:

    git remote add cgal git@github.com:CGAL/cgal.git

Then you can create your own working branch. Say, your name is *rineau* and you will work on a new read/write function for the `Triangulation` package:

    git checkout -b Triangulation-add_input_output-rineau --no-track cgal/master
    
More details on the developement of new features can be found on the [CGAL wiki](https://github.com/CGAL/cgal/wiki/Developing-Features-with-Git).

# My work
Jing Yang's GSoC 2020 Submission

Supervisor: @gdamiand 
Github Page: https://github.com/jingyangcarl
Personal Website: https://jingyangcarl.com/

## Overview
The submission includes two features: clipping plane feature and web-viewer feature.
| [Clipping Plane Project](https://github.com/CGAL/cgal-public-dev/tree/gsoc2020-basic-viewer-jyang) | [Web Viewer Project](https://github.com/CGAL/cgal-public-dev/tree/gsoc2020-web-viewer-jyang) |
| --- | --- |
| [![Clipping Plane Project](https://img.youtube.com/vi/VBLP3gglM4k/0.jpg)](https://www.youtube.com/watch?v=VBLP3gglM4k) | [![Web Viewer Project](https://img.youtube.com/vi/Yis21x23YLU/0.jpg)](https://www.youtube.com/watch?v=Yis21x23YLU) |

## Commit log
### Clipping Plane
* cgal/cgal-public-dev@4e642f1ad4 disable clipping plane in opengl 3.0
* cgal/cgal-public-dev@59911d1676 add compatibility shader
* cgal/cgal-public-dev@d672031c36 adjust logic for color output
* cgal/cgal-public-dev@d5c6c8971b re-declare attribute variables as uniform
* cgal/cgal-public-dev@17e14bb40b update
* cgal/cgal-public-dev@45dc247de9 Update key logic; disable clipping plane in 2D mode; de-activate clipping plane operation when clipping plane is off
* cgal/cgal-public-dev@e7a3e7ae9d update key logic
* cgal/cgal-public-dev@06d1a43fae update
* cgal/cgal-public-dev@92594bc8e4 clean code and added comments
* cgal/cgal-public-dev@682d31b22d add alt+c to toggle clipping plane rendering
* cgal/cgal-public-dev@0c2b6a8212 apply size factor to clipping plane translation
* cgal/cgal-public-dev@5d2178c811 adjust clipping plane size
* cgal/cgal-public-dev@4bcdb0b75f make rendering compatible to mono mode and color mode for faces, edges, and vertices
* cgal/cgal-public-dev@1e00f812fd clean code
* cgal/cgal-public-dev@84b35d3d2c apply clipping mode to vertices rendering
* cgal/cgal-public-dev@3c978c543c disable initial clipping plane and add description when enabled clipping plane
* cgal/cgal-public-dev@85e433ec9d add clipping to wireframe mode
* cgal/cgal-public-dev@0c2c14afb9 add another clipping mode for not rendering faces outside the clipping plane
* cgal/cgal-public-dev@d59ab3fc11 clean code
* cgal/cgal-public-dev@44c2b1358e make the real clipping plane adaptive to the mesh size
* cgal/cgal-public-dev@479155f7d0 add wheel event to control transparency
* cgal/cgal-public-dev@48dbf5677d change key binding to rotate and translate clipping plane and add description
* cgal/cgal-public-dev@8447045ea2 add description
* cgal/cgal-public-dev@bc5c2e2411 add translation
* cgal/cgal-public-dev@775a8fbbe0 debug
* cgal/cgal-public-dev@5b2b022467 add corresponding rotation for clipping plane in fragment shader
* cgal/cgal-public-dev@7fc9885fe0 able to rotate real clipping plane
* cgal/cgal-public-dev@b4ce4af509 draw real clipping plane at xoy plane
* cgal/cgal-public-dev@20fbee6e31 clean comments
* cgal/cgal-public-dev@dd9b1a23b0 change clipping plane from a const variable in shader to an attribute value
* cgal/cgal-public-dev@5549aae482 add d parameter to clipping plane and update shader
* cgal/cgal-public-dev@1b53109e50 try to render clipping plane
* cgal/cgal-public-dev@0a7d98ae1a add menu description
* cgal/cgal-public-dev@3db2c0e809 add key function to toggle clipping plane rendering mode
* cgal/cgal-public-dev@f175052d75 add variable to toggle clipping plane
* cgal/cgal-public-dev@96d602a0df clean up debug code
* cgal/cgal-public-dev@06c0a7f72b debug for transparent surfaces
* cgal/cgal-public-dev@e8f2b9a1c4 define macro instead of using magic numbers
* cgal/cgal-public-dev@81e8d5dddb VAO_MONO_FACES render can be possibly removed, since VAO_COLORED_FACES implements the save functions, need to be discussed
* cgal/cgal-public-dev@df3239df1c simplify shader and add comments
* cgal/cgal-public-dev@c79f1ca405 use lambda function to simplify code
* cgal/cgal-public-dev@5716df9ae0 implement one render pipeline to render both solid and transparent mesh
* cgal/cgal-public-dev@860c79eb72 add rendering mode for solid mesh
* cgal/cgal-public-dev@bbad2d5605 successfully pass rendering mode into shaders
* cgal/cgal-public-dev@e688a55d77 disable redundant face drawing
* cgal/cgal-public-dev@7b052cfed7 render bug
* cgal/cgal-public-dev@a97c7c69b9 discard faces outside the clipping is one possible solution
* cgal/cgal-public-dev@219d7e2af4 Debug for the non-working situation when toggle other render model
* cgal/cgal-public-dev@537525401e modified original alpha channel rendering to enable normal alpha blending
* cgal/cgal-public-dev@560e97d63c clipping plane on ambient  and specular
* cgal/cgal-public-dev@10f6c466a8 Merge remote-tracking branch 'refs/remotes/cgal-public-dev/gsoc2020-basic-viewer-jyang' into gsoc2020-basic-viewer-jyang
* cgal/cgal-public-dev@2207b5105d update render formula for accuracy
* cgal/cgal-public-dev@53c7c5b747 clipping plane feature on alpha channel, but the current alpha blending will be disabled when switch render mode, need to be fixed
* cgal/cgal-public-dev@3b2b648450 clipping plane feature on blue channel only
* cgal/cgal-public-dev@b446a25c61 Update cmake
* cgal/cgal-public-dev@5ccec340ee edit shaders to enable clipping along clipping plane
* cgal/cgal-public-dev@57a68fdf2b add vertical clipping plane
* cgal/cgal-public-dev@e9dfef38ac add add object3D class to clipping plane demo
* cgal/cgal-public-dev@6a2b6145b5 enable alpha blending
* cgal/cgal-public-dev@a0acf94169 enable phong shading with clipping plane on
* cgal/cgal-public-dev@872c483584 split by a virtual clipping plane
* cgal/cgal-public-dev@cd45b9202d Merge branch 'gsoc2020-basic-viewer-jyang' of github.com:CGAL/cgal-public-dev into gsoc2020-basic-viewer-jyang
* cgal/cgal-public-dev@c709018d97 Add phong shading and view matrix
* cgal/cgal-public-dev@ed2489684c Add CMakelists
* cgal/cgal-public-dev@efdabf179a able to render a red cube
* cgal/cgal-public-dev@21c226aa3e clean up
* cgal/cgal-public-dev@b6f0a71377 init a demo for shader
* cgal/cgal-public-dev@d238fc112c update gitignore
* cgal/cgal-public-dev@38f16d827e add gitignore to ingore Visual Studio Code produced files

### Web Viewer

* cgal/cgal-public-dev@31f8a6aed3 update
* cgal/cgal-public-dev@f828cb4c3d clean up
* cgal/cgal-public-dev@0ac93690c8 able to draw edges
* cgal/cgal-public-dev@ecef5589c7 able to render the mesh in any size
* cgal/cgal-public-dev@73b51cf771 try rendering face
* cgal/cgal-public-dev@f2682f04d7 update
* cgal/cgal-public-dev@a41fdc927a able to render larger data
* cgal/cgal-public-dev@aa0e8e3a63 the branch is able to render vertices, but the limitation of tcp packeet is 64K, which is 65535 bytes and need to be considered for larger mesh
* cgal/cgal-public-dev@669333d47d merge socket connection to branch from demo
* cgal/cgal-public-dev@68c38ad2c7 able to render lines
* cgal/cgal-public-dev@142536a33f able to render triangles
* cgal/cgal-public-dev@661854ebb9 clean
* cgal/cgal-public-dev@c136513db2 able to render point cloud from c++
* cgal/cgal-public-dev@811922a385 add reaction to default geometry
* cgal/cgal-public-dev@b5c5d40e5f Express backend is able to receive message from React frontend
* cgal/cgal-public-dev@5dbb73d913 clean code
* cgal/cgal-public-dev@e8ef0cc75d c++ application is able to send message to frontend via backend
* cgal/cgal-public-dev@d5a0d0820b React Frontend is able to receive message from Express backend
* cgal/cgal-public-dev@00ddb7f009 clean code and improve output
* cgal/cgal-public-dev@b11bad646b react frontend is able to connect express backend
* cgal/cgal-public-dev@eaa19d9358 clean
* cgal/cgal-public-dev@f4e7221031 clients is able to recieve the response from the server
* cgal/cgal-public-dev@c15d0d7122 the client is able to send message and the server is able to receive the message
* cgal/cgal-public-dev@ce1e111bad server is able to recognize the connection
* cgal/cgal-public-dev@8c3e5928ce the client app is able to connet the server.
* cgal/cgal-public-dev@1564cce755 add socket.io to the backend
* cgal/cgal-public-dev@bf5c7bbe7a add backend server
* cgal/cgal-public-dev@747b5f4944 small modification
* cgal/cgal-public-dev@23f76ba54d Initialize three widget
* cgal/cgal-public-dev@27cd7644d3 react app initialization
* /cgal-public-dev@6d541e5441 update gitignore
