# GLUT is deprecated since OS X 10.9, but the warnings are annoying
# Drawing to a GLUT window should be replaced by a frame buffer
# See:
#   https://developer.apple.com/library/mac/documentation/GraphicsImaging/Conceptual/OpenGL-MacProgGuide/opengl_offscreen/opengl_offscreen.html
#   http://www.swiftless.com/tutorials/opengl/framebuffer.html

OS=$(uname -s) || exit 1

case $OS in
'Darwin')
	gcc -Wall -lz -Wno-deprecated meshgeometry.c -o meshgeometry_mac -framework Carbon -framework OpenGL -framework GLUT
	;;
'Linux')
	# The '-Wno-unknown-pragmas' is there to suppress warnings about #pragma mark
	gcc -Wall -Wno-unknown-pragmas meshgeometry.c -o meshgeometry_unix -lGL -lGLU -lglut -lm -lz
	;;
*)
	echo "Unknown OS $OS, if it is some kind of unix, try:" 1>&2
	echo "  gcc -Wall meshgeometry.c -o meshgeometry_unix -lGL -lGLU -lglut -lm -lz" 1>&2
esac
