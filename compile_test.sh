# GLUT is deprecated since OS X 10.9, but the warnings are annoying
# Drawing to a GLUT window should be replaced by a frame buffer
# See:
#   https://developer.apple.com/library/mac/documentation/GraphicsImaging/Conceptual/OpenGL-MacProgGuide/opengl_offscreen/opengl_offscreen.html
#   http://www.swiftless.com/tutorials/opengl/framebuffer.html

OS="$(uname -s)" || exit 1

case "${OS}" in
'Darwin')
    echo "Compiling tests for Mac"
    g++ -std=c++11 -Wall -lz -Wno-deprecated meshgeometry_test.cpp -o meshgeometry_test_mac -framework Carbon -framework OpenGL -framework GLUT
    ;;
*)
    echo "Unknown OS $OS, if it is some kind of unix, try:" 1>&2
    echo "  gcc -Wall meshgeometry.c -o meshgeometry_unix -lGL -lGLU -lglut -lm -lz" 1>&2
    ;;
esac
