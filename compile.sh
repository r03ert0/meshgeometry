# GLUT is deprecated since OS X 10.9, but the warnings are annoying
# Drawing to a GLUT window should be replaced by a frame buffer
# See:
#   https://developer.apple.com/library/mac/documentation/GraphicsImaging/Conceptual/OpenGL-MacProgGuide/opengl_offscreen/opengl_offscreen.html
#   http://www.swiftless.com/tutorials/opengl/framebuffer.html

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     machine=Linux;;
    Darwin*)    machine=Mac;;
    CYGWIN*)    machine=Cygwin;;
    MINGW*)     machine=MinGw;;
    *)          machine="UNKNOWN:${unameOut}"
esac

if [ $machine == 'Mac' ]; then
    echo "Compiling for Mac"
    gcc -Wall -lm -lz -Wno-deprecated meshgeometry.c -o meshgeometry_mac -framework Carbon -framework OpenGL -framework GLUT
fi

if [ $machine == 'Linux' ]; then
    echo "Compiling for Linux"
    gcc -Wall  -lm -lz meshgeometry.c -o meshgeometry_unix -lGL -lGLU -lglut
fi

if [ $machine == 'Cygwin' ]; then
    echo "Compiling for Cygwin"
    gcc -Wall meshgeometry.c -o meshgeometry_win.exe -lopengl32 -lglut32
fi

