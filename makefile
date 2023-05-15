main: lodepng image
	gcc image_v5.o lodepng.o -o result -lm
lodepng:
	gcc lodepng.c -c lodepng.o
image:
	gcc image_v5.c -c image_v5.o
 
