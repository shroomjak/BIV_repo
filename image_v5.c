#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lodepng.h"

#define HIGH_LEVEL 200
#define LOW_LEVEL 90

unsigned char * load_PNG(const char* filename, unsigned * width, unsigned * height);
void write_PNG(const char* filename, unsigned char* image, unsigned width, unsigned height);

void get_pixel(int i, int j, unsigned char *r, unsigned char *g, unsigned char *b, unsigned char *a, unsigned char* image, int width);
void get_bw(unsigned char * src, unsigned char * dst, unsigned width, unsigned height);
void preparation(unsigned char * src, unsigned width, unsigned height);
void Gauss_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height);
void Sobel_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height);
void color(unsigned char * src, unsigned char * dst, unsigned width, unsigned height);

int main(void){
	const char * filename = "skull.png";
	unsigned w, h;
	unsigned char * input = load_PNG(filename, &w, &h);
	if(!input){
        printf("I can not read the picture from the file %s. Error.\n", filename);
        return -1;
    }
    
    unsigned char * raw_bw = malloc(w*h*sizeof(unsigned char)); 
    unsigned char * filtered_bw = malloc(w*h*sizeof(unsigned char));
    unsigned char * output = malloc(w*h*4*sizeof(unsigned char));
    
    get_bw(input, raw_bw, w, h);
    preparation(raw_bw, w, h);
    //Sobel_filter(raw_bw, filtered_bw, w, h);
    //Gauss_filter(filtered_bw, filtered_bw, w, h);
    Gauss_filter(raw_bw, filtered_bw, w, h);
    color(filtered_bw, output, w, h);
    
    const char * new_image = "skull-modified_gauss_2.png";
    write_PNG(new_image, output, w, h);
	return 0;
}

void color(unsigned char * src, unsigned char * dst, unsigned width, unsigned height){
	for(int i = 1; i < width*height-1; i++){
		dst[4*i] = (src[i]+0.5*src[i-1])/1.5;
		dst[4*i+1] = (80+src[i]);
		dst[4*i+2] = (100+src[i]);
		dst[4*i+3] = 255;
	}
	return;
}


void Gauss_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height){
	double Gauss_mat[3][3] = {{0.09, 0.12, 0.09},
							  {0.12, 0.12, 0.12},
							  {0.09, 0.12, 0.09}};
	double c; 
	for(int i = 1; i < height-1; i++)
		for(int j = 1; j < width-1; j++){
			c = 0.0;
			for(int k = -1; k <= 1; k++){
				for(int w = -1; w <= 1; w++){
					c += Gauss_mat[k+1][w+1]*src[(i+k)*width+(j+w)];
				}
			}
			dst[i*width+j] = c;
		}
	return;
}

void Sobel_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height){
	double x_grad = 0, y_grad = 0;
	double Sobel_x[3][3] = {{-1.0, 0.0, 1.0},
							{-sqrt(2), 0.0, +sqrt(2)},
							{-1.0, 0.0, 1.0}};
	double Sobel_y[3][3] = {{-1.0, -sqrt(2), -1.0},
							{0.0, 0.0, 0.0},
							{1.0, +sqrt(2), 1.0}};
	for(int i = 1; i < height-1; i++)
		for(int j = 1; j < width-1; j++){
			x_grad = 0.0;
			y_grad = 0.0;
			for(int k = -1; k <= 1; k++){
				for(int w = -1; w <= 1; w++){
					x_grad += Sobel_x[k+1][w+1]*src[(i+k)*width+(j+w)];
					y_grad += Sobel_y[k+1][w+1]*src[(i+k)*width+(j+w)];
				}
			}
			dst[i*width+j] = sqrt(x_grad*x_grad+y_grad*y_grad);
		}
	return;
}

void preparation(unsigned char * src, unsigned width, unsigned height){
	for(int i = 0; i < width*height; i++){
		if(src[i] > HIGH_LEVEL)	src[i] = 255;
		if(src[i] < LOW_LEVEL)	src[i] = 0;
	}
	return;
}

void get_bw(unsigned char * src, unsigned char * dst, unsigned width, unsigned height){
	for(int i = 0; i < width*height; i++)
		dst[i] = 0.299*src[4*i]+0.587*src[4*i+1]+0.114*src[4*i+2];
	return;
}

void get_pixel(int i, int j, unsigned char *r, unsigned char *g, unsigned char *b, unsigned char *a, unsigned char* image, int width){
   *r =  image[4*(width*j + i) + 0]; 
   *g =  image[4*(width*j + i) + 1]; 
   *b =  image[4*(width*j + i) + 2]; 
   *a =  image[4*(width*j + i) + 3];  
   return;
} 

unsigned char * load_PNG(const char* filename, unsigned * width, unsigned * height) {
	unsigned char * image = NULL;
	unsigned error = lodepng_decode32_file(&image, width, height, filename);
	if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
	return image;
}

void write_PNG(const char* filename, unsigned char* image, unsigned width, unsigned height) {
	unsigned char * png;
	size_t pngsize;

	unsigned error = lodepng_encode32(&png, &pngsize, image, width, height);
	if(!error) lodepng_save_file(png, pngsize, filename);
	/*if there's an error, display it*/
	if(error) printf("error %u: %s\n", error, lodepng_error_text(error));
	return;
}
