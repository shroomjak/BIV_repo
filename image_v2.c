#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "lodepng.h"
#define NO_LOOK -1

typedef struct Graph {
	int ** V;
}Graph;

void dfs(Graph * G, int N, int col[N], int s, int color);
unsigned char * load_PNG(const char* filename, unsigned * width, unsigned * height);
void write_PNG(const char* filename, unsigned char* image, unsigned width, unsigned height);
void get_pixel(int x, int y, int *r, int *g, int *b, int *a, unsigned char* image, int width) ;
int is_close(int r1, int g1, int b1, int a1, int r2, int g2, int b2, int a2);
Graph* init_graph(int N);
void add_edge(Graph *G, int i, int j , int x, int y, int width);


void dfs(Graph * G, int N, int col[N], int s, int color){
	int i;
	col[s] = 1;
	for(i = 1; i <= G->V[s][1]; i++)
		if(col[G->V[s][i]] == -1)
			dfs(G, N, col, G->V[s][i], color);
	col[s] = color;
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

void get_pixel(int i, int j, int *r, int *g, int *b, int *a, unsigned char* image, int width) 
{
   *r =  image[4*(width*j + i) + 0]; 
   *g =  image[4*(width*j + i) + 1]; 
   *b =  image[4*(width*j + i) + 2]; 
   *a =  image[4*(width*j + i) + 3];  
   return;
}

int is_close(int r1, int g1, int b1, int a1, int r2, int g2, int b2, int a2) {
    int e_r = 10, e_g = 10, e_b = 10, e_a = 255;
    return fabs(r1 - r2) < e_r &&  fabs(g1 - g2) < e_g &&  fabs(b1 - b2) < e_b &&  fabs(a1 - a2) < e_a;
 }


Graph* init_graph(int N) {
	Graph * G = malloc(sizeof(Graph));
	G->V = malloc(N*sizeof(int*));
	//if(!G || !G->list) printf("Error");
	for (int i = 0 ; i < N; i++){
		G->V[i] = malloc(2*sizeof(int));
		(G->V[i])[1] = 0;
		(G->V[i])[0] = NO_LOOK;
	}
	return G;
}


void add_edge(Graph * G, int i, int j , int x, int y, int width) {
	int n = j*width + i;
	int m = y*width + x;
	G->V[n] = realloc(G->V[n], ((G->V[n])[1]+3)*sizeof(int)); 
	G->V[m] = realloc(G->V[m], ((G->V[m])[1]+3)*sizeof(int));
	
	(G->V[n])[1]++; (G->V[m])[1]++;
	(G->V[n])[(G->V[n])[1]+1] = m;
	(G->V[m])[(G->V[m])[1]+1] = n; 
	
	return;
}

int main(void) {
    const char * filename = "watermelon.png";
    unsigned w, h;
    int r, g, b, a;
    int r1, g1, b1, a1;
    int i, j;
    unsigned colors = 2;
	unsigned char * picture = load_PNG(filename, &w, &h);
	/*for(i = 0; i < 4*w*h; i++){
		printf("r:%d g:%d b:%d a:%d\n", picture[i], picture[i+1], picture[i+2], picture[i+3]);
	}*/
    if (picture == NULL){
        printf("I can not read the picture from the file %s. Error.\n", filename);
        return -1;
    }
    Graph * G = init_graph(w*h);
    if(!G) {
      printf("Can not allocate memory for Graph\n");
      return -1;
    }


    // read file and convert it to 2D array
        // function get_pixel is simple
    for (i = 1; i < w-1; i++){
        for (j = 1; j < h-1; j++){
            get_pixel(i, j, &r, &g, &b, &a, picture, w);
            
            get_pixel(i-1, j, &r1, &g1, &b1, &a1, picture, w );
            if (is_close(r,  g,  b,  a,
                         r1, g1, b1, a1  )) {
                             add_edge(G, i, j , i-1, j, w);
                         }

            get_pixel(i+1, j, &r1, &g1, &b1, &a1, picture, w );
            if (is_close(r,  g,  b,  a,
                         r1, g1, b1, a1  )) {
                             add_edge(G, i, j, i+1, j, w);
                         }
            get_pixel(i, j-1, &r1, &g1, &b1, &a1, picture, w );
            if (is_close(r,  g,  b,  a,
                         r1, g1, b1, a1)) {
                        add_edge(G, i, j, i, j-1, w);
            }
            get_pixel(i, j+1, &r1, &g1, &b1, &a1, picture, w );
            if (is_close(r,  g,  b,  a,
                         r1, g1, b1, a1 )) {
                        add_edge(G, i, j, i, j+1, w);
            }
			get_pixel(i+1, j+1, &r1, &g1, &b1, &a1, picture, w );
            if (is_close(r,  g,  b,  a,
                         r1, g1, b1, a1 )) {
                        add_edge(G, i, j, i+1, j+1, w);
            }
            get_pixel(i-1, j-1, &r1, &g1, &b1, &a1, picture, w );
            if (is_close(r,  g,  b,  a,
                         r1, g1, b1, a1 )) {
                        add_edge(G, i, j, i-1, j-1, w);
            }
            get_pixel(i+1, j-1, &r1, &g1, &b1, &a1, picture, w );
            if (is_close(r,  g,  b,  a,
                         r1, g1, b1, a1 )) {
                        add_edge(G, i, j, i+1, j-1, w);
            }
            get_pixel(i-1, j+1, &r1, &g1, &b1, &a1, picture, w );
            if (is_close(r,  g,  b,  a,
                         r1, g1, b1, a1 )) {
                        add_edge(G, i, j, i-1, j+1, w);
            }

        }
    }

    // analyze 2D array
        // use graph connectivity algorithm for connectivity areas  
	int col[w*h];
    printf("w:%d h:%d\n", (int)w, (int)h);
    for(i = 0; i < w*h; i++) col[i] = -1;
	for(i = 0; i < w*h; i++){
		if(col[i] == -1){
			dfs(G, w*h, col, i, colors);
			colors++;
		}
	}
	printf("\ncolors: %d\n", colors); 
	for(i = 0; i < 4*w*h; i+=4){
		j = (col[i/4]*255.0/(double)(colors));
		//printf("%d %d %d\n", i/4, j, col[i/4]);
		picture[i] = j;
		picture[i+1] = 0;
		picture[i+2] = 0;
		picture[i+3] = 255;
	}
    // convert 2D array to file and write it
    char * new_image = "skull-modified.png";
    write_PNG(new_image, picture, w, h);
    return 0;
}


