#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "lodepng.h"

#define INF 1000000

typedef struct node{
	int v;
	long rank;
	int max;
	struct node * p;
}node;

unsigned char * load_PNG(const char* filename, unsigned * width, unsigned * height);
void write_PNG(const char* filename, unsigned char* image, unsigned width, unsigned height);
void get_pixel(int i, int j, int *r, int * g, int *b, int *a, unsigned char* image, int width);
int dist(int i, int j, int x, int y, unsigned char* image, unsigned width);

node * make_set(int x);
node * find_set(node * a);
node * union_set(node * a, node * b);
node ** kruskal(int n, int m, int * inc_mat[n]);

void sharp_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height);
void Gauss_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height);

void radix_sort(int n, int * arr[n], int, int);	
int get_digit(int, int, int);							
int find_max(int n, int * arr[n], int, int);
int find_min(int n, int * arr[n], int, int);
int get_max_digit(int, int);


int min(int a, int b){
	return (a>b)?b:a;
}
int max(int a, int b){
	return (a>b)?a:b;
}

int find_max(int n, int * arr[n], int l, int r){							
	int res, i;
	res = arr[l][2];
	
	for(i = l+1; i <= r; i++) if(res < arr[i][2]) res = arr[i][2];
	return res;
}

int get_digit(int x, int k, int i){ // k is the base of numeral system
	int j = 1;
	while(i > j++) x /= k;
	return x%k;
}

int get_max_pow(int x, int k){ // k is the base of numeral system
	int i;
	for(i = 0; x /= k; i++);
	return i;
}

void radix_sort(int n, int * arr[n], int l, int r){					
	int i, j, k;
	int N = 10;
	int max = find_max(n, arr, l, r);
	short shift = 9;
	int size = N+shift;
	int ** sort_arr = malloc(n*sizeof(int*));
	for(i = 0; i < n; i++) sort_arr[i] = malloc(3*sizeof(int));

	int q = get_max_pow(max, N)+1;		
	int cnt_digit[size];		
	for(k = 1; k <= q; k++){ 
		for (j = 0; j < size; j++) cnt_digit[j] = 0;		
		for(i = l; i <= r; i++)							
			cnt_digit[get_digit(arr[i][2], N, k)+shift]++;
		for(j = 1; j < size; j++)						 
			cnt_digit[j] += cnt_digit[j-1];					
		for(i = r; i >= l; i--)							
			sort_arr[--cnt_digit[get_digit(arr[i][2], N, k)+shift]] = arr[i];	
		
		for(i = l; i <= r; i++) arr[i] = sort_arr[i];
	}	
	return;	
}

node * make_set(int x){
	node * set = malloc(sizeof(node));
	set->v = x;			//value
	set->p = set;		//points to himself
	set->max = INF;
	set->rank = 0;		//zero rank if there is not any vertex which are point to this one
	return set;
}

node * find_set(node * a){
	if(a->p != a)
		a->p = find_set(a->p);		//compressing path heuristic
	return a->p;	
}

node * union_set (node * a, node * b) {
	a = find_set(a);
    b = find_set(b);
    if(a->rank < b->rank) {
		a->p = b;
        return a;
    }
    else{
		b->p = a;
		if(b->rank == a->rank)
			a->rank++;
		return a;
    }
    return NULL;
}


node ** kruskal(int n, int m, int * inc_mat[n]){
	int i;
	int k = 100;
	node ** arr = malloc(n*sizeof(node*));
	for(i = 0; i < n; i++) arr[i] = make_set(i);
	for(i = 0; i < m; i++){
		int credit = min(find_set(arr[(int)inc_mat[i][0]])->max + (int)(k/(find_set(arr[(int)inc_mat[i][0]])->rank+1)), find_set(arr[(int)inc_mat[i][1]])->max+(int)(k/(find_set(arr[(int)inc_mat[i][1]])->rank+1)));
		//int credit = min(find_set(arr[(int)inc_mat[i][0]])->max, find_set(arr[(int)inc_mat[i][1]])->max);
		if(find_set(arr[(int)inc_mat[i][0]]) != find_set(arr[(int)inc_mat[i][1]]) && (credit > inc_mat[i][2]))
		{
			union_set(arr[(int)inc_mat[i][0]], arr[(int)inc_mat[i][1]]);
			find_set(arr[(int)inc_mat[i][1]])->max = find_set(arr[(int)inc_mat[i][0]])->max = credit - inc_mat[i][2];
		}
	}	
	return arr;
}


void get_pixel(int i, int j, int *r, int *g, int *b, int *a, unsigned char* image, int width){
   *r =  image[4*(width*i + j) + 0]; 
   *g =  image[4*(width*i + j) + 1]; 
   *b =  image[4*(width*i + j) + 2]; 
   *a =  image[4*(width*i + j) + 3];  
   return;
} 

int dist(int i, int j, int x, int y, unsigned char* image, unsigned width){
	int r1, r2, g1, g2, b1, b2, a1, a2;
	get_pixel(i, j, &r1, &g1, &b1, &a1, image, width);
	get_pixel(x, y, &r2, &g2, &b2, &a2, image, width);
    return ((double)((r1-r2)*(r1-r2)+(b1-b2)*(b1-b2)+(g1-g2)*(g1-g2)+(a1-a2)*(a1-a2)));	//high contast
}

void sharp_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height){
	double norm = 0.0;
	int i, j;
	int n = 3;
	double mat[3][3]=	{{0, -1.0, 0},
						{-1.0, 5.0, -1.0},
						{0, -1.0, 0}};	
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			norm += (mat[i][j]);
	double r, g, b, a; 
	for(i = 0; i < height; i++)
		for(j = 0; j < width; j++){
			r = 0.0;
			g = 0.0;
			b = 0.0;
			a = 0.0;
			for(int k = 0; k < 3; k++){
				for(int w = 0; w < 3; w++){
					if((i+k-1) >= 0 && (i+k-1) < height && (j+w-1) >= 0 && (j+w-1) < width){
						r += mat[k][w]*(double)src[4*((i+k-1)*width+(j+w-1))];
						g += mat[k][w]*(double)src[4*((i+k-1)*width+(j+w-1))+1];
						b += mat[k][w]*(double)src[4*((i+k-1)*width+(j+w-1))+2];
						a += mat[k][w]*(double)src[4*((i+k-1)*width+(j+w-1))+3];
					}
				}
			}
			dst[4*(i*width+j)] = r/norm;
			dst[4*(i*width+j)+1] = g/norm;
			dst[4*(i*width+j)+2] = b/norm;
			dst[4*(i*width+j)+3] = a/norm;
		}
	return;
}

void Sobel_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height){
	double x_grad = 0, y_grad = 0;
	double Sobel_x[3][3] = {{-1.0, 0.0, 1.0},
							{-2.0, 0.0, 2.0},
							{-1.0, 0.0, 1.0}};
	double Sobel_y[3][3] = {{-1.0, -2.0, -1.0},
							{0.0, 0.0, 0.0},
							{1.0, 2.0, 1.0}};
	double r, g, b, a;
	for(int i = 1; i < height-1; i++)
		for(int j = 1; j < width-1; j++){
			r = g = b = a = 0.0;
			x_grad = 0.0;
			y_grad = 0.0;
			for(int k = -1; k <= 1; k++)
				for(int w = -1; w <= 1; w++){
					x_grad += Sobel_x[k+1][w+1]*src[4*((i+k)*width+(j+w))];
					y_grad += Sobel_y[k+1][w+1]*src[4*((i+k)*width+(j+w))];
			}r = sqrt(x_grad*x_grad+y_grad*y_grad);
			for(int k = -1; k <= 1; k++)	
				for(int w = -1; w <= 1; w++){
					x_grad += Sobel_x[k+1][w+1]*src[4*((i+k)*width+(j+w))+1];
					y_grad += Sobel_y[k+1][w+1]*src[4*((i+k)*width+(j+w))+1];
			}g = sqrt(x_grad*x_grad+y_grad*y_grad);
			for(int k = -1; k <= 1; k++)
				for(int w = -1; w <= 1; w++){	
					x_grad += Sobel_x[k+1][w+1]*src[4*((i+k)*width+(j+w))+2];
					y_grad += Sobel_y[k+1][w+1]*src[4*((i+k)*width+(j+w))+2];
			}b = sqrt(x_grad*x_grad+y_grad*y_grad);
			for(int k = -1; k <= 1; k++)	
				for(int w = -1; w <= 1; w++){	
				x_grad += Sobel_x[k+1][w+1]*src[4*((i+k)*width+(j+w))+3];
				y_grad += Sobel_y[k+1][w+1]*src[4*((i+k)*width+(j+w))+3];
			}a = sqrt(x_grad*x_grad+y_grad*y_grad);
			dst[4*(i*width+j)] = r;
			dst[4*(i*width+j)+1] = g;
			dst[4*(i*width+j)+2] = b;
			dst[4*(i*width+j)+3] = a;
		}
	return;
}

void Gauss_filter(unsigned char * src, unsigned char * dst, unsigned width, unsigned height){
	double norm = 0.0;
	int i, j;
	int n = 3;
	double mat[3][3]=	{{0.5, 1.0, 0.5},
						 {1.0, 2.5, 1.0},
						 {0.5, 1.0, 0.5}};	
	for(i = 0; i < n; i++)
		for(j = 0; j < n; j++)
			norm += (mat[i][j]);
	double r, g, b, a; 
	for(i = 0; i < height; i++)
		for(j = 0; j < width; j++){
			r = 0.0;
			g = 0.0;
			b = 0.0;
			a = 0.0;
			for(int k = 0; k < 3; k++){
				for(int w = 0; w < 3; w++){
					if((i+k-1) >= 0 && (i+k-1) < height && (j+w-1) >= 0 && (j+w-1) < width){
						r += mat[k][w]*(double)src[4*((i+k-1)*width+(j+w-1))];
						g += mat[k][w]*(double)src[4*((i+k-1)*width+(j+w-1))+1];
						b += mat[k][w]*(double)src[4*((i+k-1)*width+(j+w-1))+2];
						a += mat[k][w]*(double)src[4*((i+k-1)*width+(j+w-1))+3];
					}
				}
			}
			dst[4*(i*width+j)] = r/norm;
			dst[4*(i*width+j)+1] = g/norm;
			dst[4*(i*width+j)+2] = b/norm;
			dst[4*(i*width+j)+3] = a/norm;
		}
	return;
}

/*void median_filter(unsigned char * src, unsigned char * dst, int width, int height){
	int pixels[9][3]
	return;
}*/

int main(void) {
    const char * filename = "skull.png";
    unsigned w, h;
    int i, j, k, p;
	unsigned char * raw = load_PNG(filename, &w, &h);
	unsigned char * filtered = malloc(4*w*h*sizeof(unsigned char));
	unsigned char * colored = malloc(4*w*h*sizeof(unsigned char));
    if (raw == NULL){
        printf("I can not read the picture from the file %s. Error.\n", filename);
        return -1;
    }
    
    //sharp_filter(raw, filtered, w, h);
	
	Gauss_filter(raw, filtered, w, h);
	//filtered = raw;
    //Sobel_filter(raw, filtered, w, h);
    //Gauss_filter(filtered, filtered, w, h);
    //Sobel_filter(filtered, filtered, w, h);
    
    int N = 3*w*h; 
    int ** edges = malloc(N*sizeof(int *));
    for(i = 0; i < N; i++)
		edges[i] = malloc(3*sizeof(int));
	
    int n[8][2] = {{0, 1}, {-1, 0}, {-1, 1}, {1, -1}, {1, 0}, {0, -1}, {1, 1}, {-1, -1}};
    p = 0;
    for(i = 0; i < h; i++){
		for(j = 0; j < w; j++){
			for(k = 0; k < 8; k++){
				if(i+n[k][0] >= 0 && i+n[k][0] < h && j+n[k][1] >= 0 && j+n[k][1] < w){
					//if(dist(i, j, i+n[k][0], j+n[k][1], picture, w) < delta){
						edges[p][0] = i*w+j;
						edges[p][1] = (i+n[k][0])*w+(j+n[k][1]);
						edges[p][2] = dist(i, j, i+n[k][0], j+n[k][1], filtered, w);
						//printf("%lf\n", edges[p][2]);
						p++;
					//}
				}
			}
		}
	}
	N = p;
	edges = realloc(edges, N*sizeof(int*));
	
	printf("ok\n");
	radix_sort(N, edges, 0, N-1);
	printf("ok\n");
	node ** arr = kruskal(w*h, N, edges);
	printf("ok\n");
    for(i = 0; i < h; i++){
		for(j = 0; j < w; j++){ 
			int v = find_set(arr[i*w+j])->v;
			/*
			colored[4*(i*w+j)] 	= filtered[4*v+0];
			colored[4*(i*w+j)+1]	= filtered[4*v+1];
			colored[4*(i*w+j)+2]	= filtered[4*v+2];
			colored[4*(i*w+j)+3]	= filtered[4*v+3];
			*/
			
			colored[4*(i*w+j)] 		= 255*fabs(cos(v));
			colored[4*(i*w+j)+1]	= 255*fabs(sin(2*v+1));
			colored[4*(i*w+j)+2]	= 255*fabs(cos(1.5*v+2));
			colored[4*(i*w+j)+3]	= filtered[4*v+3];
			
		}
	}
	//colored = filtered;
    const char * new_image = "skull-modified_v8.png";
    write_PNG(new_image, colored, w, h);
    return 0;
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


