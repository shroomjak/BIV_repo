#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "lodepng.h"
#define NO_LOOK -1
#define INF 1e6

typedef struct node{
	int v;
	int rank;
	int max;
	struct node * p;
}node;

unsigned char * load_PNG(const char* filename, unsigned * width, unsigned * height);
void write_PNG(const char* filename, unsigned char* image, unsigned width, unsigned height);
void get_pixel(int i, int j, unsigned *r, unsigned * g, unsigned *b, unsigned *a, unsigned char* image, int width);
int dist(unsigned i, unsigned j, unsigned x, unsigned y, unsigned char* image, int width);

node * make_set(int x);
node * find_set(node * a);
node * union_set(node * a, node * b);

void swap(int *, int *);
void quick_sort(int n, int * arr[n], int, int);
int partition(int n, int * arr[n], int, int);

node ** kruskal(int n, int m, int * inc_mat[n]);
int min(int a, int b){
	return (a>b)?b:a;
}
int max(int a, int b){
	return (a>b)?a:b;
}

void swap(int * a, int * b){
	int tmp = *a;
	*a = *b;
	*b = tmp;
	return;
}

void sort(int n, int * arr[n], int l, int r){
	int key, i, j;
	for(j = l; j <= r-1; j++){
		key = j;
		for(i = j+1;i <= r;i++){				
			if(arr[key][2] > arr[i][2]) key = i;
		}
		swap(arr[j], arr[key]);
		swap(arr[j]+1, arr[key]+1);
		swap(arr[j]+2, arr[key]+2);
	}
	return;
}

void quick_sort(int n, int * arr[n], int l, int r){
	srand(time(NULL));
	if(l < r){
		int index = partition(n, arr, l, r);
		quick_sort(n, arr, l, index-1);
		quick_sort(n, arr, index+1, r); 
	}
	return;
}


int partition(int n, int * arr[n], int l, int r){
	int index = l + rand()%(r-l);
	int value = arr[index][2];
	if(index != r){
		swap(arr[r], arr[index]);
		swap(arr[r]+1, arr[index]+1);
		swap(arr[r]+2, arr[index]+2);
	}
	int i = l;
	int j;
	for(j = l; j < r; j++)		//without last element, because we compare other with it
		if(arr[j][2] <= value){
			swap(arr[i], arr[j]);
			swap(arr[i]+1, arr[j]+1);
			swap(arr[i]+2, arr[j]+2);
			i++;
		}
	swap(arr[i], arr[r]);
	swap(arr[i]+1, arr[r]+1);
	swap(arr[i]+2, arr[r]+2);
	return i;
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
	node ** arr = malloc(n*sizeof(node*));
	for(i = 0; i < n; i++) arr[i] = make_set(i);
	for(i = 0; i < m; i++){
		if(find_set(arr[(int)inc_mat[i][0]]) != find_set(arr[(int)inc_mat[i][1]]) &&
		(min(find_set(arr[(int)inc_mat[i][0]])->max, find_set(arr[(int)inc_mat[i][1]])->max) >= inc_mat[i][2] || 
		max(find_set(arr[(int)inc_mat[i][0]])->max, find_set(arr[(int)inc_mat[i][1]])->max) == INF ))
		//)
		{
			union_set(arr[(int)inc_mat[i][0]], arr[(int)inc_mat[i][1]]);
			find_set(arr[(int)inc_mat[i][0]])->max = find_set(arr[(int)inc_mat[i][1]])->max = inc_mat[i][2];
		}
	}	
	return arr;
}


void get_pixel(int i, int j, unsigned *r, unsigned *g, unsigned *b, unsigned *a, unsigned char* image, int width){
   *r =  image[4*(width*i + j) + 0]; 
   *g =  image[4*(width*i + j) + 1]; 
   *b =  image[4*(width*i + j) + 2]; 
   *a =  image[4*(width*i + j) + 3];  
   return;
} 

int dist(unsigned i, unsigned j, unsigned x, unsigned y, unsigned char* image, int width){
	unsigned r1, r2, g1, g2, b1, b2, a1, a2;
	get_pixel(i, j, &r1, &g1, &b1, &a1, image, width);
	get_pixel(x, y, &r2, &g2, &b2, &a2, image, width);
    return sqrt((r1-r2)*(r1-r2)+(b1-b2)*(b1-b2)+(g1-g2)*(g1-g2)+(a1-a2)*(a1-a2));
}

int main(void) {
    const char * filename = "smile.png";
    unsigned w, h;
    int i, j, k, p;
	unsigned char * picture = load_PNG(filename, &w, &h);
	
    if (picture == NULL){
        printf("I can not read the picture from the file %s. Error.\n", filename);
        return -1;
    }
    int N = 3*w*h; 
    int ** edges = malloc(N*sizeof(int *));
    for(i = 0; i < N; i++)
		edges[i] = malloc(3*sizeof(int));
	
    int n[3][2] = {{-1, 0}, {-1, 1}, {0, 1}};
    p = 0;
    for(i = 0; i < h; i++){
		for(j = 0; j < w; j++){
			for(k = 0; k < 3; k++){
				if(i+n[k][0] >= 0 && i+n[k][0] < h && j+n[k][1] >= 0 && j+n[k][1] < w){
					//if(dist(i, j, i+n[k][0], j+n[k][1], picture, w) < delta){
						edges[p][0] = i*w+j;
						edges[p][1] = (i+n[k][0])*w+(j+n[k][1]);
						edges[p][2] = dist(i, j, i+n[k][0], j+n[k][1], picture, w);
						//printf("%d\n", edges[p][2]);
						p++;
					//}
				}
			}
		}
	}
	N = p;
	edges = realloc(edges, N*sizeof(int*));
	/*
	for(i = 0; i < N; i++){
		printf("%d: %d %d %d", i, edges[i][0], edges[i][1], edges[i][2]);
		printf("\n");
	}*/
	printf("ok\n");
	sort(N, edges, 0, N-1);
	printf("ok\n");
	node ** arr = kruskal(w*h, N, edges);
	printf("ok\n");
    for(i = 0; i < w*h; i++){
		picture[4*i] 	= 255*fabs(sin(find_set(arr[i])->v));
		picture[4*i+1]	= 255*fabs(sin(2*find_set(arr[i])->v));
		picture[4*i+2]	= 255*fabs(sin(3*find_set(arr[i])->v));
		picture[4*i+3]	= 255;
	}
    const char * new_image = "smile-modified_v4.png";
    write_PNG(new_image, picture, w, h);
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


