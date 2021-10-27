#include "global.h"

#define Array(TYPE) struct { TYPE* data; int size; int used; }
#define initArray(A, TYPE, I) A.data = malloc(I * sizeof(TYPE)); A.used = 0; A.size = I
#define initArrayP(A, TYPE, I) A.data = malloc(I * sizeof(TYPE)); A.used = 0; A.size = I; for(int i = 0; i < A.size; i++){TYPE arr = {NULL,0,0};A.data[i]=arr;}
#define insertArray(A, TYPE, E) if (A.used == A.size) { A.size *= 2; A.data = realloc(A.data, A.size * sizeof(TYPE)); } A.data[A.used++] = E
#define freeArray(A) free(A.data); A.data = NULL; A.used = A.size = 0
#define freeArrayP(A) for(int i = 0; i < A.size; i++){ if(A.data[i].data){freeArray(A.data[i]);}} freeArray(A)
#define deleteAllArray(A) A.used = 0

#ifndef ARRAY_TYPES
#define ARRAY_TYPES

typedef Array(char) Array_char;
typedef Array(Array_char) Array_chars;
typedef Array(Array_chars) Array_chars2;
typedef Array(int) Array_int;
typedef Array(double) Array_double;
typedef Array(Array_double) Array_doubles;

#endif
/**
 * This code was inspired by the question in: 
 * https://stackoverflow.com/questions/42293192/making-a-dynamic-array-that-accepts-any-type-in-c
 */