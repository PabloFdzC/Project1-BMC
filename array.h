#include "global.h"

#define Array(TYPE) struct { TYPE* data; int size; int used; }
#define initArray(A, TYPE, I) A.data = malloc(I * sizeof(TYPE)); A.used = 0; A.size = I
#define initArrayP(A, TYPE, I) initArray(A, TYPE, I); for(A.used = 0; A.used < A.size; A.used++){TYPE arr = {NULL,0,0};A.data[A.used]=arr;} A.used=0
#define resizeArray(A, TYPE, S) A.data = realloc(A.data, (A.size+S) * sizeof(TYPE)); A.size += S
#define appendArray(A, TYPE, E) if (A.used == A.size) { A.size *= 2; resizeArray(A, TYPE, 0); } A.data[A.used++] = E
#define appendArrayP(A, TYPE, E) appendArray(A, TYPE, E); for(int _ = A.used+1; _ < A.size; _++){TYPE arr = {NULL,0,0};A.data[_]=arr;}
#define freeArray(A) if(A.data){free(A.data);} A.data = NULL; A.used = A.size = 0//; printf(#A);printf("\n")
#define freeArrayP(A) for(int _ = 0; _ < A.size; _++){ if(A.data[_].data){printf("i:%d\n",_);freeArray(A.data[_]);}} freeArray(A);
#define deleteAllArray(A) A.used = 0

#ifndef ARRAY_TYPES
#define ARRAY_TYPES

typedef Array(char) Array_char;
typedef Array(Array_char) Array_chars;
typedef Array(Array_chars) Array_chars2;
typedef Array(int) Array_int;
typedef Array(Array_int) Array_ints;
typedef Array(float) Array_float;
typedef Array(Array_float) Array_floats;

#endif
/**
 * This code was inspired by the question in: 
 * https://stackoverflow.com/questions/42293192/making-a-dynamic-array-that-accepts-any-type-in-c
 */