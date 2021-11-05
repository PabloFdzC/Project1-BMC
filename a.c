#include "global.h"
#include "array.h"



void swap(Array_int a, int k, int l){
  if(k < a.used && k >= 0 && l < a.used && l >= 0){
    int tmp = a.data[k];
    a.data[k] = a.data[l];
    a.data[l] = tmp;
  }
}

void reverse(Array_int a, int k, int l){
  while (k < l){
    swap(a, k, l);
    l--;
    k++;
  }
}

void permutation(Array_int a){
  do{
    int k = -1, l = -1;
    for(int i = a.used-2; i >= 0; i--){
      if(a.data[i] < a.data[i+1]){
        k = i;
        break;
      }
    }
    for(int i = a.used-1; i >= 0; i--){
      if(a.data[k] < a.data[i]){
        l = i;
        break;
      }
    }
    swap(a, k, l);
    reverse(a, k+1, a.used-1);
  }while(a.data[0] > a.data[a.used-1]);
}

int main(){
  for(){
    Array_float *d;
    initArray(d, float, 2);
    for(int i = 0; i < 100000; i++){
      appendArray(d, float, 2.0);
    }
  }
  freeArray(d);
} 