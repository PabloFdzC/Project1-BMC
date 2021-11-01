#include "permutationNoReverse.h"

// The following code was created from the explanation at:
// https://www.quora.com/How-would-I-generate-lexically-ordered-permutations-of-a-string
// and tweaked so that reverse permutations were omitted
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

void permutationNoReverse(Array_int a){
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