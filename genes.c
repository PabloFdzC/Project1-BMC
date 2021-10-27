#include "genes.h"

Array_chars geneNames;
Array_doubles probabilities;

void initializeGenes(){
  initArrayP(geneNames, Array_char, 5);
  initArrayP(probabilities, Array_double, 5);
}

void createGeneNames(int total){
  for(int i = 0; i <= total; i++){
    Array_char name;
    if(geneNames.data[i].data != NULL && geneNames.size < i){
      continue;
    } else {
      initArray(name, char, 4);
      if(i < 10){
        sprintf(name.data, "G0%d", i);
      } else {
        sprintf(name.data, "G%d", i);
      }
      insertArray(geneNames, Array_char, name);
    }
  }
  geneNames.used = total;
}

int getTotalGenes(){
  return geneNames.used;
}

Array_char getGeneName(int i){
  return geneNames.data[i];
}

void freeGenes(){
  freeArrayP(geneNames);
  freeArrayP(probabilities);
}