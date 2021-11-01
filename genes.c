#include "genes.h"
#include "permutationNoReverse.h"

Array_chars geneNames;
Array_doubles probabilities;
Array_ints orders;
Array_int genesIn;
Array_double relativeDistances;
int totalPossibleOrders = 0;
int totalGenesIn = 0;
int totalGenes = 0;
bool genesInitialized = false;
bool ordersCreated = false;
bool newDistancesNeeded = false;

int factorial(int n);
void fillGenesIn();
void fillProbabilityMatrix(int totalGenes);
bool orderIsPossible(Array_int order);

int factorial(int n){
  int fact = 1;
  for(int i = 1; i <= n; i++){
    fact = fact * i;
  }
  return fact;
}

void initializeGenes(int total){
  totalGenes = total;
  int totalOrders = factorial(totalGenes)/2;
  if(!genesInitialized){
    genesInitialized = true;
    initArrayP(geneNames, Array_char, totalGenes);
    initArrayP(probabilities, Array_double, totalGenes);
    initArrayP(orders, Array_int, totalOrders);
    initArray(genesIn, int, totalGenes);
    initArray(relativeDistances, double, totalGenes);
  }
  fillGenesIn();
  fillProbabilityMatrix(totalGenes);
}

void fillGenesIn(){
  for(int i = 0; i < genesIn.size; i++){
    appendArray(genesIn, int, 0);
  }
}

void fillProbabilityMatrix(int totalGenes){
  int pSize = probabilities.used;
  for(int i = 0; i < totalGenes; i++){
    Array_double probs;
    if(genesInitialized && pSize > i){
      probs = probabilities.data[i];
      deleteAllArray(probs);
    } else {
      initArray(probs, double, totalGenes);
    }
    for(int j = 0; j < totalGenes; j++){
      if(i == j){
        appendArray(probs, double, 0.0);
      } else {
        appendArray(probs, double, 2.0);
      }
    }
    if(genesInitialized && pSize > i){
      probabilities.data[i] = probs;
    } else {
      appendArrayP(probabilities, Array_double, probs);
    }
  }
}

void createGeneNames(int totalGenes){
  for(int i = 0; i < totalGenes; i++){
    Array_char name;
    if(geneNames.used > i && geneNames.data[i].size != 0){
      continue;
    } else {
      initArray(name, char, 4);
      if(i+1 < 10){
        sprintf(name.data, "G0%d", i+1);
      } else {
        sprintf(name.data, "G%d", i+1);
      }
      appendArrayP(geneNames, Array_char, name);
    }
  }
  geneNames.used = totalGenes;
}

void freeGenes(){
  freeArrayP(geneNames);
  freeArrayP(probabilities);
  freeArrayP(orders);
}

void createOrders(){
  if(!ordersCreated){
    newDistancesNeeded = ordersCreated = true;
    totalPossibleOrders = 0;
    int totalOrders = factorial(totalGenesIn)/2;
    bool alreadyIn, stop = false, possible;
    int i2 = 0;
    for(int i = 0; i < totalOrders; i++){
      alreadyIn = false;
      Array_int order;
      if(orders.data[i].used > 0){
        deleteAllArray(orders.data[i]);
        alreadyIn = true;
      } else {
        initArray(order, int, totalGenes);
      }
      for(int j = 0; j < totalGenesIn; j++){
        if(i == 0){
          while(genesIn.data[i2] == 0){
            i2++;
          }
          appendArray(order, int, i2++);
        } else {
          appendArray(order, int, orders.data[i-1].data[j]);
        }
      }
      if(i > 0){
        possible = false;
        while(!possible){
          permutationNoReverse(order);
          possible = orderIsPossible(order);
          if(!possible){
            totalOrders--;
            if(i >= totalOrders){
              stop = true;
              break;
            }
          }
        }
      }
      if(!stop){
        totalPossibleOrders++;
        if(alreadyIn){
          orders.data[i] = order;
        } else {
          appendArrayP(orders, Array_int, order);
        }
      }
    }
  }
}

void insertGenes(int g1, int g2){
  ordersCreated = false;
  bool notIn1 = genesIn.data[g1] == 0;
  bool notIn2 = genesIn.data[g2] == 0;
  if(notIn1){
    totalGenesIn++;
    genesIn.data[g1] = 1;
  }
  if(notIn2){
    totalGenesIn++;
    genesIn.data[g2] = 1;
  }
}

void inferProbabilities(){
  for(int i = 0; i < probabilities.used; i++){
    for(int j = 0; j < probabilities.used; j++){
      for(int k = 0; k < probabilities.used; k++){
        double sum = probabilities.data[i].data[k] + probabilities.data[k].data[j];
        if(probabilities.data[i].data[j] > sum){
          probabilities.data[i].data[j] = sum;
          insertGenes(i, j);
        }
      }
    }
  }
}

bool areGenesInitialized(){
  return genesInitialized;
}

void insertProbability(double probability, int row, int col){
  probabilities.data[row].data[col] = probability;
  probabilities.data[col].data[row] = probability;
  if(row > genesIn.size && row > col){
    resizeArray(genesIn, int, row);
  } else if(col > genesIn.size) {
    resizeArray(genesIn, int, row);
  }
  insertGenes(row, col);
}

bool orderIsPossible(Array_int order){
  int first, second, third;
  double firstWithSecond, firstWithThird;
  bool possible = true;
  if(order.used > 2){
    for(int i = 0; i < order.used-2; i++){
      first = order.data[i];
      for(int j = i+1; j < order.used-1; j++){
        second = order.data[j];
        firstWithSecond = probabilities.data[first].data[second];
        for(int k = j+1; k < order.used; k++){
          third = order.data[k];
          firstWithThird = probabilities.data[first].data[third];
          if(firstWithSecond < 0.5 && firstWithThird < 0.5 && firstWithSecond > 0.0 &&
            firstWithThird > 0.0 && firstWithSecond >= firstWithThird){
            return possible = false;
          }
        }
      }
      
    }
  }
  return possible;
}

int getTotalGenes(){
  return geneNames.used;
}

Array_char getGeneName(int i){
  return geneNames.data[i];
}

int getTotalPossibleOrders(){
  return totalPossibleOrders;
}

Array_ints getOrders(){
  return orders;
}

Array_double getRelativeDistances(Array_int order){
  int first, second, i, j, between;
  double borders, left, right;
  for(int r = 0; r < probabilities.used; r++){
    for(int c = 0; c < probabilities.used; c++){
      printf("%lf,",probabilities.data[r].data[c]);
    }
    printf("\n");
  }
  for(i = 0; i < relativeDistances.size; i++){
    relativeDistances.data[i] = 0;
  }
  for(i = 0; i < order.used-1; i++){
    first = order.data[i];
    between = 0;
    for(j = i+1; j < order.used; j++){
      second = order.data[j];
      between++;
      if(probabilities.data[first].data[second] > 0.0 && probabilities.data[first].data[second] < 0.5){
        break;
      }
    }
    for(int k = i; k < j; k++){
      relativeDistances.data[k] = probabilities.data[i].data[j]/between;
    }
    if(between == 1 && i > 0){
      borders = left = right = 0.0;
      if(i+2 < order.used){
        first = order.data[i-1];
        second = order.data[i+2];
        if(probabilities.data[first].data[second] > 0.0 && probabilities.data[first].data[second] < 0.5){
          borders = probabilities.data[first].data[second];
        }
        first = order.data[i+1];
        second = order.data[i+2];
        if(probabilities.data[first].data[second] > 0.0 && probabilities.data[first].data[second] < 0.5){
          right = probabilities.data[first].data[second];
        }
      }
      first = order.data[i-1];
      second = order.data[i];
      if(probabilities.data[first].data[second] > 0.0 && probabilities.data[first].data[second] < 0.5){
        left = probabilities.data[first].data[second];
      }
      if(borders > 0.0){
        if(left > 0.0){
          if(right == 0.0){
            relativeDistances.data[i+1] = borders - left;
          }
        } else if (right > 0.0){
          relativeDistances.data[i-1] = borders - right;
          relativeDistances.data[i+1] = right;
        }
      }
      
    }
  }
  relativeDistances.used = order.used-1;
  
  return relativeDistances;
}