#include "genes.h"
#include "permutationNoReverse.h"
#define ALLOWEDERROR 0.05
int countApend = 0,countApend2 = 0;

Array_chars geneNames;
Array_floats probabilities;
Array_ints orders;
Array_int genesIn;
Array_floats distancesInBetween;
int totalPossibleOrders = 0;
int totalGenesIn = 0;
int totalGenes = 0;
int totalProbabilities = 0;
bool genesInitialized = false;
bool ordersCreated = false;
bool newDistancesNeeded = false;
bool validProbs = true, infered = false, loaded = false;
int asdf = 0;

int factorial(int n);
void fillGenesIn(int total);
void fillProbabilityMatrix(int totalGenes);
bool orderIsPossible(Array_int order);
bool isValidTrio(float prob1, float prob2, float prob3);
Array_float getdistancesInBetween(Array_int order, Array_float distances);
void checkValidity(int g1, int g2);
void checkValidityAll();

int factorial(int n){
  int fact = 1;
  for(int i = 1; i <= n; i++){
    fact = fact * i;
  }
  return fact;
}

void initializeGenes(int total){
  int totalOrders = factorial(total)/2;
  if(!genesInitialized){
    genesInitialized = true;
    initArrayP(geneNames, Array_char, total);
    initArrayP(probabilities, Array_float, total);
    initArrayP(orders, Array_int, totalOrders);
    initArray(genesIn, int, total);
    initArrayP(distancesInBetween, Array_float, total);
  }
  createGeneNames(total);
  fillGenesIn(total);
  fillProbabilityMatrix(total);
  totalProbabilities = totalGenes = total;
}

void fillGenesIn(int total){
  for(int i = 0; i < total; i++){
    if(genesIn.used <= i){
      appendArray(genesIn, int, 0);
    }
  }
}

void fillProbabilityMatrix(int totalGenes){
  int pSize = probabilities.used;
  for(int i = 0; i < totalGenes; i++){
    Array_float probs = {NULL, 0, 0};
    if(genesInitialized && pSize > i){
      probs = probabilities.data[i];
      if(totalGenes > probabilities.data[i].used){
        for(int j = probabilities.data[i].used; j < totalGenes; j++){
          if(i == j){
            appendArray(probs, float, 0.0);
          } else {
            appendArray(probs, float, 2.0);
          }
        }
      }
    } else {
      initArray(probs, float, totalGenes);
      for(int j = 0; j < totalGenes; j++){
        if(i == j){
          appendArray(probs, float, 0.0);
        } else {
          appendArray(probs, float, 2.0);
        }
      }
    }
    if(!genesInitialized || pSize <= i){
      countApend++;
      appendArrayP(probabilities, Array_float, probs);
    } else {
      probabilities.data[i] = probs;
    }
  }
}

void createGeneNames(int totalGenes){
  for(int i = 0; i < totalGenes; i++){
    Array_char name;
    if(geneNames.used > i || (geneNames.size > i && geneNames.data[i].size != 0)){
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
  //printf("\napend1:%d, append2:%d\n\n",countApend,countApend2);
  freeArrayP(geneNames);
  freeArrayP(probabilities);
  freeArrayP(orders);
  freeArray(genesIn);
  freeArrayP(distancesInBetween);
}

void createOrders(){
  if(!ordersCreated && genesInitialized && validProbs){
    printf("\nPROBABILITIES\n\n");
    for(int i = 0; i < totalProbabilities; i++){
      for(int j = 0; j < totalProbabilities; j++){
        printf("%f,", probabilities.data[i].data[j]);
      }
      printf("\n");
    }
    bool first = true;
    newDistancesNeeded = true;
    totalPossibleOrders = 0;
    int totalOrders = factorial(totalGenesIn)/2;
    bool alreadyIn, stop = false, possible;
    int i2 = 0;
    for(int i = 0; i < totalOrders && !stop; i++){
      alreadyIn = false;
      Array_int order;
      Array_float distances;
      //printf("os:%d\n", orders.size);
      if(orders.used > i){
        //printf("GOOD:%d\n", orders.used);
        deleteAllArray(orders.data[i]);
        order = orders.data[i];
        alreadyIn = true;
      } else {
        //printf("BAD:%d\n", orders.used);
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
      possible = false;
      while(!possible){
        //printf("i: %d, odsU: %d\n",i, orders.used);
        //printf("ORDERU: %d, %p, %d\n",order.used, order.data, totalGenesIn);
        if(first){
          first = false;
        } else {
          permutationNoReverse(order);
        }
        if(alreadyIn){
          //printf("AQUI i:%d tO:%d dibU:%d dibS:%d\n", i, totalOrders, distancesInBetween.used, distancesInBetween.size);
          distances = getdistancesInBetween(order, distancesInBetween.data[i]);
        } else {
          int t = order.used-1;
          initArray(distances, float, t);
          distances = getdistancesInBetween(order, distances);
        }
        float sum = 0;
        for(int k = 0; k < distances.used; k++){
          //printf("PASA 0.1.1 k:%d\n",k);
          sum += distances.data[k];
          if(distances.data[k] < 0 || sum > 0.5){
            possible = false;
            break;
          } else {
            possible = true;
          }
        }
        //printf("PASA 0.1.2\n");
        if(!possible){
          totalOrders--;
          if(i >= totalOrders){
            stop = true;
            //printf("PASA 0.1.2.3\n");
            break;
          }
        }
        //printf("PASA 0.1.3 i:%d tO:%d\n", i, totalOrders);
      }
      
      if(!stop){
        if(possible){
          totalPossibleOrders++;
          //printf("PASA 0.1.3.1 dibu:%d dibs:%d\n", distancesInBetween.used,distancesInBetween.size);
          if(alreadyIn){
            distancesInBetween.data[i] = distances;
            orders.data[i] = order;
          } else {
            appendArrayP(orders, Array_int, order);
            appendArrayP(distancesInBetween, Array_float, distances);
          }
        }
      } else {
        if(!alreadyIn){
          freeArray(order);
          freeArray(distances);
        }
      }
    }
    ordersCreated = true;
  }
}

void insertGenes(int g1, int g2){
  ordersCreated = false;
  bool oneIn = genesIn.data[g1] == 1;
  bool twoIn = genesIn.data[g2] == 1;
  if(!oneIn){
    totalGenesIn++;
    genesIn.data[g1] = 1;
  }
  if(!twoIn){
    totalGenesIn++;
    genesIn.data[g2] = 1;
  }
  if(oneIn && twoIn && !loaded){
    checkValidity(g1, g2);
  }
  if(validProbs && totalGenesIn > 0 && !infered && !loaded){
    ordersCreated = false;
    createOrders();
  }
}

void checkValidity(int g1, int g2){
  float prob1, prob2, prob3;
  bool validTrio;
  for(int i = 0; i < totalProbabilities; i++){
    prob1 = probabilities.data[g1].data[i];
    if(isValidProbability(prob1)){
      prob2 = probabilities.data[i].data[g2];
      if(isValidProbability(prob2)){
        prob3 = probabilities.data[g1].data[g2];
        if(isValidProbability(prob3)){
          validTrio = isValidTrio(prob1, prob2, prob3);
          if(!validTrio){
            validProbs = false;
          }
          if(!validProbs || !validTrio){
            checkValidityAll();
          }
        }
      }
    }
  }
}

void checkValidityAll(){
  float prob1, prob2, prob3;
  bool foundBad = false;
  for(int first = 0; first < totalProbabilities && !foundBad; first++){
    for(int second = first; second < totalProbabilities && !foundBad; second++){
      prob1 = probabilities.data[first].data[second];
      if(isValidProbability(prob1)){
        for(int third = second; third < totalProbabilities && !foundBad; third++){
          prob2 = probabilities.data[second].data[third];
          if(isValidProbability(prob2)){
            prob3 = probabilities.data[first].data[third];
            if(isValidProbability(prob3)){
              foundBad = !isValidTrio(prob1, prob2, prob3);
            }
          }
        }
      }
    }
  }
  if(!foundBad){
    validProbs = true;
  }
}

bool isValidTrio(float prob1, float prob2, float prob3){
  return ((prob3 > prob1 && prob3 > prob2 && 
  ((prob3 >= prob2+prob1 && prob3-ALLOWEDERROR < prob2+prob1) ||
  (prob3 <= prob2+prob1 && prob3+ALLOWEDERROR > prob2+prob1))) || 
          (prob2 > prob1 && prob2 > prob3 && 
  ((prob2 >= prob3+prob1 && prob2-ALLOWEDERROR < prob3+prob1) ||
  (prob2 <= prob3+prob1 && prob2+ALLOWEDERROR > prob3+prob1))) ||
          (prob1 > prob2 && prob1 > prob3 && 
  ((prob1 >= prob3+prob2 && prob1-ALLOWEDERROR < prob3+prob2) ||
  (prob1 <= prob3+prob2 && prob1+ALLOWEDERROR > prob3+prob2)))
  );
}

bool inferProbabilities(){
  float sum, first, second;
  int countInfereds = 0;
  if(validProbs){
    for(int i = 0; i < totalProbabilities; i++){
      for(int j = 0; j < totalProbabilities; j++){
        for(int k = 0; k < totalProbabilities; k++){
          if(i != j && j != k && k != i){
            first = probabilities.data[i].data[k];
            if(isValidProbability(first)){
              second = probabilities.data[k].data[j];
              if(isValidProbability(second)){
                if(!isValidProbability(probabilities.data[i].data[j])){
                  sum = first + second;
                  if(isValidProbability(sum)){
                    if(isValidTrio(first, second, sum)){
                      probabilities.data[i].data[j] = sum;
                      infered = true;
                      insertGenes(i, j);
                      countInfereds++;
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    if(countInfereds > 0){
      ordersCreated = false;
      createOrders();
      return true;
    }
  }
  return false;
}

bool areGenesInitialized(){
  return genesInitialized;
}

void insertProbability(float probability, int row, int col){
  if(row < totalProbabilities && col < totalProbabilities && row >= 0 && col >= 0){
    if(probability <= 0 || probability >= 0.5){
      probabilities.data[row].data[col] = 2.0;
      probabilities.data[col].data[row] = 2.0;
    } else {
      probabilities.data[row].data[col] = probability;
      probabilities.data[col].data[row] = probability;
    }
    if(row > genesIn.size && row > col){
      resizeArray(genesIn, int, row);
    } else if(col > genesIn.size) {
      resizeArray(genesIn, int, row);
    }
    if(probability <= 0 || probability >= 0.5){
      
    } else {
      insertGenes(row, col);
    }
  }
}

bool orderIsPossible(Array_int order){
  int first, second, third;
  float firstWithSecond, firstWithThird;
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

float getProbability(int row, int col){
  printf("r:%d, c:%d, pu:%d\n",row, col, probabilities.used);
  if(row < totalProbabilities && col < totalProbabilities){
    return probabilities.data[row].data[col];
  }
  return -2.0;
}

int getTotalGenes(){
  return totalGenes;
}

Array_char getGeneName(int i){
  if(i < geneNames.used){
    return geneNames.data[i];
  }
  Array_char a = {NULL, 0, 0};
  return a;
}

int getTotalPossibleOrders(){
  return totalPossibleOrders;
}

bool areOrdersCreated(){
  return ordersCreated;
}

Array_ints getOrders(){
  return orders;
}

Array_float getDistances(int i){
  if(i < distancesInBetween.used){
    return distancesInBetween.data[i];
  }
  Array_float a = {NULL, 0, 0};
  return a;
}

bool isValidProbability(float prob){
  return prob > 0.0 && prob < 0.5;
}

Array_float getdistancesInBetween(Array_int order, Array_float distances){
  int first, second, i, j, first0 =-1, total0s;
  float greatestKnown, maxProb = 0.5;
  bool solution = false;
  // printf("---------------------------------------------------\nOrder: ");
  // for(i = 0; i < order.used; i++){
  //   printf("%d,",order.data[i]);
  // }
  // printf("\nDistances: ");
  deleteAllArray(distances);
  for(i = 0; i < order.used; i++){
    if(i+1 < totalProbabilities){
      first = order.data[i];
      second = order.data[i+1];
      if(isValidProbability(probabilities.data[first].data[second])){
        appendArray(distances, float, probabilities.data[first].data[second]);
      } else {
        if(first0 == -1){
          first0 = i;
        }
        appendArray(distances, float, 0);
      }
    }
  }
  // for(i = 0; i < distances.size; i++){
  //   printf("%f,",distances.data[i]);
  // }
  // printf("\n");
  while(first0 != -1){
    total0s = 0;
    maxProb = 0.5;
    j = i = first0;
    for(i = first0; i >= 0; i--){
      greatestKnown = 0;
      first = order.data[i];
      for(j = first0+1; j < order.used; j++){
        second = order.data[j];
        solution = isValidProbability(probabilities.data[first].data[second]);
        //printf("first:%d, second:%d, p:%f, s%d\n", first, second, probabilities.data[first].data[second], solution);
        if(solution){
          greatestKnown = probabilities.data[first].data[second];
          for(int k = i; k < j && k < distances.used; k++){
            greatestKnown = greatestKnown - distances.data[k];
            //printf("gk:%f,k:%d, du:%d\n", greatestKnown,k, distances.used);
            if(greatestKnown < 0){
              distances.data[0] = -1;
              //printf("PASA 0.1\n");
              return distances;
            }
            if(distances.data[k] == 0.0){
              total0s++;
            }
          }
          //printf("PASA 1\n");
          for(int k = i; k < j && k < distances.used; k++){
            if(distances.data[k] == 0.0){
              distances.data[k] = greatestKnown / total0s;
            }
            
          }
          //printf("PASA 2\n");
          break;
        }
      }
      //printf("PASA 3\n");
      if(solution){
        break;
      }
    }
    total0s = 0;
    //printf("PASA 4\n");
    for(int k = i; k < j && k < distances.used; k++){
      if(distances.data[k] == 0.0){
        total0s++;
      } else {
        maxProb -= distances.data[k];
      }
    }
    if(total0s > 0){
      for(int k = i; k < j && k < distances.used; k++){
        if(distances.data[k] == 0.0){
          distances.data[k] = maxProb/total0s;
        }
      }
    }
    first0 = -1;
    for(int k = j; k < distances.used; k++){
      if(distances.data[k] == 0){
        first0 = k;
        break;
      }
    }
  }
  distances.used = order.used-1;
  //printf("---------------------------------------------------\n");
  return distances;
}

void bye0sRight(char num[9]){
    int i = 8;
    char commaOrDot;
    float a;
    sscanf("0.3", "%f", &a);
    if(a != 0.3){
      commaOrDot = ',';
    } else {
      commaOrDot = '.';
    }
    for(; i >= 0; i--){
      if(num[i] == '\0'){
        continue;
      }
      if(num[i] == '0'){
        num[i] = '\0';
      } else if(num[i] == commaOrDot){
        break;
      } else {
        break;
      }
    }
}

bool loadFile(char fileName[]){
  Array_chars geneNamesFile;
  Array_floats probabilitiesFile;
  Array_char buffer;
  initArray(buffer, char, 4);
  bool fileArrays = false;
  bool namesFilled = false, probsFilled = false;
  bool badFormat = false;
  int size = 0, probRow = 0, probCol = 1;
  float probability;

  FILE * file = fopen (fileName,"r");
  if (file != NULL){
    for(char c = fgetc(file); feof(file) == 0; c = fgetc(file)){
      //printf("c:%c, pC:%d, pR:%d, fA:%d, nf:%d, pf:%d\n",c, probCol, probRow, fileArrays, namesFilled, probsFilled);
      if(c == '|' || c == '\n'){
        if(geneNamesFile.used == size){
          namesFilled = true;
        }
        if(namesFilled){
          if(probRow == size){
            probsFilled = true;
          }
        }
        deleteAllArray(buffer);
        continue;
      }
      if(isspace(c)){
        continue;
      }
      if(!fileArrays && isdigit(c)){
        do{
          appendArray(buffer, char, c);
          c = fgetc(file);
        }while(isdigit(c));
        ungetc(c, file);
        appendArray(buffer, char, '\0');
        sscanf(buffer.data, "%d", &size);
        if(size > 0){
          initArray(geneNamesFile, Array_char, size);
          initArray(probabilitiesFile, Array_char, size);
          for(int i = 0; i < size; i++){
            Array_float probs = {NULL, 0, 0};
            initArray(probs, float, size);
            for(int j = 0; j < size; j++){
              if(i == j){
                appendArray(probs, float, 0.0);
              } else {
                appendArray(probs, float, 2.0);
              }
            }
            countApend2++;
            appendArrayP(probabilitiesFile, Array_float, probs);
          }
          fileArrays = true;
        }
      } else if(!namesFilled && isalnum(c)){
        if(!fileArrays || geneNamesFile.used > size || probRow >= size || probCol >= size){
          //printf("MAL 1\n");
          badFormat = true;
        } else {
          do{
            appendArray(buffer, char, c);
            c = fgetc(file);
          }while(isalnum(c));
          ungetc(c, file);
          appendArray(buffer, char, '\0');
          Array_char name;
          initArray(name, char, buffer.used);
          for(int i = 0; i < buffer.used; i++){
            appendArray(name, char, buffer.data[i]);
          }
          appendArrayP(geneNamesFile, Array_char, name);
        }
      } else if(isdigit(c)){
        if(!fileArrays || !namesFilled || geneNamesFile.used != size){
          //printf("MAL 2\n");
          badFormat = true;
        } else {
          do{
            appendArray(buffer, char, c);
            c = fgetc(file);
          }while(isdigit(c) || c == ',');
          ungetc(c, file);
          appendArray(buffer, char, '\0');
          sscanf(buffer.data, "%f", &probability);
          probabilitiesFile.data[probRow].data[probCol] = probability;
          probabilitiesFile.data[probCol].data[probRow] = probability;
          probCol++;
          if(probCol == size){
            probRow++;
            probCol = probRow+1;
          }
        }
      } else {
        //printf("MAL 4\n");
        badFormat = true;
      }
      if(badFormat){
        break;
      }
    }
    fclose(file);
    if(badFormat){
      if(namesFilled){
        freeArrayP(geneNamesFile);
      }
      if(probsFilled){
        freeArrayP(probabilitiesFile);
      }
      return false;
    } else {
      if(areGenesInitialized()){
        freeArrayP(geneNames);
        freeArrayP(probabilities);
      } else {
        initArrayP(orders, Array_int, factorial(size)/2);
        initArray(genesIn, int, size);
        initArrayP(distancesInBetween, Array_float, size);
        fillGenesIn(size);
        totalProbabilities = totalGenes = size;
        genesInitialized = true;
      }
      geneNames = geneNamesFile;
      probabilities = probabilitiesFile;
      for(int r = 0; r < size; r++){
        for(int c = r+1; c < size; c++){
          loaded = true;
          insertGenes(r, c);
        }
      }
      checkValidityAll();
      createOrders();
      loaded = false;
    }
    return true;
  }
  return false;
}

bool saveFile(char fileName[]){
  char strFloat[9];
  char strInt[9];
  FILE * file = fopen(fileName,"w");
  if(file != NULL){
    sprintf(strInt, "%d", totalGenes);
    fputs(strInt, file);
    fputc('\n', file);
    for(int i = 0; i < geneNames.used; i++){
      fputs(geneNames.data[i].data, file);
      if(i != geneNames.used-1){
        fputc('|', file);
      } else {
        fputc('\n', file);
      }
    }
    for(int row = 0; row < totalProbabilities; row++){
      for(int col = row+1; col < totalProbabilities; col++){
        sprintf(strFloat, "%f", probabilities.data[row].data[col]);
        bye0sRight(strFloat);
        fputs(strFloat, file);
        if(col != totalProbabilities-1){
          fputc('|', file);
        } else {
          fputc('\n', file);
        }
      }
    }
    fclose(file);
    return true;
  }
  return false;
}