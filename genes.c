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
int totalGenesIn = 0, totalGenes = 0;
int totalProbabilities = 0;
bool genesInitialized = false;
bool ordersCreated = false;
bool validProbs = true, infered = false, loaded = false;

int factorial(int n);
void fillGenesIn(int total);
void fillProbabilityMatrix(int totalGenes);
bool isValidTrio(float prob1, float prob2, float prob3);
Array_float getdistancesInBetween(Array_int order, Array_float *distances, bool alreadyIn);
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
  ordersCreated = false;
  if(total < totalGenes){
    for(int i = total; i < totalGenes;i++){
      if(genesIn.data[i]){
        genesIn.data[i] = 0;
        totalGenesIn--;
      }
    }
    //printf("tot:%d %d\n", total, totalGenes);
    for(int i = 0; i < totalGenes; i++){
      for(int j = i; j < totalGenes; j++){
        if(i >= total || j >= total){
          if(i == j){
            probabilities.data[i].data[j] = 0.0;
          } else {
            probabilities.data[i].data[j] = 2.0;
            probabilities.data[j].data[i] = 2.0;
          }
        }
      }
    }
  }
  totalProbabilities = totalGenes = total;
  checkValidityAll();
  createOrders();
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
    if(i >= geneNames.used){
      initArray(name, char, 4);
      if(i+1 < 10){
        sprintf(name.data, "G0%d", i+1);
      } else {
        sprintf(name.data, "G%d", i+1);
      }
      printf("i:%d\n",i);
      appendArrayP(geneNames, Array_char, name);
    }
  }
}

void freeGenes(){
  totalPossibleOrders = 0;
  totalGenesIn = 0;
  totalGenes = 0;
  totalProbabilities = 0;
  genesInitialized = false;
  ordersCreated = false;
  validProbs = true;
  infered = false;
  loaded = false;
  freeArrayP(geneNames);
  freeArrayP(probabilities);
  freeArrayP(orders);
  freeArray(genesIn);
  freeArrayP(distancesInBetween);
}

void createOrders(){
  //printf("%d\n",validProbs);
  if(!ordersCreated && genesInitialized && validProbs){
    //printf("\nPROBABILITIES\n\n");
    // for(int i = 0; i < totalProbabilities; i++){
    //   for(int j = 0; j < totalProbabilities; j++){
    //     printf("%f,", probabilities.data[i].data[j]);
    //   }
    //   printf("\n");
    // }
    bool first = true;
    totalPossibleOrders = 0;
    int totalOrders = factorial(totalGenesIn)/2;
    bool alreadyIn, stop = false, possible;
    int i2 = 0;
    //printf("totalGI: %d\n",totalGenesIn);
    if(totalGenesIn > 0){
      for(int i = 0; i < totalOrders && !stop; i++){
        alreadyIn = false;
        Array_int order;
        Array_float distances;
        if(orders.used > i){
          deleteAllArray(orders.data[i]);
          order = orders.data[i];
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
        possible = false;
        while(!possible){
          if(first){
            first = false;
          } else {
            permutationNoReverse(order);
          }
          
          if(alreadyIn){
            distances = getdistancesInBetween(order, &(distancesInBetween.data[i]), alreadyIn);
          } else {
            int t = order.used-1;
            initArray(distances, float, t);
            distances = getdistancesInBetween(order, &distances, alreadyIn);
          }
          float sum = 0;
          for(int k = 0; k < distances.used; k++){
            sum += distances.data[k];
            if(distances.data[k] < 0 || sum > 0.5){
              possible = false;
              break;
            } else {
              possible = true;
            }
          }
          if(!possible){
            totalOrders--;
            if(i >= totalOrders){
              stop = true;
              break;
            }
          }
        }
        
        if(!stop){
          if(possible){
            totalPossibleOrders++;
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
      bool rowExists = false;
      bool colExists = false;
      for(int i = 0; i < totalProbabilities; i++){
        if(i != row){
          if(isValidProbability(probabilities.data[row].data[i])){
            rowExists = true;
            break;
          }
        }
      }
      for(int i = 0; i < totalProbabilities; i++){
        if(i != col){
          if(isValidProbability(probabilities.data[col].data[i])){
            colExists = true;
            break;
          }
        }
      }
      if(!rowExists){
        if(genesIn.data[row] == 1){
          totalGenesIn--;
        }
        genesIn.data[row] = 0;
      }
      if(!colExists){
        if(genesIn.data[col] == 1){
          totalGenesIn--;
        }
        genesIn.data[col] = 0;
      }
      ordersCreated = false;
      checkValidityAll();
      createOrders();
    } else {
      insertGenes(row, col);
    }
  }
}

float getProbability(int row, int col){
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
  if(i < totalPossibleOrders){
    return distancesInBetween.data[i];
  }
  Array_float a = {NULL, 0, 0};
  return a;
}

bool isValidProbability(float prob){
  return prob > 0.0 && prob < 0.5;
}

Array_float getdistancesInBetween(Array_int order, Array_float *distances, bool alreadyIn){
  int first, second, i, j, first0 =-1, total0s;
  float greatestKnown, maxProb = 0.5, sumOfDistances = 0.0, probability = 0.0;
  bool solution = false, fixed = false;
  // printf("---------------------------------------------------\nAlreadyIN:%d\ndistances:%p\nOrder: \n",alreadyIn, distances->data);
  // for(i = 0; i < probabilities.used; i++){
  //   for(j=0;j< probabilities.used;j++){
  //     printf("%f  ",probabilities.data[i].data[j]);
  //   }
  //   printf("\n");
  // }
  // printf("\n");
  // for(i = 0; i < order.used; i++){
  //   printf("%d,",order.data[i]);
  // }
  // printf("\n");
  deleteAllArray((*distances));
  for(i = 0; i < order.used; i++){
    if(sumOfDistances <= maxProb){
      if(i+1 < order.used){
        first = order.data[i];
        second = order.data[i+1];
        probability = probabilities.data[first].data[second];
        if(isValidProbability(probability)){
          sumOfDistances += probability;
          appendArray((*distances), float, probability);
        } else {
          if(first0 == -1){
            first0 = i;
          }
          appendArray((*distances), float, 0);
        }
      }
    } else {
      distances->data[0] = -2;
      return *distances;
    }
  }
  for(first0 = 0; first0 < distances->used; first0++){
    if(distances->data[first0] == 0.0){
      fixed = true;
      total0s = 0;
      for(i = first0; i >= 0; i--){
        greatestKnown = 0;
        first = order.data[i];
        for(j = first0+1; j < order.used; j++){
          second = order.data[j];
          solution = isValidProbability(probabilities.data[first].data[second]);
          if(solution){
            greatestKnown = probabilities.data[first].data[second];
            for(int k = i; k < j && k < distances->used; k++){
              greatestKnown = greatestKnown - distances->data[k];
              if(greatestKnown < 0){
                distances->data[0] = -1;
                return *distances;
              }
              if(distances->data[k] == 0.0){
                total0s++;
              }
            }
            //printf("PASA 1\n");
            for(int k = i; k < j && k < distances->used; k++){
              if(distances->data[k] == 0.0){
                distances->data[k] = greatestKnown / total0s;
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
    }
    if(distances->data[first0] == 0.0){
      distances->data[0] = -1;
      return *distances;
    }
  }
  for(first0 = 0; first0 < distances->used; first0++){
    for(i = first0; i >= 0; i--){
      first = order.data[i];
      for(j = first0+1; j < order.used; j++){
        second = order.data[j];
        solution = isValidProbability(probabilities.data[first].data[second]);
        //printf("look at: (%d, %d) (%d, %d) s:%d p:%f\n",first,second, i,j,solution, probabilities.data[first].data[second]);
        probability = 0.0;
        for(int k = i; k < j && k < distances->used;k++){
          //printf("%f+",distances->data[k]);
          probability += distances->data[k];
        }
        if((solution && probability > probabilities.data[first].data[second]) || probability >= 0.5){
          distances->data[0] = -1;
          return *distances;
        }
      }
    }
  }
  //printf("---------------------------------------------------\n");
  return *distances;
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
        freeGenes();
      }
      initArrayP(orders, Array_int, factorial(size)/2);
      initArray(genesIn, int, size);
      initArrayP(distancesInBetween, Array_float, size);
      fillGenesIn(size);
      totalProbabilities = totalGenes = size;
      genesInitialized = true;
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