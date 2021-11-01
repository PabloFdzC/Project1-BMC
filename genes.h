#include "array.h"

void initializeGenes();

void createGeneNames(int total);

void freeGenes();

void createOrders();

void inferProbabilities();

bool areGenesInitialized();

void insertProbability(double probability, int row, int col);

int getTotalGenes();

Array_char getGeneName(int i);

int getTotalPossibleOrders();

Array_ints getOrders();

Array_double getRelativeDistances(Array_int order);