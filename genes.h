#include "array.h"

void initializeGenes();

void createGeneNames(int total);

void freeGenes();

void createOrders();

bool inferProbabilities();

bool areGenesInitialized();

void insertProbability(float probability, int row, int col);

float getProbability(int row, int col);

int getTotalGenes();

Array_char getGeneName(int i);

int getTotalPossibleOrders();

bool areOrdersCreated();

Array_ints getOrders();

Array_float getDistances(int i);

bool isValidProbability(float prob);

void bye0sRight(char num[9]);

bool loadFile(char fileName[]);

bool saveFile(char fileName[]);