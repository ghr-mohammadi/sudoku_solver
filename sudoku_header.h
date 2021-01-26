#ifndef SUDOKU_H
#define SUDOKU_H

#include <iostream>
#include <cstdlib>
#include <random>
#include <ctime>
#include <thread>

#include "numbers.h"

#define N2 81
#define N 9
#define n 3
#define tir 5
#define NP 500 // number of population => must be even
#define sNP 0.95 // saved number of population => 0 < sNP < 1
#define phi 10
#define deghat 0.88

int intRand(const int & min, const int & max);

void printFunc(int *t);

int compare_(const void * a, const void * b);

int init(int *t, int *t_);

void initFunc(int *t, int *t_);

bool numTest(char ch);

bool lineTest(char * lchar, int * lint);

void getFunc(int *t);

int fitness(int *t);

//int compare(const void * a, const void * b);

void preSolve(int *t);

int *solveFunc(int *t);

#endif //SUDOKU_H
