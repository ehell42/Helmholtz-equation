#pragma once
#include <iostream>
#include <fstream>
#include "omp.h"

const double	PI = 3.14159265358979323846;
const double	EPS = 1.e-11;
const double	a = 0.;
const double	b = 1.;
const int		n = 100;
const double	h = (b - a) / n;
const double	k = 1. / (h * h);

double	fRight(double x, double y);

void	copyVec(double* A, double* B, int m);
double	norm(double* A, double* B, int m);

void	schemeCrossJacob(double* sol);
void	schemeCrossZeid(double* sol);
void	schemeCrossRedBlack(double* sol);

void	parallelSchemeCrossJacob(double* sol);
void	parallelSchemeCrossRedBlack(double* sol);