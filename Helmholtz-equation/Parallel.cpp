#include "Equation.h"

void	parallelSchemeCrossJacob(double* sol)
{
	double* solPrev;
	double	h2 = h * h;
	double	koeff = 1. / (4. + h2 * k * k);
	int		iter = 0;

	solPrev = new double[n * n];

	//inition conditions
	for (int i = 0; i < n * n; i++)
		sol[i] = 0;
	do
	{
		copyVec(solPrev, sol, n * n);
		iter++;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 1; i < n - 1; i++)
				for (int j = 1; j < n - 1; j++)
					sol[i * n + j] = (h2 * fRight(i * h, j * h) + (solPrev[i * n + j - 1] + solPrev[i * n + j + 1] +
						solPrev[(i - 1) * n + j] + solPrev[(i + 1) * n + j])) * koeff;
		}
	} while (norm(sol, solPrev, n * n) > EPS);
	//	while (iter < 6000);
	std::cout << iter << " iterations in parallel Jacoby\n";
	delete[] solPrev;
}

void	parallelSchemeCrossRedBlack(double* sol)
{
	double* solPrev;
	double	h2 = h * h;
	double	koeff = 1. / (4. + h2 * k * k);
	int		iter = 0;

	solPrev = new double[n * n];

	//inition conditions
	for (int i = 0; i < n * n; i++)
		sol[i] = 0;
	do
	{
		copyVec(solPrev, sol, n * n);
		iter++;
#pragma omp parallel
		{
#pragma omp for
			for (int i = 1; i < n - 1; i++)
			{
				for (int j = (i % 2) + 1; j < n - 1; j += 2)
					sol[i * n + j] = (h2 * fRight(i * h, j * h) + (sol[i * n + j - 1] + solPrev[i * n + j + 1] +
						sol[(i - 1) * n + j] + solPrev[(i + 1) * n + j])) * koeff;
				for (int j = i % 2; j < n - 1; j += 2)
					sol[i * n + j] = (h2 * fRight(i * h, j * h) + (sol[i * n + j - 1] + solPrev[i * n + j + 1] +
						sol[(i - 1) * n + j] + solPrev[(i + 1) * n + j])) * koeff;
			}
		}
	} while (norm(sol, solPrev, n * n) > EPS);
	//	while (iter < 6000);
	std::cout << iter << " iterations in parallel Red Black\n";
	delete[] solPrev;
}