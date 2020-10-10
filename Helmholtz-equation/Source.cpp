#include "Equation.h"

int	main()
{
	double*			sol;
	std::ofstream	fout;
	double			t0, t0P, tJ, tJP, tZ, tRB, tRBP;

	sol = new double[n * n];

	t0 = omp_get_wtime();
	schemeCrossJacob(sol);
	tJ = omp_get_wtime();
	std::cout << "Time Jacoby = " << tJ - t0 << std::endl;

	t0P = omp_get_wtime();
	parallelSchemeCrossJacob(sol);
	tJP = omp_get_wtime();
	std::cout << "Time parallel Jacoby = " << tJP - t0P << std::endl;

	std::cout << "Times with parallel " << (tJ - t0) / (tJP - t0P) << std::endl;

/*	t0 = omp_get_wtime();
	schemeCrossZeid(sol);
	tZ = omp_get_wtime();
	std::cout << "Time Zeidel = " << tZ - t0 << std::endl;
	*/

	t0 = omp_get_wtime();
	schemeCrossRedBlack(sol);
	tRB = omp_get_wtime();
	std::cout << "Time Red and Black = " << tRB - t0 << std::endl;

	t0P = omp_get_wtime();
	parallelSchemeCrossRedBlack(sol);
	tRBP = omp_get_wtime();
	std::cout << "Time parallel Red and Black = " << tRBP - t0P << std::endl;

	std::cout << "Times with parallel " << (tRB - t0) / (tRBP - t0P) << std::endl;

	fout.open("Helmholtz.txt");
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++)
			fout << i << '\t' << j << '\t' << sol[i * n + j] << '\n';
	}
	fout.close();

	delete[] sol;
	return 0;
}