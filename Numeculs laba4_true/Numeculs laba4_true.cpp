// Numeculs laba4_true.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <cmath>
#include <vector>


double* Gauss(double** matr, int n, int m) {

	//прямой
	double elem;
	for (int j = 0; j < n; j++) {

		double max = 0;
		int strm = 0;

		for (int i = j; i < n; i++) {
			if (abs(matr[i][j]) > max) {
				max = abs(matr[i][j]);
				strm = i;
			}
		}

		if (max > abs(matr[j][j])) {
			double* sto = matr[j];
			matr[j] = matr[strm];
			matr[strm] = sto;
		}

		elem = matr[j][j];

		for (int i = j; i < m; i++) {
			matr[j][i] /= elem;
		}


		for (int i = j + 1; i < n; i++) {
			elem = matr[i][j];
			for (int k = j; k < m; k++) {
				matr[i][k] -= elem * matr[j][k];
			}
		}

	}

	//обратный
	double* vec = new double[m];
	vec[n - 1] = matr[n - 1][n];

	for (int i = n - 2; i >= 0; i--) {
		vec[i] = matr[i][n];
		for (int j = i + 1; j < n; j++) {
			vec[i] -= matr[i][j] * vec[j];
		}
	}


	delete[] matr;

	return vec;
}

void Print(double* vec, int size) {
	for (int i = 0; i < size; i++) {
		std::cout << vec[i] << "  ";
	}
	std::cout << "\n";
}

void Print(double** matrix, int nrow, int ncol) {
	for (int i = 0; i < nrow; i++) {
		for(int j = 0; j < ncol; j++){
			printf("%10.3f", matrix[i][j]);
		}
		std::cout << "\n";
	}
	std::cout << "\n";
}

int main(void){
	int N = 6;
	int m = 5;

	double* pressure = new double[6] { 0.164, 0.328, 0.656, 0.984, 1.312, 1.640 };
	double* expiration = new double[6] { 0.448, 0.432, 0.421, 0.417, 0.414, 0.412 };

	double* POWERX = new double[2 * m + 1]{};
	for (int i = 0; i < m * 2 + 1; i++) {
		double value = 0.f;
		for (int j = 0; j < N; j++) {
			value += std::pow(pressure[j], i);
		}
		POWERX[i] = value;
	}

	std::cout << "\nPowerx: ";
	Print(POWERX, 2 * m + 1);

	double** SUM41 = new double* [m + 1];
	for (int i = 0; i < m + 1; i++) {
		SUM41[i] = new double[m + 2]{};
	}
	SUM41[0][0] = N;
	for (int l = 0; l < m + 1; l++) {
		for (int j = 0; j < m + 1; j++) {
			if (!l && !j) {
				continue;
			}
			int k = l + j;
			SUM41[l][j] = POWERX[k];
		}
	}
	std::cout << "\nMatrix:\n";
	Print(SUM41, m + 1, m + 2);

	double* PRAW = new double[m + 1]{};
	for (int l = 0; l < m + 1; l++) {
		double value = 0.f;
		for (int i = 0; i < N; i++) {
			value += expiration[i] * std::pow(pressure[i], l);
		}
		PRAW[l] = value;
	}
	std::cout << "\nPRAW:";
	Print(PRAW, m + 1);

	for (int i = 0; i < m + 1; i++) {
		SUM41[i][m + 1] = PRAW[i];
	}
	Print(SUM41, m + 1, m + 2);

	double* SOLUTION = Gauss(SUM41, m + 1, m + 2);
	std::cout << "\SOLUTION: ";
	Print(SOLUTION, m + 1);

	double S = 0.f;
	for (int i = 0; i < N; i++) {
		double value = expiration[i];
		for (int j = 0; j < m + 1; j++) {
			value -= (SOLUTION[j] * std::pow(pressure[i], j));
		}
		value = std::pow(value, 2);
		S += value;
	}

	S *= (1 / (N - m));
	S = std::sqrt(S);
	std::cout << "\nS:" << S << "\n";

	return 0;
}


