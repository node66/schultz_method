#include<iostream>
#include<fstream>
#include<vector>
#include <iomanip>

using namespace std;

const double EPS = 0.01;

struct csr_matrix
{
	unsigned n, nnz;
	unsigned *ANC, *ANL;
	double *AV;

}matrix_csr;


void Read(csr_matrix & matrix_csr); //Считывание матрицы из файла 
inline void print(vector<vector<double>> s); //Вывод матрицы в консоль
inline double norm(vector<vector<double>> s); //Возвращает норму матрицы
void find_psi(csr_matrix matrix_csr, vector<vector<double>> U, vector<vector<double>> &psi);//Считает матрицу psi
void find_U(csr_matrix  matrix_csr, vector<vector<double>> &U, vector<vector<double>> psi);//Считаем матрицу U

int main() {
	setlocale(0, "");
	/*Считываем матрицу*/
	Read(matrix_csr);
	/*Задаем матрицу начального приближения*/
	vector<vector<double>> U(matrix_csr.n, vector<double>(matrix_csr.n, 0.0));
	U[0][0] = 0.6; U[0][1] = -0.5; U[0][2] = -0.4;
	U[1][0] = 0.1; U[1][1] = 0.6;  U[1][2] = 0.2;
	U[2][0] = 0.1; U[2][1] = -0.5; U[2][2] = 0.5;

	/*Инициализируем матрицу psi*/
	vector<vector<double>> psi(matrix_csr.n, vector<double>(matrix_csr.n, 0.0));
	unsigned i(0);
	find_psi(matrix_csr, U, psi);
	/*Сам итерационный процесс*/
	if (norm(psi) < 1) {
		do
		{
			find_psi(matrix_csr, U, psi);
			find_U(matrix_csr, U, psi);
		} while (norm(psi) > EPS && ++i);
	}
	else
		cerr << " Условие теоремы не выполнено" << endl;
	cout << "Обратная матрица = " << endl;
	print(U);
	cout << "Проверка " << endl;
	/*Умножение найденой обратной матрицы на исходную матрицу */
	vector<vector<double>> tmp(matrix_csr.n, vector<double>(matrix_csr.n, 0.0));
	for (size_t k = 0; k < matrix_csr.n; k++)
	{
		for (size_t i = 0; i < matrix_csr.n; i++)
		{
			double f = 0.0;
			const size_t lb = matrix_csr.ANL[i];
			const size_t ub = matrix_csr.ANL[i + 1];

			for (size_t j = lb; j < ub; j++)
				f += matrix_csr.AV[j] * U[matrix_csr.ANC[j]][k];
			tmp[i][k] += f;
		}
	}
	print(tmp);
	delete[] matrix_csr.ANC;
	delete[] matrix_csr.ANL;
	delete[] matrix_csr.AV;
	system("pause");
	return 0;
}

void Read(csr_matrix & matrix_csr)
{
	double *AV1;
	unsigned int  *ANR1, *ANC1;

	ifstream data;
	data.open("m.mtx");
	char s[255];
	do
	{
		data.getline(s, 255);
	} while ((s[0] == '%') && (!data.eof()));

	unsigned size[3];
	size_t k = 0;
	for (size_t i = 0; i < 3; i++) {
		size[i] = 0;
		while ((s[k] != ' ') && (s[k] != '\0'))
		{
			size[i] = size[i] * 10 + s[k] - 0x30;
			k++;
		}
		k++;
	};

	matrix_csr.n = size[0];
	matrix_csr.nnz = size[2];
	AV1 = new double[matrix_csr.nnz];
	ANC1 = new unsigned[matrix_csr.nnz];
	ANR1 = new unsigned[matrix_csr.nnz];

	k = 0;
	while ((k < matrix_csr.nnz) && (!data.eof()))
	{
		double val;
		unsigned val_r, val_c;
		data >> val_r >> val_c >> val;
		AV1[k] = val;
		ANR1[k] = val_r - 1;
		ANC1[k] = val_c - 1;
		k++;
	}
	data.close();
	matrix_csr.AV = AV1;
	matrix_csr.ANL = ANR1;
	matrix_csr.ANC = ANC1;
	/*Соритируем массивы*/
	for (size_t i = 0; i < matrix_csr.nnz; i++)
		for (size_t j = i + 1; j < matrix_csr.nnz; j++)
			if (matrix_csr.ANL[j] < matrix_csr.ANL[i]) {
				swap(matrix_csr.ANL[i], matrix_csr.ANL[j]);
				swap(matrix_csr.ANC[i], matrix_csr.ANC[j]);
				swap(matrix_csr.AV[i], matrix_csr.AV[j]);
			}
	/*Формируем массив ANL */
	unsigned count = 0;
	vector <unsigned int> ANL = { 0 };
	for (size_t i = 1; i < matrix_csr.nnz; i++)
	{
		count = 0;
		do {
			++count;
		} while (matrix_csr.ANL[i - 1] == matrix_csr.ANL[i] && ++i <= matrix_csr.nnz);
		ANL.push_back(count);
	}

	for (size_t i = 0; i < ANL.size(); i++)
		matrix_csr.ANL[i] = ANL[i];
	/*Сдвигаем массив ANL*/
	for (unsigned i = 1; i < ANL.size(); i++) {
		unsigned t;
		t = matrix_csr.ANL[i - 1] + matrix_csr.ANL[i];
		matrix_csr.ANL[i] = t;
	}

}

void find_psi(csr_matrix  matrix_csr, vector<vector<double>> U, vector<vector<double>> &psi)
{
	vector<vector<double>> tmp(matrix_csr.n, vector<double>(matrix_csr.n, 0.0));
	for (size_t k = 0; k < matrix_csr.n; k++)
	{
		for (size_t i = 0; i < matrix_csr.n; i++)
		{
			double f = 0.0;
			const size_t lb = matrix_csr.ANL[i];
			const size_t ub = matrix_csr.ANL[i + 1];

			for (size_t j = lb; j < ub; j++)
				f += matrix_csr.AV[j] * U[matrix_csr.ANC[j]][k];
			tmp[i][k] += f;
		}
	}
	for (size_t i = 0; i < matrix_csr.n; i++)
		for (size_t j = 0; j < matrix_csr.n; j++)
			psi[i][j] = i == j ? 1 - tmp[i][j] : 0 - tmp[i][j];
}

void find_U(csr_matrix  matrix_csr, vector<vector<double>> &U, vector<vector<double>> psi)
{
	for (size_t i = 0; i < matrix_csr.n; i++)
		for (size_t j = 0; j < matrix_csr.n; j++)
			psi[i][j] = i == j ? 1 + psi[i][j] : psi[i][j];

	double tmp = 0.0;
	for (size_t i = 0; i < matrix_csr.n; i++) {
		for (size_t j = 0; j < matrix_csr.n; j++) {
			tmp = 0.0;
			for (size_t k = 0; k < matrix_csr.n; k++)
				tmp += U[i][k] * psi[k][j];
			U[i][j] = tmp;
		}
	}
}

void print(vector<vector<double>> s) {
	for (vector<double> x : s) {
		for (double y : x)
			cout << setw(8) << setprecision(4) << setiosflags(ios::fixed | ios::showpoint) << y << " ";
		cout << endl;
	}
}

inline double norm(vector<vector<double>> s)
{
	double tmp(0.0);
	for (vector<double> x : s)
		for (double y : x)
			tmp += y * y;
	return sqrt(tmp);
}
