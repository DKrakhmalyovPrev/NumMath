/* 
Пояснение к реализации:
В глобальной области определены граничные условия задачи Дирихле (Они могут изменяться, если угодно), значения в правой части
и другие математические особенности задачи.
В функции PissmanRachford реализован метод Писсмана-Рэчфорда решения матрицы.
Прямое аналитическое решение использует баблиотеку math.h
В ходе работы программы пользователь задаёт два числа -- количество элементов, на которое будет разбита сетка приближения.
Полученные результаты (Аналитического и численного решения) записываются в отдельный файл.
*/
#define _CRT_SECURE_NO_DEPRECATE
#include <new>
#include <cstdlib>
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <algorithm>
# define M_PI 3.14159265358979323846


using namespace std;
FILE *str;

/***************************************************************************************************************/
double u1(double y) { return -(y-0.5)*(y-0.5);}//y*(1-y);} return 5.0; }//y;}// //граничные условия 
double u2(double y) { return -(y-0.5)*(y-0.5);}//y*(1-y);}5.0; }//2+y;}//
double u3(double x) { return fabs(sin(M_PI*x)); }//}5.0; }//x;}//0.75-(x-1)*(x-1);}//
double u4(double x) { return 0.75-(x-1)*(x-1);}//fabs(sin(Pi*x))*exp(x);} 5.0; }//1+x;}//
double f(double x, double y) { return fabs(x - y);}//0; }//-4.0;}//fabs(x-y);}//правая часть
double v(double x, double y) { return 1 - (x - 1)*(x - 1) - (y - 0.5)*(y - 0.5);} //5.0; }//x+y;}//1-(x-1)*(x-1)-(y-0.5)*(y-0.5);} //решение модельной задачи
double L1(double **u, int i, int j, double h1) { return (u[i - 1][j] - 2 * u[i][j] + u[i + 1][j]) / (h1*h1); } //разностные операторы
double L2(double **u, int i, int j, double h2) { return (u[i][j - 1] - 2 * u[i][j] + u[i][j + 1]) / (h2*h2); }
double error(double **u, double **uk, int M, int N) {
	double err = 0.0;
	for (int i = 0;i <= M;i++) {
		for (int j = 0;j<N;j++) {
			err = max(fabs(u[i][j] - uk[i][j]), err);
		}
	}
	return err;
}
/***************************************************************************************************************/
//решение трёхдиагональной СЛАУ методом Писсмана--Рэчфорда
double *PismanRachford(int NN, double *A, double *B, double *C, double *F, double k1, double k2, double w1, double w2)
{
	double *Y; double *Alpha; double *Beta; int i;
	Y = new double[NN + 1];
	Alpha = new double[NN + 1];
	Beta = new double[NN + 1];
	Alpha[1] = k1; Beta[1] = w1;

	for (i = 1; i<NN; i++)
	{
		Alpha[i + 1] = B[i] / (C[i] - A[i] * Alpha[i]);
		Beta[i + 1] = (A[i] * Beta[i] + F[i]) / (C[i] - Alpha[i] * A[i]);
	}

	Y[NN] = (k2*Beta[NN] + w2) / (1 - k2*Alpha[NN]);
	for (i = (NN - 1); i >= 0; i--)
	{
		Y[i] = Alpha[i + 1] * Y[i + 1] + Beta[i + 1];
	}
	delete[] Alpha;
	delete[] Beta;

	return Y;

}
/***************************************************************************************************************/
int main(void)
{

	int M, N;
	double XSize = 1.0;//длина прямоугольника
	double YSize = 1.0;//высота прямоугольника
	double tao, *tmp, *tmp1, *pr, *X, *Y, *A, *B, *C, *A1, *B1, *C1, **F, **FK, **u0, **UK, **U, h1, h2, k1, k2, w1, w2;//UK=Uk+1/2;U=Uk+1;F=Fk+1/2;FK=Fk+1
	int i, j, k;
	/***************************************************************************************************************/
	cout << "Enter M:";
	cin >> M;
	cout << "Enter N:";
	cin >> N;
	str = fopen("res.txt", "w");//инициализируем файл вывода
	X = new double[M + 1]; //инициализируем сетку
	Y = new double[N + 1];
	//cout << "0";
	//инициализируем все остальное   
	A = new double[M + 1];
	B = new double[M + 1];
	C = new double[M + 1];
	A1 = new double[N + 1];
	B1 = new double[N + 1];
	C1 = new double[N + 1];
	pr = new double[M + 1];
	tmp = new double[M + 1];
	F = new double *[M + 1]; for (i = 0;i <= M;i++) { F[i] = new double[N + 1]; }
	U = new double *[M + 1]; for (i = 0;i <= M;i++) { U[i] = new double[N + 1]; }
	FK = new double *[M + 1]; for (i = 0;i <= M;i++) { FK[i] = new double[N + 1]; }
	u0 = new double *[M + 1]; for (i = 0;i <= M;i++) { u0[i] = new double[N + 1]; }
	UK = new double *[M + 1]; for (i = 0;i <= M;i++) { UK[i] = new double[N + 1]; }



	/***************************************************************************************************************/
	tao = 0.001;
	h1 = XSize / M;
	h2 = YSize / N;
	for (i = 0;i <= M;i++) { X[i] = i*h1; } //заполняем сетку координатами
	for (j = 0;j <= N;j++) { Y[j] = j*h2; }
	//cout << "1";

	for (i = 0;i <= M;i++) {
		for (j = 0;j <= N;j++)
		{
			u0[i][j] = /*v(X[i],Y[j])*/0.0;
		}
	}//заполняем  произвольное начальное условие
	//cout << "2";
	/***************************************************************************************************************/
	/*for(i=0;i<=M;i++){ //решение в начальный момент времени
	for (j=0;j<=N;j++){
	U[i][j]=u0[i][j];
	}
	}*/
	U = u0;
	//cout << "3";
	for (i = 0;i <= M;i++) {
		U[i][0] = u3(X[i]);
		U[i][N] = u4(X[i]);
		UK[i][0] = u3(X[i]);
		UK[i][N] = u4(X[i]);
	}
	//cout << "4";
	for (j = 0;j <= N;j++) {
		U[0][j] = u1(Y[j]);
		U[M][j] = u2(Y[j]);
		UK[0][j] = u1(Y[j]);
		UK[M][j] = u2(Y[j]);
	}
	//cout << "5";
	for (i = 1;i <= M;i++) { //заполняем коэффициенты для прогонки по строкам
		A[i] = tao / (2 * (h1*h1));
		C[i] = tao / (h1*h1) + 1;
		B[i] = tao / (2 * (h1*h1));
	}
	//cout << "6";
	for (i = 1;i <= N;i++) { //заполняем коэффициенты для прогонки по столбцам
		A1[i] = tao / (2 * (h2*h2));
		C1[i] = tao / (h2*h2) + 1;
		B1[i] = tao / (2 * (h2*h2));
	}
	//UK=U;
	//do { 


	for (i = 1;i<M;i++) {//заполняем правую часть для первой прогонки
		for (j = 1;j<N;j++) {
			F[i][j] = U[i][j] + 0.5*tao*(L2(U, i, j, h2) + f(X[i], Y[j]));
		}
	}
	for (j = 1;j<N;j++) {//первая прогонка - по строкам
		for (i = 1;i<M;i++) {
			pr[i] = F[i][j];
		}
		tmp = PismanRachford(M, A, B, C, pr, 0, 0, u1(Y[j]), u2(Y[j]));
		for (i = 1;i<M;i++) {
			UK[i][j] = tmp[i];
		}
	}
	for (i = 1;i<M;i++) {//заполняем правую часть для второй прогонки
		for (j = 1;j<N;j++) {
			FK[i][j] = UK[i][j] + 0.5*tao*(L1(UK, i, j, h1) + f(X[i], Y[j]));
		}
	}
	for (i = 1;i<M;i++) {//вторая прогонка - по столбцам
		U[i] = PismanRachford(N, A1, B1, C1, FK[i], 0, 0, u3(X[i]), u4(X[i]));
	}
	//}
	//while (fabs(error(U,UK,M,N))>0.01);

	/***************************************************************************************************************/

	for (j = N;j >= 0;j--) {
		for (i = 0;i <= M;i++) {
			if (U[i][j]<0) fprintf(str, "%.5f   ", U[i][j]);
			else fprintf(str, " %.5f   ", U[i][j]);
		}
		fprintf(str, "\n");
		for (i = 0;i <= M;i++) {
			if (v(X[i], Y[j])<0) fprintf(str, "%.5f   ", v(X[i], Y[j]));
			else fprintf(str, " %.5f   ", v(X[i], Y[j]));
		}
		fprintf(str, "\n\n");
	}
	fprintf(str, "\n\n");
	for (j = N;j >= 0;j--) {
		for (i = 0;i <= M;i++) {
			if (UK[i][j]<0) fprintf(str, "%.5f   ", UK[i][j]);
			else fprintf(str, " %.5f   ", UK[i][j]);
		}
		fprintf(str, "\n");
		for (i = 0;i <= M;i++) {
			if (v(X[i], Y[j])<0) fprintf(str, "%.5f   ", v(X[i], Y[j]));
			else fprintf(str, " %.5f   ", v(X[i], Y[j]));
		}
		fprintf(str, "\n\n");
	}
	return 0;
}
