#include <iostream>
#include <fstream>
#include <array>
#include <vector>
#include "math.h"
//#include "Balance.h"

int nomertesta = 1;
const int dim = 2;
const double k = 20; // Коэф. жесткости.
const double m = 0.3;
double F(int i, double t, double* y);
double Function_system(size_t i, double tau, double* y, double* temp);
double** Grid_for_Newton(int N, double tau);
double* Runge_Kutta_4 (double t0, double T, double tau);
double* Evident_Euler (double t0, double T, double tau);
double* Runge_Kutta_2 (double t0, double T, double tau);
double* Not_Evident_Euler (double t0, double T, double tau);


int main(int argc, char* argv[])
{
	for (size_t i = 0; i < 2; ++i)
		std::cout << Evident_Euler(0., 5., 0.01)[i] << std::endl;

	for (size_t i = 0; i < 2; ++i)
		std::cout << Runge_Kutta_2(0., 5., 0.01)[i] << std::endl;

	return 0;
}

double F(int i, double t, double* y) // Systema yravnenii razmera razm, i - nomer yravneniya v systeme
{
	double F;
	if (i == 0) // Pervoe yravnenie systemi
	{
		switch (nomertesta) {
		case 1:
			F = y[1];
			return F;
		case 2:
			F = 2 * y[0] + y[1] * y[1] - 1;
			return F;
		case 3:
			F = 1 - y[0] * y[0] - y[1] * y[1];
			return F;
		case 4:
			F = (0.4 + 0.3 * sin(0.5 * t)) * y[1] + 0.1 * y[0] - y[0] * y[0] * y[0] - y[0] * y[1] * y[1];
			return F;
		case 5:
			F = 10 * (y[1] - y[0]);
			return F;
		case 6:
			F = y[1];
			return F;
		default:
			std::cout << "Invalid number!";
			F = 1;
			return F;
		}
	}
	else if (i == 1) // Vtoroe yravnenie systemi
	{
		switch (nomertesta) {
		case 1:
			F = -k / m * y[0];
			return F;
		case 2:
			F = 6 * y[0] - y[1] * y[1] + 1;
			return F;
		case 3:
			F = 2 * y[0];
			return F;
		case 4:
			F = -(0.4 + 0.3 * sin(0.5 * t)) * y[0] + 0.2 * y[1] - y[1] * y[1] * y[1] - y[1] * y[0] * y[0];
			return F;
		case 5:
			F = y[0] * (28 - y[2]) - y[1];
			return F;
		case 6:
			F = -4 * y[1] - 5 * y[0];
			return F;
		default:
			std::cout << "Nevernii nomer testa!!!!!";
			F = 1;
			return F;
		}
	}
	else if (i == 2)
	{
		F = y[0] * y[1] - 8 / 3 * y[2];
		return F;
		//	cout << "Ne zadana systema pod tekyshii razmer = " << razm << endl;
	}

}

double* Runge_Kutta_4(double t0, double T, double tau)
{
	double* y;
	double* temp;
	double* k1;
	double* k2;
	double* k3;
	double* k4;
	for (size_t i = 0; i < 2; ++i) {
		y = new double[2];
		temp = new double[2];
		k1 = new double[2];
		k2 = new double[2];
		k3 = new double[2];
		k4 = new double[2];
	}
	 
	y[0] = 0.; y[1] = 1.;//начальное условие

	int N = (T + t0) / tau + 1;

	for (size_t i=0; i < N; ++i) {
	
		for(size_t j = 0; i < 2; ++j)
			temp[j] = y[j];

		for (size_t j = 0; i < 2; ++j)
			k1[j] = F(j, i*t0, y);

		for (size_t j = 0; i < 2; ++j)
			temp[j] = y[j] + tau / 2 * k1[j];

		for (size_t j = 0; i < 2; ++j)
			k2[j] = F(i, i*t0 + tau/2, temp);

		for (size_t j = 0; i < 2; ++j)
			temp[j] = y[j] + tau / 2 * k2[j];

		for (size_t j = 0; i < 2; ++j)
			k3[j] = F(i, i * t0 + tau / 2, temp);

		for (size_t j = 0; i < 2; ++j)
			temp[j] = y[j] + tau * k3[j];

		for (size_t j = 0; i < 2; ++j)
			k4[j] = F(i, i*t0 + tau, temp);

		for (size_t j = 0; i < 2; ++j)
			y[j] = y[j] + tau / 6 * (k1[j] + 2 * k2[j] + 2 * k3[j] + k4[j]);
		
	}

	return y;
}

double* Evident_Euler(double t0, double T, double tau)
{
	double* y;
	double* temp;
	for (size_t i = 0; i < 2; ++i) {
		temp = new double[2];
		y = new double[2];
	}

	y[0] = 1.; y[1] = 0.;//начальное условие
	int N = round((t0 + T) / tau) + 1;

	std::ofstream file;
	file.open("1.dat");
	if (file.is_open())
	{
		for (size_t i = 0; i < N; ++i) //i отвечает за шаг
		{
			for (size_t j = 0; j < 2; ++j)
				temp[j] = y[j];

			for (size_t j = 0; j < 2; ++j)
				y[j] = temp[j] + tau * F(j, i * tau, temp);

			file << y[0] << " " << y[1] << std::endl;
		}

	}
	file.close();
	delete[] temp;
	return y;
}

double* Not_Evident_Euler(double t0, double T, double tau)
{
	double* y;
	double* temp;
	for (size_t i = 0; i < 2; ++i) {
		temp = new double[2];
		y = new double[2];
	}

	y[0] = 1.; y[1] = 0.; //начальное условие
	int N = round((t0 + T) / tau) + 1;

	std::ofstream file;
	file.open("2.dat");
	if (file.is_open())
	{

	}
	file.close();
	delete[] temp;
	return y;
}

double* Runge_Kutta_2(double t0, double T, double tau)
{
	double* y;
	double* temp;
	double* k1;
	double* k2;

	for (size_t i = 0; i < 2; ++i) {
		y = new double[2];
		temp = new double[2];
		k1 = new double[2];
		k2 = new double[2];
	}

	y[0] = 1.; y[1] = 0.;//начальное условие

	int N = (T + t0) / tau + 1;

	std::ofstream file;
	file.open("RK_2.dat");

	if (file.is_open()) {
		for (size_t i = 0; i < N; ++i) {

			for (size_t j = 0; j < 2; ++j)
				temp[j] = y[j];

			for (size_t j = 0; j < 2; ++j)
				k1[j] = F(j, i * t0, temp);

			for (size_t j = 0; j < 2; ++j)
				y[j] = temp[j] + tau * k1[j];

			for (size_t j = 0; j < 2; ++j)
				k2[j] = F(j, i * t0 + tau, y);

			for (size_t j = 0; j < 2; ++j)
				y[j] = temp[j] + tau * (0.5 * k1[j] + 0.5 * k2[j]);

			file << y[0] << " " << y[1] << std::endl;
		}
	}
	file.close();

	return y;
}

double Function_system(size_t i, double tau, double* y, double* temp)
{
	double Y;
	Y = y[i] - temp[i] - tau*F(i,tau,y);//задача для маятника, t_n не учитывается функцией
	return Y;
}

double** Grid_for_Newton(int N, double tau)
{
	double** Grid = new double* [2];
		for (size_t i = 0; i < 2; ++i)
			Grid[i] = new double[N];


}