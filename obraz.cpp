#include "obraz.h"
#include "CImg.h"
#include <cmath>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

using namespace std;


Obraz::Obraz(char *argv) {

	if (cimg_library::cimg::is_file(argv)) {
		cimg_library::CImg<unsigned char> image(argv);
		x = image.width();
		y = image.height();
		r = new int *[y];
		g = new int *[y];
		b = new int *[y];
		for (int i = 0; i < y; i++)
		{
			r[i] = new int[x];
			g[i] = new int[x];
			b[i] = new int[x];
			for (int j = 0; j < x; j++) {
				r[i][j] = image.atXYZC(i, j, 0, 0);
				g[i][j] = image.atXYZC(i, j, 0, 1);
				b[i][j] = image.atXYZC(i, j, 0, 2);
			}
		}
	}
	else {
		cout << "ERROR! WRONG ARGUMENTS - SAVED FILE WILL BE CORRUPTED\n";
		x = 1;
		y = 1;
		r = new int *[y];
		g = new int *[y];
		b = new int *[y];
		for (int i = 0; i < y; i++)
		{
			r[i] = new int[x];
			g[i] = new int[x];
			b[i] = new int[x];
			for (int j = 0; j < x; j++) {
				r[i][j] = 255;
				g[i][j] = 255;
				b[i][j] = 255;
			}
		}


	}
}

void Obraz::save() {
	cimg_library::CImg<unsigned char> image(x, y, 1, 3);
	for (int i = 0; i<y; i++)
	{
		for (int j = 0; j<x; j++) {
			image.set_linear_atXY(r[i][j], i, j, 0, 0);
			image.set_linear_atXY(g[i][j], i, j, 0, 1);
			image.set_linear_atXY(b[i][j], i, j, 0, 2);
		}
	}
	image.save("edited.bmp");
	cout << "Saved as \"edited.bmp\"\n";
}

void Obraz::update(int x1, int y1, int **r1, int **g1, int **b1)
{
	x = x1;
	y = y1;
	r = new int *[y];
	g = new int *[y];
	b = new int *[y];
	for (int i = 0; i < y; i++)
	{
		r[i] = new int[x];
		g[i] = new int[x];
		b[i] = new int[x];
		for (int j = 0; j < x; j++) {
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];
		}
	}

}
// --------------------------------------------TASK 1 ------------------------------------------------------
void Obraz::brightness(int c)
{
	for (int i = 0; i < y; i++)
	{
		for (int j = 0; j < x; j++)
		{
			int ur = r[i][j] + c;
			int ug = g[i][j] + c;
			int ub = b[i][j] + c;
			if (ur > 255) r[i][j] = 255;
			else if (ur < 0) r[i][j] = 0;
			else r[i][j] = ur;
			if (ug > 255) g[i][j] = 255;
			else if (ug < 0) g[i][j] = 0;
			else g[i][j] = ug;
			if (ub > 255) b[i][j] = 255;
			else if (ub < 0) b[i][j] = 0;
			else b[i][j] = ub;
		}
	}

}


void Obraz::contrast(float a)
{
	float gg = r[23][23] / 31;
	cout << (float)gg;
	cout << "\n";
	cout << (long double)pow(gg, 1);
	for (int i = 0; i < y; i++)
	{
		for (int j = 0; j < x; j++)
		{

			double ur = pow(r[i][j] / 255, a) * 255;
			double ug = pow(g[i][j] / 255, a) * 255;
			double ub = pow(b[i][j] / 255, a) * 255;

			r[i][j] = ur;
			g[i][j] = ug;
			b[i][j] = ub;
		}
	}

}

void Obraz::negative()
{
	for (int i = 0; i < y; i++)
	{
		for (int j = 0; j < x; j++)
		{
			r[i][j] = 255 - r[i][j];
			g[i][j] = 255 - g[i][j];
			b[i][j] = 255 - b[i][j];
		}
	}
}

void Obraz::h_flip()
{
	int k = x;
	for (int i = 0; i < k; i++)
	{
		k--;
		for (int j = 0; j < y; j++)
		{
			int r1 = r[i][j];
			r[i][j] = r[k][j];
			r[k][j] = r1;
			int g1 = g[i][j];
			g[i][j] = g[k][j];
			g[k][j] = g1;
			int b1 = b[i][j];
			b[i][j] = b[k][j];
			b[k][j] = b1;
		}
	}
}

void Obraz::v_flip()
{
	int k = y;
	for (int j = 0; j < k; j++)
	{
		k--;
		for (int i = 0; i < x; i++)
		{
			int r1 = r[i][j];
			r[i][j] = r[i][k];
			r[i][k] = r1;
			int g1 = g[i][j];
			g[i][j] = g[i][k];
			g[i][k] = g1;
			int b1 = b[i][j];
			b[i][j] = b[i][k];
			b[i][k] = b1;
		}
	}
}

void Obraz::d_flip()
{
	this->v_flip();
	this->h_flip();
}

void Obraz::enlarge(int c)
{
	int x1 = x*c;
	int y1 = y*c;
	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y1];
	g1 = new int *[y1];
	b1 = new int *[y1];
	for (int i = 0; i < y1; i++)
	{
		r1[i] = new int[x1];
		g1[i] = new int[x1];
		b1[i] = new int[x1];
	}

	int ii = 0;
	for (int i = 0; i < y; i++)
	{
		int jj = 0;

		for (int j = 0; j < x; j++) {
			int ii2 = ii;
			for (int k = 0; k < c; k++)
			{
				int jj2 = jj;
				for (int l = 0; l < c; l++)
				{
					r1[ii][jj] = r[i][j];
					g1[ii][jj] = g[i][j];
					b1[ii][jj] = b[i][j];
					jj++;
					if (l == c - 1) jj = jj2;
				}
				ii++;
				if (k == c - 1) ii = ii2;
			}
			jj = jj + c;
		}
		ii = ii + c;
	}

	update(x1, y1, r1, g1, b1);
}

void Obraz::shrink(int c)
{
	int x1 = x / c;
	int y1 = y / c;
	int **r1;
	int **g1;
	int **b1;


	r1 = new int *[y1];
	g1 = new int *[y1];
	b1 = new int *[y1];
	int ii = 0;
	for (int i = 0; i < y1; i++)
	{
		r1[i] = new int[x1];
		g1[i] = new int[x1];
		b1[i] = new int[x1];
		int jj = 0;
		for (int j = 0; j < x1; j++) {
			r1[i][j] = r[ii][jj];
			g1[i][j] = g[ii][jj];
			b1[i][j] = b[ii][jj];
			jj = jj + c;
		}
		ii = ii + c;
	}
	update(x1, y1, r1, g1, b1);
}


int Obraz::mid(int i, int j, int **r)
{
	int arr[9];
	arr[0] = r[i - 1][j - 1];
	arr[1] = r[i - 1][j];
	arr[2] = r[i - 1][j + 1];
	arr[3] = r[i][j - 1];
	arr[4] = r[i][j];
	arr[5] = r[i][j + 1];
	arr[6] = r[i + 1][j - 1];
	arr[7] = r[i + 1][j];
	arr[8] = r[i + 1][j + 1];

	for (int k = 0; k < 9; k++) {
		for (int l = k + 1; l < 9; l++)
		{

			if (arr[k] > arr[l]) {
				int u = arr[k];
				arr[k] = arr[l];
				arr[l] = u;
			}
		}
	}
	return arr[4];
}

void Obraz::median()
{
	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];
	for (int i = 1; i < y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 1; j < x - 1; j++) {
			r1[i][j] = mid(i, j, r);
			g1[i][j] = mid(i, j, g);
			b1[i][j] = mid(i, j, b);
		}
	}

	for (int i = 1; i < y - 1; i++) {
		for (int j = 1; j < x - 1; j++) {
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];

		}
	}
}

int Obraz::harmonic(int i, int j, int **r)
{
	float k = 0;
	float arr[9];
	arr[0] = r[i - 1][j - 1];
	arr[1] = r[i - 1][j];
	arr[2] = r[i - 1][j + 1];
	arr[3] = r[i][j - 1];
	arr[4] = r[i][j];
	arr[5] = r[i][j + 1];
	arr[6] = r[i + 1][j - 1];
	arr[7] = r[i + 1][j];
	arr[8] = r[i + 1][j + 1];

	for (int i = 0; i < 9; i++) {
		arr[i] = 1 / arr[i];
		k += arr[i];
	}
	int l = 9 / k;
	return l;
}


void Obraz::hmean() {
	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];
	for (int i = 1; i < y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 1; j < x - 1; j++) {
			//zmiana
			r1[i][j] = harmonic(i, j, r);
			g1[i][j] = harmonic(i, j, g);
			b1[i][j] = harmonic(i, j, b);
		}
	}

	for (int i = 1; i < y - 1; i++) {
		for (int j = 1; j < x - 1; j++) {
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];
		}
	}

}

double Obraz::maxPixel()
{
	double tempMaximum = 0;
	for (int i = 0; i < x; i++)
		for (int j = 0; j < y; j++)
		{
			if ((r[i][j] + g[i][j] + b[i][j]) >= tempMaximum)
			{
				tempMaximum = r[i][j] + g[i][j] + b[i][j];
			}
		}
	return tempMaximum / 3;
}

double Obraz::square(int number)
{
	return number*number;
}

double Obraz::meanSquareError(Obraz image2)
{
	double numberOfElements = x*y;
	double sum = 0;
	if (x == image2.x&&y == image2.y)
	{
		for (int i = 0; i < x; i++)
			for (int j = 0; j < y; j++)
			{
				sum += square(b[i][j] - image2.b[i][j] + g[i][j] - image2.g[i][j] + r[i][j] - image2.r[i][j]);

			}
	}
	else { cout << "Different image size" << endl; return 0; }

	return sum / numberOfElements;
}

double Obraz::peakMeanSquareError(Obraz image2)
{
	double pixel = maxPixel();

	return	meanSquareError(image2) / square(pixel);
}


double Obraz::signalToNoiseRatio(Obraz image2)
{
	double numberOfElements = x*y;
	double denominator = meanSquareError(image2) * numberOfElements;
	double sumOfIntensities = 0;
	for (int i = 0; i < x; i++)
		for (int j = 0; j < y; j++)
		{
			sumOfIntensities += square(r[i][j] + g[i][j] + b[i][j]);
		}
	return 10 * log10(sumOfIntensities / denominator);
}

double Obraz::peakSignalToNoiseRatio(Obraz image2)
{
	double numberOfElements = x*y;
	return 10 * log10(maxPixel() / (meanSquareError(image2)*numberOfElements));
}

double Obraz::maxDiff(Obraz image2)
{
	double tempMaximum = 0;
	for (int i = 0; i < x; i++)
		for (int j = 0; j < y; j++)
		{
			if ((abs(r[i][j] - image2.r[i][j] + g[i][j] - image2.g[i][j] + b[i][j] - image2.b[i][j])) > tempMaximum)
			{
				tempMaximum = abs(r[i][j] - image2.r[i][j] + g[i][j] - image2.g[i][j] + b[i][j]) - image2.b[i][j];
			}
		}
	return tempMaximum;
}

//---------------------------------------------TASK 2-------------------------------------------------------
void Obraz::histogram(char c) {
	switch (c) {
	case 'r':
	{
		for (int k = 0; k<256; k++) {
			numbers[k] = 0;
			for (int i = 0; i<y; i++) {
				for (int j = 0; j<x; j++) {
					if (k == r[i][j])numbers[k]++;
				}
			}
		}

		int maxVal = 1;
		for (int i = 0; i<256; i++) {
			if (numbers[i]>maxVal)maxVal = numbers[i];
		}
		cimg_library::CImg<unsigned char> image(256, maxVal, 1, 3);
		for (int i = 0; i<256; i++) {
			int tempMax = maxVal - numbers[i] - 1;
			for (int j = maxVal - 1; j >= tempMax; j--) {
				image.set_linear_atXY(255, i, j, 0, 0);
				image.set_linear_atXY(255, i, j, 0, 1);
				image.set_linear_atXY(255, i, j, 0, 2);
			}
		}
		image.save("histogram.bmp");
		cout << "Saved as \"histogram.bmp\"\n";
		break;
	}
	case 'g':
	{
		for (int k = 0; k<256; k++) {
			numbers[k] = 0;
			for (int i = 0; i<y; i++) {
				for (int j = 0; j<x; j++) {
					if (k == g[i][j])numbers[k]++;
				}
			}
		}

		int maxVal = 1;
		for (int i = 0; i<256; i++) {
			if (numbers[i]>maxVal)maxVal = numbers[i];
		}
		cimg_library::CImg<unsigned char> image(256, maxVal, 1, 3);
		for (int i = 0; i<256; i++) {
			int tempMax = maxVal - numbers[i] - 1;
			for (int j = maxVal - 1; j >= tempMax; j--) {
				image.set_linear_atXY(255, i, j, 0, 0);
				image.set_linear_atXY(255, i, j, 0, 1);
				image.set_linear_atXY(255, i, j, 0, 2);
			}
		}
		image.save("histogram.bmp");
		cout << "Saved as \"histogram.bmp\"\n";
		break;
	}
	case 'b':
	{
		for (int k = 0; k<256; k++) {
			numbers[k] = 0;
			for (int i = 0; i<y; i++) {
				for (int j = 0; j<x; j++) {
					if (k == b[i][j])numbers[k]++;
				}
			}
		}

		int maxVal = 1;
		for (int i = 0; i<256; i++) {
			if (numbers[i]>maxVal)maxVal = numbers[i];
		}
		cimg_library::CImg<unsigned char> image(256, maxVal, 1, 3);
		for (int i = 0; i<256; i++) {
			int tempMax = maxVal - numbers[i] - 1;
			for (int j = maxVal - 1; j >= tempMax; j--) {
				image.set_linear_atXY(255, i, j, 0, 0);
				image.set_linear_atXY(255, i, j, 0, 1);
				image.set_linear_atXY(255, i, j, 0, 2);
			}
		}
		image.save("histogram.bmp");
		cout << "Saved as \"histogram.bmp\"\n";
		break;
	}
	default:
	{
		cout << "ERROR: wrong argument!\n";
		break;
	}
	}
}


void Obraz::hraleigh(int gmin, float alpha) {
	this->histogram('r');
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			int sumHistogram = 0;
			for (int m = 255; m != r[i][j]; m--)sumHistogram += numbers[m];
			float N = x*y;
			N *= 1.0;
			float result = N / sumHistogram;
			result = log(result);
			result = result *alpha*alpha * 2;
			result = sqrt(result);
			result += gmin;
			r[i][j] = result;
		}
	}

	this->histogram('g');
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			int sumHistogram = 0;
			for (int m = 255; m != g[i][j]; m--)sumHistogram += numbers[m];

			float N = x*y;
			N *= 1.0;
			float result = N / sumHistogram;
			result = log(result);
			result = result *alpha*alpha * 2;
			result = sqrt(result);
			result += gmin;
			g[i][j] = result;
		}
	}

	this->histogram('b');
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			int sumHistogram = 0;
			for (int m = 255; m != b[i][j]; m--)sumHistogram += numbers[m];
			float N = x*y;
			N *= 1.0;
			float result = N / sumHistogram;
			result = log(result);
			result = result *alpha*alpha * 2;
			result = sqrt(result);
			result += gmin;
			b[i][j] = result;
		}
	}
}

int maskSum(int mask[9], int i, int j, int **r) {
	int arr[9];
	arr[0] = r[i - 1][j - 1] * mask[0];
	arr[1] = r[i - 1][j] * mask[1];
	arr[2] = r[i - 1][j + 1] * mask[2];
	arr[3] = r[i][j - 1] * mask[3];
	arr[4] = r[i][j] * mask[4];
	arr[5] = r[i][j + 1] * mask[5];
	arr[6] = r[i + 1][j - 1] * mask[6];
	arr[7] = r[i + 1][j] * mask[7];
	arr[8] = r[i + 1][j + 1] * mask[8];
	int sum = 0;
	for (int k = 0; k<9; k++) sum += arr[k];
	return sum;
}

void Obraz::convolution(int mask[9]) {

	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];
	for (int i = 1; i<y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 1; j<x - 1; j++) {
			r1[i][j] = maskSum(mask, i, j, r);
			g1[i][j] = maskSum(mask, i, j, g);
			b1[i][j] = maskSum(mask, i, j, b);
		}
	}

	for (int i = 1; i<y - 1; i++) {
		for (int j = 1; j<x - 1; j++) {
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];

		}
	}this->save();
}

int orobertsiCalculate(int i, int j, int **r) {
	float result = r[i + 1][j + 1] / 1.0;
	result = r[i][j] - result;
	result *= result;
	float result2 = r[i][j + 1] - r[i + 1][j];
	result2 /= 1.0;
	result2 *= result2;
	result += result2;
	result = sqrt(result);
	return result;
}

void Obraz::orobertsi() {

	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];

	for (int i = 0; i<y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 0; j<x - 1; j++) {
			r1[i][j] = orobertsiCalculate(i, j, r);
			g1[i][j] = orobertsiCalculate(i, j, g);
			b1[i][j] = orobertsiCalculate(i, j, b);
		}
	}

	for (int i = 0; i<y - 1; i++) {
		for (int j = 0; j<x - 1; j++) {
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];

		}
	}this->save();

}

int Obraz::kirschOpCalc(int i, int j, int **r) {

	int arr[15];
	arr[0] = arr[8] = r[i - 1][j - 1];
	arr[7] = r[i - 1][j];
	arr[6] = arr[14] = r[i - 1][j + 1];
	arr[1] = arr[9] = r[i][j - 1];
	arr[5] = arr[13] = r[i][j + 1];
	arr[2] = arr[10] = r[i + 1][j - 1];
	arr[3] = arr[11] = r[i + 1][j];
	arr[4] = arr[12] = r[i + 1][j + 1];

	float Sresult = 0;
	float Tresult = 0;
	float result = 0;
	int tempMax = 0;

	for (int i = 0; i < 8; i++) {
		Sresult = arr[i] + arr[i + 1] + arr[i + 2];
		Tresult = arr[i + 3] + arr[i + 4] + arr[i + 5] + arr[i + 6] + arr[i + 7];

		result = abs((5 * Sresult) - (3 * Tresult));
		if (result > tempMax) tempMax = result;
	} return tempMax;
}


void Obraz::kirsh() {

	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];

	for (int i = 1; i < y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 1; j < x - 1; j++) {
			//zmiana
			r1[i][j] = max(1, kirschOpCalc(i, j, r));
			g1[i][j] = max(1, kirschOpCalc(i, j, g));
			b1[i][j] = max(1, kirschOpCalc(i, j, b));
		}
	}

	for (int i = 1; i < y - 1; i++) {
		for (int j = 1; j < x - 1; j++) {
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];
		}
	}
	this->save();
}

int Obraz::uolisValue(int i, int j, int **r)
{
	int k = 1;
	double arr[5];
	arr[0] = r[i - 1][j];
	arr[1] = r[i][j - 1];
	arr[2] = r[i][j + 1];
	arr[3] = r[i + 1][j];
	arr[4] = r[i][j];

	for (int p = 0; p < 4; p++) {
		if (arr[p] == 0) arr[p] = 1;
		k *= arr[p];
	}

	double result = 0;
	result = (1 / 4)*log10((pow(arr[4], 4)) / k) * 10;
	if (result > log10(255)) return 255;
	if (result < (-log10(255))) return 255;
	else return result;
}
void Obraz::uolis() {

	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];

	for (int i = 1; i<y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 1; j < x - 1; j++) {
			r1[i][j] = uolisValue(i, j, r);
			g1[i][j] = uolisValue(i, j, g);
			b1[i][j] = uolisValue(i, j, b);
		}
	}
	for (int i = 1; i < y - 1; i++) {
		for (int j = 1; j < x - 1; j++) {
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];
		}
	}this->save();

}

float Obraz::mean(char c)
{
	//tablica z wartociami histogramu
	float N = x*y;
	double sum = 0;
	this->makeHistogram(c);

	for (int hist = 0; hist < 256; hist++)
	{
		sum += hist * numbers[hist];
	}
	return ((1 / N) * sum);
}

float Obraz::variance(char c) {

	float sum = 0;
	float N = x*y;
	float vmean = mean(c);
	this->makeHistogram(c);

	for (int hist = 0; hist < 256; hist++)
	{
		sum += (hist - vmean) * (hist - vmean) * numbers[hist];
	}
	return (1 / (N * sum));
}

float Obraz::standardDeviation(char c) {

	return sqrt(variance(c));
}
float Obraz::variationCoefficient1(char c) {

	return standardDeviation(c) / mean(c);
}
float Obraz::asymmetryCoefficient(char c) {

	float sum = 0;
	float N = x*y;
	float amean = mean(c);
	int histogram[256];
	this->makeHistogram(c);
	for (int hist = 0; hist< 256; hist++)
	{
		sum += pow((hist - amean), 3) * numbers[hist];
	}

	return 1 / pow(standardDeviation(c), 3) * (1 / (N * sum));
}
float Obraz::flatteningCoefficient(char c) {
	float sum = 0;
	float N = x*y;
	float fmean = mean(c);
	this->makeHistogram(c);

	for (int hist = 0; hist < 256; hist++)
	{
		sum += pow((hist - fmean), 4) * numbers[hist];
	}

	return 1 / pow(standardDeviation(c), 4) * (1 / (N * sum)) - 3;

}
float Obraz::variationCoefficient2(char c) {

	float sum = 0;
	float N = x*y;
	this->makeHistogram(c);
	for (int hist = 0; hist < 256; hist++)
	{
		sum += pow(numbers[hist], 2);
	}

	return pow(1 / N, 2) * sum;
}
float Obraz::informationSourceEntropy(char c) {
	float sum = 0;
	float N = x*y;
	this->makeHistogram(c);
	for (int hist = 0; hist < 256; hist++)
	{
		if (numbers[hist] != 0)			//calculates only when hist !=0 so log2 won't be equal to -infinity
		{
			sum += numbers[hist] * log2(numbers[hist] / N);
		}
	}
	return (-1 / N) * sum;

}

void Obraz::makeHistogram(char c) {

	switch (c)
	{
	case 'r':
		for (int i = 0; i < 256; i++) numbers[i] = 0;
		for (int x1 = 0; x1 < x; x1++)
			for (int y1 = 0; y1 < y; y1++)
				numbers[r[x1][y1]]++;
		break;
	case 'g':
		for (int i = 0; i < 256; i++) numbers[i] = 0;
		for (int x1 = 0; x1 < x; x1++)
			for (int y1 = 0; y1 < y; y1++)
				numbers[g[x1][y1]]++;
		break;
	case 'b':
		for (int i = 0; i < 256; i++) numbers[i] = 0;
		for (int x1 = 0; x1 < x; x1++)
			for (int y1 = 0; y1 < y; y1++)
				numbers[b[x1][y1]]++;
		break;
	default:
		cout << "Wrong value. Must choose r/b/g";
		break;
	}
}

//===================================TASK 3================================================

void Obraz::Choose_Structural_Element(int structural_element, bool str_el[3][3])
{
	for (int x = 0; x < 3; x++)
	{
		for (int y = 0; y < 3; y++)
		{
			str_el[x][y] = 0;
		}
	}

	if (structural_element == 1)
	{
		str_el[1][1] = str_el[2][1] = 1;
	}
	else if (structural_element == 2)
	{
		str_el[1][1] = str_el[1][2] = 1;
	}
	else if (structural_element == 3)
	{
		for (int x = 0; x < 3; x++)
		{
			for (int y = 0; y < 3; y++)
			{
				str_el[x][y] = 1;
			}
		}
	}
	else if (structural_element == 4)
	{
		str_el[0][1] = str_el[1][0] = 1;
		str_el[1][1] = 1;
		str_el[1][2] = str_el[2][1] = 1;
	}
	else if (structural_element == 5)
	{
		str_el[1][1] = str_el[1][2] = str_el[2][1] = 1;
	}
	else if (structural_element == 6)
	{
		str_el[2][1] = str_el[1][2] = 1;
	}
	else if (structural_element == 7)
	{
		str_el[0][1] = str_el[1][1] = str_el[2][1] = 1;
	}
	else if (structural_element == 8)
	{
		str_el[0][1] = str_el[2][1] = 1;
	}
	else if (structural_element == 9)
	{
		str_el[1][1] = str_el[0][1] = str_el[0][2] = 1;
	}
	else if (structural_element == 10)
	{
		str_el[1][1] = str_el[1][0] = str_el[2][0] = 1;
	}
}
void Obraz::Choose_HMT_Structural_Element(int structural_element, int type, bool str_el[3][3], bool iistr_el[3][3]) //niedokoñczone, nieokrelone wartoci dla bialych elementów
{
	for (int x = 0; x < 3; x++)
	{
		for (int y = 0; y < 3; y++)
		{
			str_el[x][y] = 0;
		}
	}

	if (structural_element == 11)
	{
		if (type == 1)
		{
			str_el[0][0] = str_el[0][1] = str_el[0][2] = 1;
		}
		else if (type == 2)
		{
			str_el[0][0] = str_el[1][0] = str_el[2][0] = 1;

		}
		else if (type == 3)
		{
			str_el[2][0] = str_el[2][1] = str_el[2][2] = 1;

		}
		else if (type == 4)
		{
			str_el[0][2] = str_el[1][2] = str_el[2][2] = 1;

		}
	}
	else if (structural_element == 12)
	{
		if (type == 1)
		{
			str_el[0][2] = str_el[1][1] = str_el[1][2] = str_el[2][2] = 1;

		}
		else if (type == 2)
		{
			str_el[0][1] = str_el[0][2] = str_el[1][1] = str_el[1][2] = 1;

		}
		else if (type == 3)
		{
			str_el[0][0] = str_el[0][1] = str_el[0][2] = str_el[1][1] = 1;

		}
		else if (type == 4)
		{
			str_el[0][0] = str_el[0][1] = str_el[1][0] = str_el[1][1] = 1;

		}
		else if (type == 5)
		{
			str_el[0][0] = str_el[1][0] = str_el[1][1] = str_el[2][0] = 1;

		}
		else if (type == 6)
		{
			str_el[1][0] = str_el[1][1] = str_el[2][0] = str_el[2][1] = 1;

		}
		else if (type == 7)
		{
			str_el[1][1] = str_el[2][0] = str_el[2][1] = str_el[2][2] = 1;

		}
		else if (type == 8)
		{
			str_el[1][1] = str_el[1][2] = str_el[2][1] = str_el[2][2] = 1;

		}
	}
	//=====================================DRUGA MASKA=============================================
	for (int x = 0; x < 3; x++)
	{
		for (int y = 0; y < 3; y++)
		{
			iistr_el[x][y] = 0;
		}
	}
	if (structural_element == 11)
	{


		if (type == 1 || type == 2 || type == 3 || type == 4)
			iistr_el[1][1] = 1;
	}
	else if (structural_element == 12)
	{
		if (type == 1)
		{
			iistr_el[0][0] = iistr_el[1][0] = iistr_el[2][0] = 1;

		}
		else if (type == 2)
		{
			iistr_el[1][0] = iistr_el[2][0] = iistr_el[2][1] = 1;

		}
		else if (type == 3)
		{
			iistr_el[2][0] = iistr_el[2][1] = iistr_el[2][2] = 1;

		}
		else if (type == 4)
		{
			iistr_el[1][2] = iistr_el[2][1] = iistr_el[2][2] = 1;

		}
		else if (type == 5)
		{
			iistr_el[0][2] = iistr_el[1][2] = iistr_el[2][2] = 1;

		}
		else if (type == 6)
		{
			iistr_el[0][1] = iistr_el[0][2] = iistr_el[1][2] = 1;

		}
		else if (type == 7)
		{
			iistr_el[0][0] = iistr_el[0][1] = iistr_el[0][2] = 1;

		}
		else if (type == 8)
		{
			iistr_el[0][0] = iistr_el[0][1] = iistr_el[1][0] = 1;
		}
	}
}


bool Obraz::isWhite(int i, int j)
{
	if (r[i][j] == 255 && g[i][j] == 255 && b[i][j] == 255) return 1;
	if (r[i][j] == 0 && g[i][j] == 0 && b[i][j] == 0) return 0;
}

void Obraz::dilation(bool str_el[3][3])
{
	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];

	bool stan;

	for (int i = 1; i < y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];


		for (int j = 1; j < x - 1; j++) {
			stan = false;//false -->niezamalowywaæ
			if (isWhite(i, j)) {
				for (int p = 0; p < 3; p++) {
					for (int d = 0; d < 3; d++) {
						if (str_el[p][d] == 1 && (!isWhite(i - 1 + p, j - 1 + d)))
						{
							stan = true; //true-->zamalowaæ
						}

					}
				}
				if (stan)
				{
					r1[i][j] = 0;
					b1[i][j] = 0;
					g1[i][j] = 0;
				}
				if (!stan)
				{
					r1[i][j] = 255;
					b1[i][j] = 255;
					g1[i][j] = 255;
				}
			}
			else
			{
				r1[i][j] = 0;
				b1[i][j] = 0;
				g1[i][j] = 0;
			}
		}
	}

	for (int i = 1; i < y - 1; i++)
	{
		for (int j = 1; j < x - 1; j++)
		{
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];
		}

	} this->save();

}

void Obraz::erosion(bool str_el[3][3])
{
	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];

	bool stan;

	for (int i = 1; i < y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];


		for (int j = 1; j < x - 1; j++) {
			stan = true;//true -->zamalowane
			if (!isWhite(i, j)) {
				for (int p = 0; p < 3; p++) {
					for (int d = 0; d < 3; d++) {
						if (str_el[p][d] == 1 && (isWhite(i - 1 + p, j - 1 + d)))// (r[i - 1 + p][j - 1 + d] == 0 && g[i - 1 + p][j - 1 + d] == 0 && b[i - 1 + p][j - 1 + d]) == 0)
						{
							stan = false; //true-->zamalowaæ
						}

					}
				}
				if (stan)
				{
					r1[i][j] = 0;
					b1[i][j] = 0;
					g1[i][j] = 0;
				}
				if (!stan)
				{
					r1[i][j] = 255;
					b1[i][j] = 255;
					g1[i][j] = 255;
				}
			}
			else
			{
				r1[i][j] = 255;
				b1[i][j] = 255;
				g1[i][j] = 255;
			}
		}
	}

	for (int i = 1; i < y - 1; i++)
	{
		for (int j = 1; j < x - 1; j++)
		{
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];
		}

	} this->save();

}

/*void Obraz::HMTtransformation(int structural_element, int z)// bool str_el[3][3], bool iistr_el[3][3])
{
int **r1;
int **g1;
int **b1;

r1 = new int *[y];
g1 = new int *[y];
b1 = new int *[y];

bool stan;
bool tab1[3][3];
bool tab2[3][3];
//int typ = 0;

for (int i = 1; i < y - 1; i++)
{
r1[i] = new int[x];
g1[i] = new int[x];
b1[i] = new int[x];

for (int j = 1; j < x - 1; j++) {
stan = true;//true-->zamalowaæ
if (!isWhite(i, j)) { // czy zadziała tylko na str_el=12?
/*
if (structural_element == 11) typ= 5; //zrobi 4 razy
if (structural_element == 12) typ = 9; //zrobi 8 razy

for (int g = 1; g < typ; g++) { *
Choose_HMT_Structural_Element(structural_element, z, tab1, tab2);
for (int p = 0; p < 3; p++) {
for (int d = 0; d < 3; d++) {
if (!((tab2[p][d] == 1 && (isWhite(i - 1 + p, j - 1 + d))) || (tab1[p][d] == 1 && (!isWhite(i - 1 + p, j - 1 + d)))))//
{
stan = false; //false -->niezamalowywaæ
}
} // dać break?
}
//	}

if (stan)
{
cout << "O";
r1[i][j] = 0;// 255;
b1[i][j] = 0;// 255;
g1[i][j] = 0;// 255;
}

if (!stan)
{
cout << "X";
r1[i][j] = 255;// 0;
b1[i][j] = 255;//  0;
g1[i][j] = 255;// 0;
}
}

else
{

r1[i][j] = 255;// 0;
b1[i][j] = 255;// 0;
g1[i][j] = 255;//0;
}
}
}


for (int i = 1; i < y - 1; i++)
{
for (int j = 1; j < x - 1; j++)
{
r[i][j] = r1[i][j];
g[i][j] = g1[i][j];
b[i][j] = b1[i][j];
}
} this->save();
}
*/

void Obraz::HMTtransformation(int structural_element, int z) {
	// Structural element
	//const int 
	z = 1; // TODO: Make it into an argument
	bool hit[3][3];  // Matches 1
	bool miss[3][3]; // Matches 0
	Choose_HMT_Structural_Element(structural_element, z, hit, miss);

	// Output image
	int** r1 = new int*[y];
	int** g1 = new int*[y];
	int** b1 = new int*[y];

	// Visit each pixel
	bool white, match;
	for (int i = 1; i < y - 1; i++) {
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 1; j < x - 1; j++) {
			// Check if neighbours match st. element
			match = true;
			for (int p = 0; p < 3; p++) {
				for (int d = 0; d < 3; d++) {
					white = isWhite(i - 1 + p, j - 1 + d);
					if ((hit[p][d] && !white) || (miss[p][d] && white)) {
						//if (!(hit[p][d] == 1 && white) || !(miss[p][d] == 1 && !white)) {
						match = false;
					}
				}
			}
			// Color the pixel accordingly
			if (match) {
				r1[i][j] = 255;
				g1[i][j] = 255;
				b1[i][j] = 255;
			}
			else {
				r1[i][j] = 0;
				g1[i][j] = 0;
				b1[i][j] = 0;
			}
		}
	}

	// Save output
	for (int i = 1; i < y - 1; i++)
	{
		for (int j = 1; j < x - 1; j++)
		{
			r[i][j] = r1[i][j];
			g[i][j] = g1[i][j];
			b[i][j] = b1[i][j];
		}
	} this->save();
}


bool isHit(int **rR, int i, int j) {
	if (rR[i][j] == 0) return true;
	if (rR[i - 1][j] == 0) return true;
	if (rR[i][j - 1] == 0) return true;
	if (rR[i + 1][j] == 0) return true;
	if (rR[i][j + 1] == 0) return true;
	else return false;
}

void Obraz::mydilation(int **rR, int **gG, int **bB)
{
	int **r1;
	int **g1;
	int **b1;

	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];


	for (int i = 1; i < y - 1; i++)
	{
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 1; j < x - 1; j++) {
			if (isHit(rR, i, j) && isHit(gG, i, j) && isHit(bB, i, j))
			{
				r1[i][j] = 0; g1[i][j] = 0; b1[i][j] = 0;
			}
			else { r1[i][j] = 255; g1[i][j] = 255; b1[i][j] = 255; }
		}
	}


	for (int i = 1; i<y - 1; i++) {
		for (int j = 1; j<x - 1; j++) {
			rR[i][j] = r1[i][j];
			gG[i][j] = g1[i][j];
			bB[i][j] = b1[i][j];
		}
	}

}

bool Obraz::isMatrixEqual(int **r1, int **g1, int **b1, int **r2, int **g2, int **b2) {
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			if (r1[i][j] != r2[i][j] || g1[i][j] != g2[i][j] || b1[i][j] != b2[i][j]) return false;
		}
	}
	return true;
}

void Obraz::m2(int ii, int jj) {
	//Tworzymy biały obraz
	if (r[ii][jj] == 255) {
		int **rx1;
		int **gx1;
		int **bx1;

		rx1 = new int *[y];
		gx1 = new int *[y];
		bx1 = new int *[y];
		for (int i = 0; i < y; i++)
		{
			rx1[i] = new int[x];
			gx1[i] = new int[x];
			bx1[i] = new int[x];
			for (int j = 0; j < x; j++) {
				rx1[i][j] = 255;
				gx1[i][j] = 255;
				bx1[i][j] = 255;
			}
		}

		rx1[ii][jj] = 0;
		gx1[ii][jj] = 0;
		bx1[ii][jj] = 0;
		//obraz do porownywania
		int **rx0;
		int **gx0;
		int **bx0;

		rx0 = new int *[y];
		gx0 = new int *[y];
		bx0 = new int *[y];
		for (int i = 0; i < y; i++)
		{
			rx0[i] = new int[x];
			gx0[i] = new int[x];
			bx0[i] = new int[x];
			for (int j = 0; j < x; j++) {
				rx0[i][j] = rx1[i][j];
				gx0[i][j] = gx1[i][j];
				bx0[i][j] = bx1[i][j];
			}
		}

		//inverse obraz
		int **rA;
		int **gA;
		int **bA;

		rA = new int *[y];
		gA = new int *[y];
		bA = new int *[y];
		for (int i = 0; i < y; i++)
		{
			rA[i] = new int[x];
			gA[i] = new int[x];
			bA[i] = new int[x];
			for (int j = 0; j < x; j++) {
				if (r[i][j] == 255) { rA[i][j] = 0; gA[i][j] = 0; bA[i][j] = 0; }
				else { rA[i][j] = 255; gA[i][j] = 255; bA[i][j] = 255; }
			}
		}


		//for(int h=0;h<200;h++){
		do {

			for (int i = 0; i < y; i++)
			{
				for (int j = 0; j < x; j++) {
					rx0[i][j] = rx1[i][j];
					gx0[i][j] = gx1[i][j];
					bx0[i][j] = bx1[i][j];
				}
			}


			mydilation(rx1, gx1, bx1);

			for (int i = 0; i < y; i++)
			{
				for (int j = 0; j < x; j++) {
					if (r[i][j] == 0 && rx1[i][j] == 0) { rx1[i][j] == 255; gx1[i][j] = 255; bx1[i][j] = 255; }
				}
			}

			if (isMatrixEqual(rx1, gx1, bx1, rx0, gx0, bx0)) break;
		} while (1);


		for (int i = 0; i<y; i++) {
			for (int j = 0; j<x; j++) {
				if (rx1[i][j] == 0 && gx1[i][j] == 0 && bx1[i][j] == 0) { r[i][j] = 0; g[i][j] = 0; b[i][j] = 0; }
			}
		}
		this->save();
	}
	else cout << "Given pixel is not a part of a hole - its black \n";
}

bool Obraz::areSeedsIn(int **r1, int seedPoint) {
	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			if (r[i][j] == seedPoint && r1[i][j] == -1) return false;
		}
	}
	return true;
}

void Obraz::region(int seedPoint, int backgroundColor) {
	int **r1;
	int **g1;
	int **b1;


	r1 = new int *[y];
	g1 = new int *[y];
	b1 = new int *[y];

	for (int i = 0; i < y; i++) {
		r1[i] = new int[x];
		g1[i] = new int[x];
		b1[i] = new int[x];
		for (int j = 0; j < x; j++) {
			r1[i][j] = -1;
			g1[i][j] = -1;
			b1[i][j] = -1;
		}
	}


	do {
		for (int i = 0; i < y; i++) {
			for (int j = 0; j < x; j++) {
				if (r[i][j] == seedPoint && r1[i][j] == -1) { r1[i][j] = seedPoint; g1[i][j] = seedPoint; b1[i][j] = seedPoint; }
			}
		}
		do {
			bool change = false;
			for (int i = 1; i<y - 1; i++) {
				for (int j = 1; j<x - 1; j++) {
					if (r1[i][j] != -1) {
						if (r1[i + 1][j] == -1 && abs(r1[i][j] - r[i + 1][j]) <= 3) { r1[i + 1][j] = r[i + 1][j]; g1[i + 1][j] = g[i + 1][j]; b1[i + 1][j] = b[i + 1][j]; change = true; }
						if (r1[i - 1][j] == -1 && abs(r1[i][j] - r[i - 1][j]) <= 3) { r1[i - 1][j] = r[i - 1][j]; g1[i - 1][j] = g[i - 1][j]; b1[i - 1][j] = b[i - 1][j]; change = true; }
						if (r1[i][j + 1] == -1 && abs(r1[i][j] - r[i][j + 1]) <= 3) { r1[i][j + 1] = r[i][j + 1]; g1[i][j + 1] = g[i][j + 1]; b1[i][j + 1] = b[i][j + 1]; change = true; }
						if (r1[i][j - 1] == -1 && abs(r1[i][j] - r[i][j - 1]) <= 3) { r1[i][j - 1] = r[i][j - 1]; g1[i][j - 1] = g[i][j - 1]; b1[i][j - 1] = b[i][j - 1]; change = true; }
					}
				}
			}
			//cout<<"tu";
			if (!change) break;
		} while (1);
		//  cout<<"tam";
		if (areSeedsIn(r1, seedPoint)) break;
	} while (1);

	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			if (r1[i][j] == -1) { r[i][j] = backgroundColor; g[i][j] = backgroundColor; b[i][j] = backgroundColor; }
			else { r[i][j] = r1[i][j]; g[i][j] = g1[i][j]; b[i][j] = b1[i][j]; }
		}
	}

	this->save();
}

//TASK 4

complex<double> wyznacznik(int q, int m, int xx) {
	complex<double> result(cos(2 * M_PI*q*m / xx), sin(-1 * 2 * M_PI*q*m / xx));
	return result;
}

complex<double> wyznacznikInverse(int q, int m, int xx) {
	complex<double> result(cos(2 * M_PI*q*m / xx), sin(1 * 2 * M_PI*q*m / xx));
	return result;
}

void Obraz::slowFourier() {
	complex<double> **rC;
	complex<double> **gC;
	complex<double> **bC;

	rC = new complex<double> *[y];
	gC = new complex<double> *[y];
	bC = new complex<double> *[y];

	for (int q = 0; q < y; q++) {
		rC[q] = new complex<double>[x];
		gC[q] = new complex<double>[x];
		bC[q] = new complex<double>[x];
		for (int p = 0; p < x; p++) {

			for (int m = 0; m<y; m++) {
				for (int n = 0; n<x; n++) {
					complex<double> rrC(r[m][n], 0.0);
					complex<double> ggC(g[m][n], 0.0);
					complex<double> bbC(b[m][n], 0.0);
					rC[q][p] += (rrC*wyznacznik(q, m, y)*wyznacznik(p, n, x));
					gC[q][p] += (ggC*wyznacznik(q, m, y)*wyznacznik(p, n, x));
					bC[q][p] += (bbC*wyznacznik(q, m, y)*wyznacznik(p, n, x));
				}
			}
			rC[q][p] /= sqrt(x*y);
			gC[q][p] /= sqrt(x*y);
			bC[q][p] /= sqrt(x*y);
		}
	}


	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			r[i][j] = abs(rC[i][j]);
			g[i][j] = abs(gC[i][j]);
			b[i][j] = abs(bC[i][j]);
		}
	}
	this->save();
}

void Obraz::inverseSlowFourier() {

	complex<double> **rCF;
	complex<double> **gCF;
	complex<double> **bCF;

	rCF = new complex<double> *[y];
	gCF = new complex<double> *[y];
	bCF = new complex<double> *[y];

	for (int q = 0; q < y; q++) {
		rCF[q] = new complex<double>[x];
		gCF[q] = new complex<double>[x];
		bCF[q] = new complex<double>[x];
		for (int p = 0; p < x; p++) {

			for (int m = 0; m<y; m++) {
				for (int n = 0; n<x; n++) {
					complex<double> rrC(r[m][n], 0.0);
					complex<double> ggC(g[m][n], 0.0);
					complex<double> bbC(b[m][n], 0.0);
					rCF[q][p] += (rrC*wyznacznik(q, m, y)*wyznacznik(p, n, x));
					gCF[q][p] += (ggC*wyznacznik(q, m, y)*wyznacznik(p, n, x));
					bCF[q][p] += (bbC*wyznacznik(q, m, y)*wyznacznik(p, n, x));
				}
			}
			rCF[q][p] /= sqrt(x*y);
			gCF[q][p] /= sqrt(x*y);
			bCF[q][p] /= sqrt(x*y);
		}
	}


	complex<double> **rC;
	complex<double> **gC;
	complex<double> **bC;

	rC = new complex<double> *[y];
	gC = new complex<double> *[y];
	bC = new complex<double> *[y];

	for (int q = 0; q < y; q++) {
		rC[q] = new complex<double>[x];
		gC[q] = new complex<double>[x];
		bC[q] = new complex<double>[x];
		for (int p = 0; p < x; p++) {

			for (int m = 0; m<y; m++) {
				for (int n = 0; n<x; n++) {
					rC[q][p] += (rCF[m][n] * wyznacznikInverse(q, m, y)*wyznacznikInverse(p, n, x));
					gC[q][p] += (gCF[m][n] * wyznacznikInverse(q, m, y)*wyznacznikInverse(p, n, x));
					bC[q][p] += (bCF[m][n] * wyznacznikInverse(q, m, y)*wyznacznikInverse(p, n, x));
				}
			}
			rC[q][p] /= sqrt(x*y);
			gC[q][p] /= sqrt(x*y);
			bC[q][p] /= sqrt(x*y);
		}
	}


	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			r[i][j] = abs(rC[i][j]);
			g[i][j] = abs(gC[i][j]);
			b[i][j] = abs(bC[i][j]);
		}
	}
	this->save();
}

void Obraz::fastFourierY(complex<double> *inputR, complex<double> *inputG, complex<double> *inputB, int N, bool isNormal) { //N to dlugosc lub szerokosc obrazu, niewazne w sumie bo obraz musi byc kwadratowy
	if (N == 1) return;

	complex<double> *leftR; complex<double> *leftG; complex<double> *leftB;
	leftR = new complex<double>[N / 2]; leftG = new complex<double>[N / 2]; leftB = new complex<double>[N / 2];
	complex<double> *rightR; complex<double> *rightG; complex<double> *rightB;
	rightR = new complex<double>[N / 2]; rightG = new complex<double>[N / 2]; rightB = new complex<double>[N / 2];
	complex<double> *temp1R; complex<double> *temp1G; complex<double> *temp1B;
	temp1R = new complex<double>[N / 2]; temp1G = new complex<double>[N / 2]; temp1B = new complex<double>[N / 2];
	complex<double> *temp2R; complex<double> *temp2G; complex<double> *temp2B;
	temp2R = new complex<double>[N / 2]; temp2G = new complex<double>[N / 2]; temp2B = new complex<double>[N / 2];

	for (int i = 0; i < (N / 2); i++) {
		leftR[i] = inputR[i]; leftG[i] = inputG[i]; leftB[i] = inputB[i];
		rightR[i] = inputR[i + (N / 2)]; rightG[i] = inputG[i + (N / 2)]; rightB[i] = inputB[i + (N / 2)];
	}

	for (int i = 0; i < (N / 2); i++) {
		complex<double> W;
		if (isNormal)    W = wyznacznik(1, i, N);
		else W = wyznacznikInverse(1, i, N);
		temp1R[i] = leftR[i] + rightR[i]; temp1G[i] = leftG[i] + rightG[i]; temp1B[i] = leftB[i] + rightB[i];
		temp2R[i] = leftR[i] - rightR[i]; temp2G[i] = leftG[i] - rightG[i]; temp2B[i] = leftB[i] - rightB[i];
		temp2R[i] *= W; temp2G[i] *= W; temp2B[i] *= W;
	}
	delete leftR; delete leftG; delete leftB;
	delete rightR; delete rightG; delete rightB;

	fastFourierY(temp1R, temp1G, temp1B, N / 2, isNormal);
	fastFourierY(temp2R, temp2G, temp2B, N / 2, isNormal);

	int i = 0;
	for (int k = 0; i<N; ++k)
	{
		inputR[i] = temp1R[k]; inputG[i] = temp1G[k]; inputB[i] = temp1B[k];
		inputR[i + 1] = temp2R[k]; inputG[i + 1] = temp2G[k]; inputB[i + 1] = temp2B[k];
		i += 2;
	}
	delete temp1R; delete temp1G; delete temp1B;
	delete temp2R; delete temp2G; delete temp2B;
}

void Obraz::swap_quarters() {
	int tempR;
	int tempG;
	int tempB;
	for (int i = 0; i < y / 2; i++)
		for (int j = 0; j < x / 2; j++)
		{
			tempR = r[i][j]; tempG = g[i][j]; tempB = b[i][j];
			r[i][j] = r[i + y / 2][j + x / 2]; g[i][j] = g[i + y / 2][j + x / 2]; b[i][j] = b[i + y / 2][j + x / 2];
			r[i + y / 2][j + x / 2] = tempR; g[i + y / 2][j + x / 2] = tempG; b[i + y / 2][j + x / 2] = tempB;
			tempR = r[i + y / 2][j]; tempG = g[i + y / 2][j]; tempB = b[i + y / 2][j];
			r[i + y / 2][j] = r[i][j + x / 2]; g[i + y / 2][j] = g[i][j + x / 2]; b[i + y / 2][j] = b[i][j + x / 2];
			r[i][j + x / 2] = tempR; g[i][j + x / 2] = tempG; b[i][j + x / 2] = tempB;
		}
}

void Obraz::swap_quarters(complex<double> **inputR, complex<double> **inputG, complex<double> **inputB) {
	complex<double> tempR;
	complex<double> tempG;
	complex<double> tempB;
	for (int i = 0; i < y / 2; i++)
		for (int j = 0; j < x / 2; j++)
		{
			tempR = inputR[i][j]; tempG = inputG[i][j]; tempB = inputB[i][j];
			inputR[i][j] = inputR[i + y / 2][j + x / 2]; inputG[i][j] = inputG[i + y / 2][j + x / 2]; inputB[i][j] = inputB[i + y / 2][j + x / 2];
			inputR[i + y / 2][j + x / 2] = tempR; inputG[i + y / 2][j + x / 2] = tempG; inputB[i + y / 2][j + x / 2] = tempB;
			tempR = inputR[i + y / 2][j]; tempG = inputG[i + y / 2][j]; tempB = inputB[i + y / 2][j];
			inputR[i + y / 2][j] = inputR[i][j + x / 2]; inputG[i + y / 2][j] = inputG[i][j + x / 2]; inputB[i + y / 2][j] = inputB[i][j + x / 2];
			inputR[i][j + x / 2] = tempR; inputG[i][j + x / 2] = tempG; inputB[i][j + x / 2] = tempB;
		}
}

void Obraz::fastFourier() {
	complex<double> **inputR; complex<double> **inputG; complex<double> **inputB;
	inputR = new complex<double> *[y]; inputG = new complex<double> *[y]; inputB = new complex<double> *[y];

	for (int i = 0; i < y; i++) {
		inputR[i] = new complex<double>[x]; inputG[i] = new complex<double>[x]; inputB[i] = new complex<double>[x];
		for (int j = 0; j<x; j++)
		{
			inputR[i][j] = complex<double>(r[i][j], 0.0);
			inputG[i][j] = complex<double>(g[i][j], 0.0);
			inputB[i][j] = complex<double>(b[i][j], 0.0);
		}
	}


	for (int j = 0; j<x; j++)
	{
		complex<double> columnR[y]; complex<double> columnG[y]; complex<double> columnB[y];
		for (int i = 0; i<y; i++) { columnR[i] = inputR[i][j]; columnG[i] = inputG[i][j]; columnB[i] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, true);
		for (int i = 0; i<y; i++) { inputR[i][j] = columnR[i]; inputG[i][j] = columnG[i]; inputB[i][j] = columnB[i]; }
	}

	for (int i = 0; i<y; i++)
	{
		complex<double> columnR[x]; complex<double> columnG[x]; complex<double> columnB[x];
		for (int j = 0; j<x; j++) { columnR[j] = inputR[i][j]; columnG[j] = inputG[i][j]; columnB[j] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, true);
		for (int j = 0; j<x; j++) { inputR[i][j] = columnR[j]; inputG[i][j] = columnG[j]; inputB[i][j] = columnB[j]; }
	}

	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			r[i][j] = abs(inputR[i][j]) / sqrt(y*x);
			g[i][j] = abs(inputG[i][j]) / sqrt(y*x);
			b[i][j] = abs(inputB[i][j]) / sqrt(y*x);
		}
	}
	swap_quarters();
	this->save();
}

void Obraz::fastFourier(complex<double> **inputR, complex<double> **inputG, complex<double> **inputB) {



	for (int j = 0; j<x; j++)
	{
		complex<double> columnR[y]; complex<double> columnG[y]; complex<double> columnB[y];
		for (int i = 0; i<y; i++) { columnR[i] = inputR[i][j]; columnG[i] = inputG[i][j]; columnB[i] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, true);
		for (int i = 0; i<y; i++) { inputR[i][j] = columnR[i]; inputG[i][j] = columnG[i]; inputB[i][j] = columnB[i]; }
	}

	for (int i = 0; i<y; i++)
	{
		complex<double> columnR[x]; complex<double> columnG[x]; complex<double> columnB[x];
		for (int j = 0; j<x; j++) { columnR[j] = inputR[i][j]; columnG[j] = inputG[i][j]; columnB[j] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, true);
		for (int j = 0; j<x; j++) { inputR[i][j] = columnR[j]; inputG[i][j] = columnG[j]; inputB[i][j] = columnB[j]; }
	}

	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			inputR[i][j] /= sqrt(y*x);
			inputG[i][j] /= sqrt(y*x);
			inputB[i][j] /= sqrt(y*x);
		}
	}
	swap_quarters(inputR, inputG, inputB);

}

void Obraz::fastFourierInverse() {
	complex<double> **inputR; complex<double> **inputG; complex<double> **inputB;
	inputR = new complex<double> *[y]; inputG = new complex<double> *[y]; inputB = new complex<double> *[y];

	for (int i = 0; i < y; i++) {
		inputR[i] = new complex<double>[x]; inputG[i] = new complex<double>[x]; inputB[i] = new complex<double>[x];
		for (int j = 0; j<x; j++)
		{
			inputR[i][j] = complex<double>(r[i][j], 0.0);
			inputG[i][j] = complex<double>(g[i][j], 0.0);
			inputB[i][j] = complex<double>(b[i][j], 0.0);
		}
	}


	for (int j = 0; j<x; j++)
	{
		complex<double> columnR[y]; complex<double> columnG[y]; complex<double> columnB[y];
		for (int i = 0; i<y; i++) { columnR[i] = inputR[i][j]; columnG[i] = inputG[i][j]; columnB[i] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, true);
		for (int i = 0; i<y; i++) { inputR[i][j] = columnR[i]; inputG[i][j] = columnG[i]; inputB[i][j] = columnB[i]; }
	}

	for (int i = 0; i<y; i++)
	{
		complex<double> columnR[x]; complex<double> columnG[x]; complex<double> columnB[x];
		for (int j = 0; j<x; j++) { columnR[j] = inputR[i][j]; columnG[j] = inputG[i][j]; columnB[j] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, true);
		for (int j = 0; j<x; j++) { inputR[i][j] = columnR[j]; inputG[i][j] = columnG[j]; inputB[i][j] = columnB[j]; }
	}

	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			inputR[i][j] /= sqrt(y*x);
			inputG[i][j] /= sqrt(y*x);
			inputB[i][j] /= sqrt(y*x);
		}
	}


	for (int j = 0; j<x; j++)
	{
		complex<double> columnR[y]; complex<double> columnG[y]; complex<double> columnB[y];
		for (int i = 0; i<y; i++) { columnR[i] = inputR[i][j]; columnG[i] = inputG[i][j]; columnB[i] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, false);
		for (int i = 0; i<y; i++) { inputR[i][j] = columnR[i]; inputG[i][j] = columnG[i]; inputB[i][j] = columnB[i]; }
	}

	for (int i = 0; i<y; i++)
	{
		complex<double> columnR[x]; complex<double> columnG[x]; complex<double> columnB[x];
		for (int j = 0; j<x; j++) { columnR[j] = inputR[i][j]; columnG[j] = inputG[i][j]; columnB[j] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, false);
		for (int j = 0; j<x; j++) { inputR[i][j] = columnR[j]; inputG[i][j] = columnG[j]; inputB[i][j] = columnB[j]; }
	}

	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			r[i][j] = abs(inputR[i][j]) / sqrt(x*y);
			g[i][j] = abs(inputG[i][j]) / sqrt(x*y);
			b[i][j] = abs(inputB[i][j]) / sqrt(x*y);
		}
	}
	this->save();
}

void Obraz::fastFourierInverse(complex<double> **inputR, complex<double> **inputG, complex<double> **inputB) {


	for (int j = 0; j<x; j++)
	{
		complex<double> columnR[y]; complex<double> columnG[y]; complex<double> columnB[y];
		for (int i = 0; i<y; i++) { columnR[i] = inputR[i][j]; columnG[i] = inputG[i][j]; columnB[i] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, false);
		for (int i = 0; i<y; i++) { inputR[i][j] = columnR[i]; inputG[i][j] = columnG[i]; inputB[i][j] = columnB[i]; }
	}

	for (int i = 0; i<y; i++)
	{
		complex<double> columnR[x]; complex<double> columnG[x]; complex<double> columnB[x];
		for (int j = 0; j<x; j++) { columnR[j] = inputR[i][j]; columnG[j] = inputG[i][j]; columnB[j] = inputB[i][j]; }
		fastFourierY(columnR, columnG, columnB, y, false);
		for (int j = 0; j<x; j++) { inputR[i][j] = columnR[j]; inputG[i][j] = columnG[j]; inputB[i][j] = columnB[j]; }
	}

	for (int i = 0; i < y; i++) {
		for (int j = 0; j < x; j++) {
			r[i][j] = abs(inputR[i][j]) / sqrt(x*y);
			g[i][j] = abs(inputG[i][j]) / sqrt(x*y);
			b[i][j] = abs(inputB[i][j]) / sqrt(x*y);
		}
	}
	//swap_quarters();
	this->save();
}

void Obraz::lowPass(int cut) {
	complex<double> **inputR; complex<double> **inputG; complex<double> **inputB;
	inputR = new complex<double> *[y]; inputG = new complex<double> *[y]; inputB = new complex<double> *[y];
	for (int i = 0; i < y; i++) {
		inputR[i] = new complex<double>[x]; inputG[i] = new complex<double>[x]; inputB[i] = new complex<double>[x];
		for (int j = 0; j<x; j++)
		{
			inputR[i][j] = complex<double>(r[i][j], 0.0);
			inputG[i][j] = complex<double>(g[i][j], 0.0);
			inputB[i][j] = complex<double>(b[i][j], 0.0);
		}
	}
	fastFourier(inputR, inputG, inputB);
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			int d = sqrt((i - y / 2)*(i - y / 2) + (j - x / 2)*(j - x / 2));
			if (d>cut) { inputR[i][j] = 0.0; inputG[i][j] = 0.0; inputB[i][j] = 0.0; }
		}
	}
	fastFourierInverse(inputR, inputG, inputB);
}

void Obraz::highPass(int cut) {
	complex<double> **inputR; complex<double> **inputG; complex<double> **inputB;
	inputR = new complex<double> *[y]; inputG = new complex<double> *[y]; inputB = new complex<double> *[y];
	for (int i = 0; i < y; i++) {
		inputR[i] = new complex<double>[x]; inputG[i] = new complex<double>[x]; inputB[i] = new complex<double>[x];
		for (int j = 0; j<x; j++)
		{
			inputR[i][j] = complex<double>(r[i][j], 0.0);
			inputG[i][j] = complex<double>(g[i][j], 0.0);
			inputB[i][j] = complex<double>(b[i][j], 0.0);
		}
	}
	fastFourier(inputR, inputG, inputB);
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			int d = sqrt((i - y / 2)*(i - y / 2) + (j - x / 2)*(j - x / 2));
			if (d<cut) { inputR[i][j] = 0.0; inputG[i][j] = 0.0; inputB[i][j] = 0.0; }
		}
	}
	fastFourierInverse(inputR, inputG, inputB);
}

void Obraz::bandPass(int min, int max) {
	complex<double> **inputR; complex<double> **inputG; complex<double> **inputB;
	inputR = new complex<double> *[y]; inputG = new complex<double> *[y]; inputB = new complex<double> *[y];
	for (int i = 0; i < y; i++) {
		inputR[i] = new complex<double>[x]; inputG[i] = new complex<double>[x]; inputB[i] = new complex<double>[x];
		for (int j = 0; j<x; j++)
		{
			inputR[i][j] = complex<double>(r[i][j], 0.0);
			inputG[i][j] = complex<double>(g[i][j], 0.0);
			inputB[i][j] = complex<double>(b[i][j], 0.0);
		}
	}
	fastFourier(inputR, inputG, inputB);
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			int d = sqrt((i - y / 2)*(i - y / 2) + (j - x / 2)*(j - x / 2));
			if (d<min || d>max) { inputR[i][j] = 0.0; inputG[i][j] = 0.0; inputB[i][j] = 0.0; }
		}
	}
	fastFourierInverse(inputR, inputG, inputB);
}

void Obraz::bandCut(int min, int max) {
	complex<double> **inputR; complex<double> **inputG; complex<double> **inputB;
	inputR = new complex<double> *[y]; inputG = new complex<double> *[y]; inputB = new complex<double> *[y];
	for (int i = 0; i < y; i++) {
		inputR[i] = new complex<double>[x]; inputG[i] = new complex<double>[x]; inputB[i] = new complex<double>[x];
		for (int j = 0; j<x; j++)
		{
			inputR[i][j] = complex<double>(r[i][j], 0.0);
			inputG[i][j] = complex<double>(g[i][j], 0.0);
			inputB[i][j] = complex<double>(b[i][j], 0.0);
		}
	}
	fastFourier(inputR, inputG, inputB);
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			int d = sqrt((i - y / 2)*(i - y / 2) + (j - x / 2)*(j - x / 2));
			if (d>min && d<max) { inputR[i][j] = 0.0; inputG[i][j] = 0.0; inputB[i][j] = 0.0; }
		}
	}
	fastFourierInverse(inputR, inputG, inputB);
}

void Obraz::highPassEdge(char *argv) {
	cimg_library::CImg<unsigned char> image(argv);
	int **r1;
	r1 = new int *[y];
	for (int i = 0; i < y; i++)
	{
		r1[i] = new int[x];
		for (int j = 0; j < x; j++) {
			r1[i][j] = image.atXYZC(i, j, 0, 0);
		}
	}
	complex<double> **inputR; complex<double> **inputG; complex<double> **inputB;
	inputR = new complex<double> *[y]; inputG = new complex<double> *[y]; inputB = new complex<double> *[y];
	for (int i = 0; i < y; i++) {
		inputR[i] = new complex<double>[x]; inputG[i] = new complex<double>[x]; inputB[i] = new complex<double>[x];
		for (int j = 0; j<x; j++)
		{
			inputR[i][j] = complex<double>(r[i][j], 0.0);
			inputG[i][j] = complex<double>(g[i][j], 0.0);
			inputB[i][j] = complex<double>(b[i][j], 0.0);
		}
	}
	fastFourier(inputR, inputG, inputB);
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			int d = sqrt((i - y / 2)*(i - y / 2) + (j - x / 2)*(j - x / 2));
			if (r1[i][j] == 0) { inputR[i][j] = 0.0; inputG[i][j] = 0.0; inputB[i][j] = 0.0; }
		}
	}
	fastFourierInverse(inputR, inputG, inputB);
}

void Obraz::phaseMod(int k, int l) {
	complex<double> **inputR; complex<double> **inputG; complex<double> **inputB;
	inputR = new complex<double> *[y]; inputG = new complex<double> *[y]; inputB = new complex<double> *[y];
	for (int i = 0; i < y; i++) {
		inputR[i] = new complex<double>[x]; inputG[i] = new complex<double>[x]; inputB[i] = new complex<double>[x];
		for (int j = 0; j<x; j++)
		{
			inputR[i][j] = complex<double>(r[i][j], 0.0);
			inputG[i][j] = complex<double>(g[i][j], 0.0);
			inputB[i][j] = complex<double>(b[i][j], 0.0);
		}
	}
	fastFourier(inputR, inputG, inputB);
	for (int i = 0; i<y; i++) {
		for (int j = 0; j<x; j++) {
			complex<double> faza = exp((-2 * M_PI*i*k) / y + (-2 * M_PI*j*l) / x + (k + l)*M_PI);
			inputR[i][j] *= faza; inputG[i][j] *= faza; inputB[i][j] *= faza;
		}
	}
	fastFourierInverse(inputR, inputG, inputB);
}