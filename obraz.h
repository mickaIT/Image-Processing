#ifndef Obraz_H
#define Obraz_H
#include "CImg.h"
#include <complex>
#include <vector>
#include <cmath>
using namespace std;


class Obraz {

	int **r;
	int **g;
	int **b;
	int x;
	int y;
	int numbers[256]; //histogram
	bool wasHistogramUsed;

public:

	Obraz(char *name);
	void save();

	void update(int x1, int y1, int **r1, int **g1, int **b1);

	void brightness(int c);

	void contrast(float c);

	void negative();

	void h_flip();

	void v_flip();

	void d_flip();

	void enlarge(int c);

	void shrink(int c);

	int mid(int i, int j, int **r);

	void median();

	int harmonic(int i, int j, int **r);

	void hmean();

	double maxPixel();

	double square(int number);

	double meanSquareError(Obraz image2);

	double peakMeanSquareError(Obraz image2);

	double signalToNoiseRatio(Obraz image2);

	double peakSignalToNoiseRatio(Obraz image2);

	double maxDiff(Obraz image2);

	void histogram(char c);

	void hraleigh(int gmin, float alpha);

	void convolution(int mask[9]);

	void orobertsi();
	float mean(char c);
	float variance(char c);
	float standardDeviation(char c);
	float variationCoefficient1(char c);
	float asymmetryCoefficient(char c);
	float flatteningCoefficient(char c);
	float variationCoefficient2(char c);
	float informationSourceEntropy(char c);
	void makeHistogram(char c);
	int uolisValue(int i, int j, int **r);
	void uolis();

	void Choose_Structural_Element(int structural_element, bool str_el[3][3]);
	void Choose_HMT_Structural_Element(int structural_element, int type, bool str_el[3][3], bool iistr_el[3][3]);
	bool isWhite(int i, int j);
	void dilation(bool str_el[3][3]);
	void erosion(bool str_el[3][3]);
	void HMTtransformation(int structural_element, int z);//, bool str_el[3][3],bool iistr_el[3][3]);
	void mydilation(int **rR, int **gG, int **bB);
	bool isMatrixEqual(int **r1, int **g1, int **b1, int **r2, int **g2, int **b2);
	void m2(int ii, int jj);
	bool areSeedsIn(int **r1, int seedPoint);
	void region(int seedPoint, int backgroundColor);
	int kirschOpCalc(int i, int j, int **r);
	void kirsh();
	void slowFourier();
	void inverseSlowFourier();
	void fastFourierY(complex<double> *inputR, complex<double> *inputG, complex<double> *inputB, int N, bool isNormal);
	void fastFourier();
	void fastFourierInverse();
	void swap_quarters();
	void fastFourier(complex<double> **inputR, complex<double> **inputG, complex<double> **inputB);
	void fastFourierInverse(complex<double> **inputR, complex<double> **inputG, complex<double> **inputB);
	void swap_quarters(complex<double> **inputR, complex<double> **inputG, complex<double> **inputB);
	void lowPass(int cut);
	void highPass(int cut);
	void bandPass(int min, int max);
	void bandCut(int min, int max);
	void highPassEdge(char *argv); // PAMIETAJ ZE TUTAJ MUSISZ UZYC ODPOWIEDNIEGO OBRAZKA I MASKI - WSZYSTKO JEST W POLECENIU
	void phaseMod(int k, int l); // nie wpisuj tu raczej duzych wartosci, dla (1,1) wychodzi prawie lena. poprobuj sobie
};

#endif