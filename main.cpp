#include "CImg.h"
#include "obraz.h"
#include <iostream>
#include <string>

using namespace cimg_library;
using namespace std;

void Help();

enum string_code {
	help,
	brightness,
	contrast,
	negative,
	hflip,
	vflip,
	dflip,
	shrink,
	enlarge,
	median,
	hmean,
	mse,
	pmse,
	stnr,
	pstnr,
	maxDiff,
	histogram,
	hraleigh,
	convolution,
	slineid,
	orobertsi,
	cmean,
	cvariance,
	cstdev,
	cvarcoi,
	casyco,
	cflatco,
	cvarcoii,
	centropy,
	ouolis,
	dilation,
	erosion,
	opening,
	closing,
	HMT,
	m2,
	region,
	slowF,
	invSlowF,
	erosionTask,
	dilationTask,
	hmtTask,
	m1,
	fft,
	fftInv,
	lowPass,
	highPass,
	bandPass,
	bandCut,
	highPassEdge,
	phaseMod,
	kirsch
};

string_code hashit(string const& inString)
{

	if (inString == "--help") return help;
	if (inString == "--brightness") return brightness;
	if (inString == "--contrast") return contrast;
	if (inString == "--negative") return negative;
	if (inString == "--hflip") return hflip;
	if (inString == "--vflip") return vflip;
	if (inString == "--dflip") return dflip;
	if (inString == "--shrink") return shrink;
	if (inString == "--enlarge") return enlarge;
	if (inString == "--median") return median;
	if (inString == "--hmean") return hmean;
	if (inString == "--mse") return mse;
	if (inString == "--pmse") return pmse;
	if (inString == "--stnr") return stnr;
	if (inString == "--pstnr") return pstnr;
	if (inString == "--maxDiff") return maxDiff;
	if (inString == "--histogram") return histogram;
	if (inString == "--hraleigh") return hraleigh;
	if (inString == "--convolution") return convolution;
	if (inString == "--slineid") return slineid;
	if (inString == "--orobertsi") return orobertsi;
	if (inString == "--cmean") return cmean;
	if (inString == "--cvariance") return cvariance;
	if (inString == "--cstdev") return cstdev;
	if (inString == "--cvarcoi") return cvarcoi;
	if (inString == "--casyco") return casyco;
	if (inString == "--cflatco") return cflatco;
	if (inString == "--cvarcoii") return cvarcoii;
	if (inString == "--centropy") return centropy;
	if (inString == "--ouolis") return ouolis;
	if (inString == "--dilation") return dilation;
	if (inString == "--erosion") return erosion;
	if (inString == "--opening") return opening;
	if (inString == "--closing") return closing;
	if (inString == "--HMT") return HMT;
	if (inString == "--m2") return m2;
	if (inString == "--region") return region;
	if (inString == "--slowF") return slowF;
	if (inString == "--invSlowF") return invSlowF;
	if (inString == "--m1") return m1;
	if (inString == "--fft") return fft;
	if (inString == "--fftInv") return fftInv;
	if (inString == "--lowPass") return lowPass;
	if (inString == "--highPass") return highPass;
	if (inString == "--bandPass") return bandPass;
	if (inString == "--bandCut") return bandCut;
	if (inString == "--highPassEdge") return highPassEdge;
	if (inString == "--phaseMod") return phaseMod;
	if (inString == "--kirsch") return kirsch;

	else return help;
}
int main(int argc, char *argv[]) {

	// Obraz obraz("c.bmp");
	// obraz.region(0,255);



	if (argc == 1) {
		string h = "--help";
		*argv = new char[2];
		argv[1] = &h[0];
	}
	switch (hashit(argv[1]))
	{
	case help:

		Help();
		break;

	case brightness:
	{
		if (argc<4) cout << "ERROR!SPECIFY VALUE\n";
		else {
			Obraz obraz(argv[2]);
			obraz.brightness(atoi(argv[3]));
			obraz.save();
		}
		break; }

	case contrast:
	{
		if (argc<4) cout << "ERROR!SPECIFY VALUE\n";
		else {
			Obraz obraz(argv[2]);
			obraz.contrast(atof(argv[3]));
			obraz.save();
		}
		break; }

	case negative: {
		Obraz obraz(argv[2]);
		obraz.negative();
		obraz.save();
		break; }

	case hflip: {
		Obraz obraz(argv[2]);
		obraz.h_flip();
		obraz.save();
		break; }

	case vflip:
	{
		Obraz obraz(argv[2]);
		obraz.v_flip();
		obraz.save();
		break; }

	case dflip:
	{
		Obraz obraz(argv[2]);
		obraz.d_flip();
		obraz.save();
		break; }

	case shrink: {
		if (argc<4) cout << "ERROR!SPECIFY VALUE\n";
		else {
			Obraz obraz(argv[2]);
			obraz.shrink(atoi(argv[3]));
			obraz.save();
		}
		break; }

	case enlarge: {
		if (argc<4) cout << "ERROR!SPECIFY VALUE\n";
		else {
			Obraz obraz(argv[2]);
			obraz.enlarge(atoi(argv[3]));
			obraz.save();
		}
		break; }

	case median:
	{
		Obraz obraz(argv[2]);
		obraz.median();
		obraz.save();
		break;
	}
	case hmean:
	{
		Obraz obraz(argv[2]);
		obraz.hmean();
		obraz.save();
		break; }

	case mse:
	{Obraz obraz(argv[2]);
	Obraz noise(argv[3]);
	Obraz denoised(argv[4]);
	cout << "Mean square error:\n Original to Noise: " << obraz.meanSquareError(noise) << " Original to Denoised: " << obraz.meanSquareError(denoised) << endl;
	break; }

	case pmse:
	{Obraz obraz(argv[2]);
	Obraz noise(argv[3]);
	Obraz denoised(argv[4]);
	cout << "Peak mean square error:\n Original to Noise: " << obraz.peakMeanSquareError(noise) << " Original to Denoised: " << obraz.peakMeanSquareError(denoised) << endl;
	break;
	}
	case stnr:
	{Obraz obraz(argv[2]);
	Obraz noise(argv[3]);
	Obraz denoised(argv[4]);
	cout << "Signal to noise ratio:\n Original to Noise: " << obraz.signalToNoiseRatio(noise) << " Original to Denoised: " << obraz.signalToNoiseRatio(denoised) << endl;
	break;
	}
	case pstnr:
	{Obraz obraz(argv[2]);
	Obraz noise(argv[3]);
	Obraz denoised(argv[4]);
	cout << "Peak signal to noise ratio:\n Original to Noise: " << obraz.peakSignalToNoiseRatio(noise) << " Original to Denoised: " << obraz.peakSignalToNoiseRatio(denoised) << endl;
	break;
	}
	case maxDiff:
	{Obraz obraz(argv[2]);
	Obraz noise(argv[3]);
	Obraz denoised(argv[4]);
	cout << "Maximum difference:\n Original to Noise: " << obraz.maxDiff(noise) << " Original to Denoised: " << obraz.maxDiff(denoised) << endl;
	break;
	}

	case histogram:
	{
		if (argc<4) cout << "ERROR: SPECIFY ARGUMENTS\n";
		else {
			Obraz obraz(argv[2]);
			obraz.histogram(argv[3][0]);
		}
		break; }

	case hraleigh:
	{
		if (argc<4) cout << "ERROR: SPECIFY ARGUMENTS\n";
		else {
			Obraz obraz(argv[2]);
			obraz.hraleigh(atoi(argv[3]), 30);
			obraz.save();
		}
		break;
	}

	case convolution:
	{
		if (argc<12) cout << "ERROR: SPECIFY ARGUMENTS\n";
		else {
			int mask[9];
			for (int i = 3, j = 0; i<12; i++, j++) mask[j] = atoi(argv[i]);
			Obraz obraz(argv[2]);
			obraz.convolution(mask);
		}
		break;
	}

	case slineid:
	{
		int mask[9];
		mask[0] = 2;
		mask[1] = -1;
		mask[2] = -1;
		mask[3] = -1;
		mask[4] = 2;
		mask[5] = -1;
		mask[6] = -1;
		mask[7] = -1;
		mask[8] = 2;
		Obraz obraz(argv[2]);
		obraz.convolution(mask);
		break;
	}

	case orobertsi:
	{
		Obraz obraz(argv[2]);
		obraz.orobertsi();
		break;
	}
	case kirsch:
	{
		Obraz obraz(argv[2]);
		obraz.kirsh();
		break;
	}
	case ouolis:
	{
		Obraz obraz(argv[2]);
		obraz.uolis();
		break;
	}
	case cmean:
	{
		Obraz obraz(argv[2]);
		cout << "Red mean: " << obraz.mean('r') << endl;
		cout << "Blue mean: " << obraz.mean('b') << endl;
		cout << "Green mean: " << obraz.mean('g') << endl;
		break;
	}

	case cvariance:
	{
		Obraz obraz(argv[2]);
		cout << "Red mean: " << obraz.variance('r') << endl;
		cout << "Blue mean: " << obraz.variance('b') << endl;
		cout << "Green mean: " << obraz.variance('g') << endl;
		break;
	}

	case cstdev:
	{
		Obraz obraz(argv[2]);
		cout << "Red mean: " << obraz.standardDeviation('r') << endl;
		cout << "Blue mean: " << obraz.standardDeviation('b') << endl;
		cout << "Green mean: " << obraz.standardDeviation('g') << endl;
		break;
	}
	case cvarcoi:
	{
		Obraz obraz(argv[2]);
		cout << "Red mean: " << obraz.variationCoefficient1('r') << endl;
		cout << "Blue mean: " << obraz.variationCoefficient1('b') << endl;
		cout << "Green mean: " << obraz.variationCoefficient1('g') << endl;
		break;
	}
	case casyco:
	{
		Obraz obraz(argv[2]);
		cout << "Red mean: " << obraz.asymmetryCoefficient('r') << endl;
		cout << "Blue mean: " << obraz.asymmetryCoefficient('b') << endl;
		cout << "Green mean: " << obraz.asymmetryCoefficient('g') << endl;
		break;
	}
	case cflatco:
	{
		Obraz obraz(argv[2]);
		cout << "Red mean: " << obraz.flatteningCoefficient('r') << endl;
		cout << "Blue mean: " << obraz.flatteningCoefficient('b') << endl;
		cout << "Green mean: " << obraz.flatteningCoefficient('g') << endl;
		break;
	}
	case cvarcoii:
	{
		Obraz obraz(argv[2]);
		cout << "Red mean: " << obraz.variationCoefficient2('r') << endl;
		cout << "Blue mean: " << obraz.variationCoefficient2('b') << endl;
		cout << "Green mean: " << obraz.variationCoefficient2('g') << endl;
		break;
	}
	case centropy:
	{
		Obraz obraz(argv[2]);
		cout << "Red mean: " << obraz.informationSourceEntropy('r') << endl;
		cout << "Blue mean: " << obraz.informationSourceEntropy('b') << endl;
		cout << "Green mean: " << obraz.informationSourceEntropy('g') << endl;
		break;
	}
	//===================================TASK 3================================================
	case dilation:
	{
		if (argc<4) cout << "ERROR: SPECIFY ARGUMENTS\n";
		if (atoi(argv[3])<1 || atoi(argv[3])>10) cout << "ERROR: wrong value. Choose 1-10\n";

		Obraz obraz(argv[2]);
		bool tab[3][3];
		obraz.Choose_Structural_Element(atoi(argv[3]), tab);

		obraz.dilation(tab);
		break;
	}

	case erosion:
	{
		if (argc<4) cout << "ERROR: SPECIFY ARGUMENTS\n";
		if (atoi(argv[3])<1 || atoi(argv[3])>10) cout << "ERROR: wrong value. Choose 1-10\n";

		Obraz obraz(argv[2]);
		bool tab[3][3];
		obraz.Choose_Structural_Element(atoi(argv[3]), tab);
		obraz.erosion(tab);
		cout << "FINISHED";
		break;
	}
	case opening:
	{
		if (argc<4) cout << "ERROR: SPECIFY ARGUMENTS\n";
		if (atoi(argv[3])<1 || atoi(argv[3])>10) cout << "ERROR: wrong value. Choose 1-10\n";

		Obraz obraz(argv[2]);
		bool tab[3][3];
		obraz.Choose_Structural_Element(atoi(argv[3]), tab);
		obraz.erosion(tab);
		Obraz edited("edited.bmp");
		edited.dilation(tab);
		cout << "FINISHED";

		break;
	}
	case closing:
	{
		if (argc<4) cout << "ERROR: SPECIFY ARGUMENTS\n";
		if (atoi(argv[3])<1 || atoi(argv[3])>10) cout << "ERROR: wrong value. Choose 1-10\n";

		Obraz obraz(argv[2]);
		bool tab[3][3];
		obraz.Choose_Structural_Element(atoi(argv[3]), tab);
		obraz.dilation(tab);
		Obraz edited("edited.bmp");
		edited.erosion(tab);

		break;
	}
	case HMT:
	{
		if (argc<5) cout << "ERROR: SPECIFY ARGUMENTS\n";
		//	if (atoi(argv[3]) != 11 && atoi(argv[3]) != 12) cout << "ERROR: wrong value. Choose 11 or 12\n";
		//	if ((atoi(argv[3]) == 11 && (atoi(argv[4])<1 || atoi(argv[4])>4)) && (atoi(argv[3]) == 12 && (atoi(argv[4])<1 || atoi(argv[4])>8))) cout << "ERROR: wrong value. Choose adequate number\n";

		Obraz obraz(argv[2]);
		//	obraz.negative();
		//obraz.save();
		//Obraz edited("edited.bmp");

		//bool tab1[3][3];
		//bool tab2[3][3];	
		//obraz.Choose_HMT_Structural_Element(, atoi(argv[4]), tab1, tab2);
		obraz.HMTtransformation(atoi(argv[3]), atoi(argv[4]));// , tab1, tab2);

		break;
	}

	case m2:
	{
		Obraz obraz(argv[2]);
		obraz.m2(atoi(argv[3]), atoi(argv[4]));
		break;
	}
	case region:
	{
		Obraz obraz(argv[2]);
		obraz.region(atoi(argv[3]), atoi(argv[4]));
		break;
	}

	case slowF:
	{
		Obraz obraz(argv[2]);
		obraz.slowFourier();
		break;
	}

	case invSlowF:
	{
		Obraz obraz(argv[2]);
		obraz.inverseSlowFourier();
		break;
	}
	case fft:
	{
		Obraz obraz(argv[2]);
		obraz.fastFourier();
		break;
	}

	case fftInv:
	{
		Obraz obraz(argv[2]);
		obraz.fastFourierInverse();
		break;
	}

	case lowPass:
	{
		Obraz obraz(argv[2]);
		obraz.lowPass(atoi(argv[3]));
		break;
	}

	case highPass:
	{
		Obraz obraz(argv[2]);
		obraz.highPass(atoi(argv[3]));
		break;
	}

	case bandPass:
	{
		Obraz obraz(argv[2]);
		obraz.bandPass(atoi(argv[3]), atoi(argv[4]));
		break;
	}

	case bandCut:
	{
		Obraz obraz(argv[2]);
		obraz.bandCut(atoi(argv[3]), atoi(argv[4]));
		break;
	}

	case highPassEdge:
	{
		Obraz obraz(argv[2]);
		obraz.highPassEdge(argv[3]);
		break;
	}

	case phaseMod:
	{
		Obraz obraz(argv[2]);
		obraz.phaseMod(atoi(argv[3]), atoi(argv[4]));
		break;
	}

	default:
	{
		Help();
	}
}

return 0;
}

void Help()
{
	cout << "Commands:\n";
	cout << "--help \n";
	cout << "--brightness [source file] [value]\n";
	cout << "--contrast [source file] [value]\n";
	cout << "--negative [source file]\n";
	cout << "--hflip [source file] \n";
	cout << "--vflip [source file] \n";
	cout << "--dflip [source file] \n";
	cout << "--shrink [source file] [value] \n";
	cout << "--enlarge [source file] [value] \n";
	cout << "--median [source file] \n";
	cout << "--hmean [source file] \n";
	cout << "--mse [original file] [noised file] [denoised file]\n";
	cout << "--pmse [original file] [noised file] [denoised file]\n";
	cout << "--stnr [original file] [noised file] [denoised file]\n";
	cout << "--pstnr [original file] [noised file] [denoised file]\n";
	cout << "--maxDiff [original file] [noised file] [denoised file]\n";
	cout << "--histogram [original file] [channel (r/g/b)]\n";
	cout << "--hraleigh [original file] [Gmin] [Gmax]\n";
	cout << "--convolution [original file] [M0] [M1] [M2] [M3] [M4] [M5] [M6] [M7] [M8]\n";
	cout << "--slineid [original file]\n";
	cout << "--orobertsi [original file]\n";
	cout << "--cmean [source file]\n";
	cout << "--cvariance [source file]\n";
	cout << "--cstdev [source file]\n";
	cout << "--cvarcoi [source file]\n";
	cout << "--casyco [source file]\n";
	cout << "--cflatco [source file]\n";
	cout << "--cvarcoii [source file]\n";
	cout << "--centropy [source file]\n";

	cout << "--dilation [source file]\n";
	cout << "--erosion [source file] \n";
	cout << "--opening [source file]\n";
	cout << "--closing [source file] \n";
	cout << "--HMT [source file] \n";
	cout << "--m2 [source file] [Yp] [Xp]\n";
	cout << "--slowF [source file]\n";
	cout << "--invSlowF [source file]\n";
}