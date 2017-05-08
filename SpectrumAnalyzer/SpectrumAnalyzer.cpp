#include "stdafx.h"

#define _USE_MATH_DEFINES

#include <sndfile.h>
#include <glut.h>
#include <cmath>
#include <vector>
#include <map>
#include <complex>
#include <iostream>
#include <iomanip>

using namespace std;

const int SPECTRUM_FREQ_STEP = 100;
const int WINDOW_WIDTH = 600;
const int WINDOW_HIEGHT = 400;
const int MAX_CHANNELS = 1;

map<int, float>* spectrumMap = nullptr;

void fft(vector<complex<float>> & samples) {
	int n = samples.size();

	if (n <= 1)
	{
		return;
	}

	vector<complex<float>> a0(n/2), a1(n/2);

	for (int i = 0, j = 0; i < n; i += 2, ++j)
	{
		a0[j] = samples[i];
		a1[j] = samples[i + 1];
	}

	fft(a0);
	fft(a1);

	float angle = 2 * M_PI / n;
	complex<float> w(1), wn(cos(angle), sin(angle));

	for (int i = 0; i < n / 2; ++i)
	{
		samples[i] = a0[i] + w * a1[i];
		samples[i + n / 2] = a0[i] - w * a1[i];
		w *= wn;
	}
}

unsigned int findNearestBiggerPowerOf2(unsigned int x)
{
	x = x - 1;
	x = x | (x >> 1);
	x = x | (x >> 2);
	x = x | (x >> 4);
	x = x | (x >> 8);
	x = x | (x >> 16);

	return x + 1;
}

void fixSamplesVectorLength(vector<complex<float>> & samplesVector)
{
	int initSize = samplesVector.size();
	int nearestBiggerPowerOf2 = findNearestBiggerPowerOf2(initSize);

	for (int i = initSize; i < nearestBiggerPowerOf2; ++i)
	{
		samplesVector.push_back(0);
	}
}

float findMaxAmplitude()
{
	if (spectrumMap == nullptr || spectrumMap->size() < 1)
	{
		throw "Spectrum map is absent";
	}

	float max = spectrumMap->begin()->second;
	for (map<int, float>::iterator it = spectrumMap->begin(); it != spectrumMap->end(); ++it)
	{
		if (it->second > max)
		{
			max = it->second;
		}
	}

	return max;
}

void printTextOnScreen(const float x, const float y, const string s)
{	
	glRasterPos2f(x, y);
	
	for (unsigned int i = 0; i < s.size(); ++i)
	{
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, s[i]);
	}
};

string toStringWithPrecision(const float value, const int precision = 6)
{
	ostringstream out;
	out << setprecision(precision) << value;

	return out.str();
}

const float canvasMaxCoord = 1.0;
const int fullWindowScaleCoef = canvasMaxCoord * 2;

float toPixesCoordX(const int x, const int mapSize)
{
	float stepSize = WINDOW_WIDTH / mapSize;

	return x * stepSize / SPECTRUM_FREQ_STEP;
}

float toCanvasCoordX(const float pixerCoordX, const float shiftX)
{
	float normalizationCoef = canvasMaxCoord / WINDOW_WIDTH;

	return pixerCoordX * normalizationCoef * fullWindowScaleCoef - canvasMaxCoord + shiftX;
}

float normalizeAmplitude(const float ampl, const float maxAmpl)
{
	return ampl / maxAmpl;
}

float toCanvasCoordY(const float normalizedCoordY, const float shiftY)
{
	return normalizedCoordY * (fullWindowScaleCoef - shiftY)  - canvasMaxCoord + shiftY;
}

float tokHz(float hz)
{
	return hz / 1000;
}

void drawXAxis(const float shiftX, const float shiftY)
{
	glBegin(GL_LINES);
	glVertex2f(-canvasMaxCoord, toCanvasCoordY(0, shiftY));
	glVertex2f(canvasMaxCoord, toCanvasCoordY(0, shiftY));
	glEnd();

	float lastCoordX = 0.0;

	int spectrumMapSize = spectrumMap->size();

	const int xAxisFreqSubPointStep = 500;
	const int xAxisFreqMainPointStep = xAxisFreqSubPointStep * 2;

	for (map<int, float>::iterator it = ++spectrumMap->begin(); it != spectrumMap->end(); ++it)
	{
		float coordX = toCanvasCoordX(toPixesCoordX(it->first, spectrumMapSize), shiftX);

		if (it->first % xAxisFreqSubPointStep == 0 
			&& it->first % xAxisFreqMainPointStep != 0)
		{
			glBegin(GL_LINES);
			glVertex2f(coordX, toCanvasCoordY(0, shiftY));
			glVertex2f(coordX, toCanvasCoordY(0, shiftY) - 0.01);
			glEnd();
		}

		if (it->first % xAxisFreqMainPointStep == 0)
		{
			glBegin(GL_LINES);
			glVertex2f(coordX, toCanvasCoordY(0, shiftY));
			glVertex2f(coordX, toCanvasCoordY(0, shiftY) - 0.02);
			glEnd();
			printTextOnScreen(coordX - 0.01, toCanvasCoordY(0, shiftY) - 0.08, toStringWithPrecision(tokHz(it->first), 0));
		}

		lastCoordX = coordX;
	}

	printTextOnScreen(lastCoordX + 0.05, toCanvasCoordY(0, shiftY) - 0.08, "kHz");
}

void drawYAxis(const float shiftX, const float shiftY)
{
	glBegin(GL_LINES);
	glVertex2f(toCanvasCoordX(0, shiftX), canvasMaxCoord);
	glVertex2f(toCanvasCoordX(0, shiftX), -canvasMaxCoord);
	glEnd();

	printTextOnScreen(toCanvasCoordX(0, shiftX) - 0.04, toCanvasCoordY(0, shiftY) - 0.08, toStringWithPrecision(0.0, 2));
	printTextOnScreen(toCanvasCoordX(0, shiftX) - 0.08, canvasMaxCoord - 0.05, "An");

	const float axisMarkStep = 0.2;
	const float axisStep = (canvasMaxCoord * 2 - shiftY) / 5;

	for (float y = -canvasMaxCoord + shiftY + axisStep, mark = axisMarkStep; y < canvasMaxCoord; y += axisStep, mark += axisMarkStep)
	{
		glBegin(GL_LINES);
		glVertex2f(toCanvasCoordX(0, shiftX), y);
		glVertex2f(toCanvasCoordX(0, shiftX) - 0.02, y);
		glEnd();

		printTextOnScreen(toCanvasCoordX(0, shiftX) - 0.08, y - 0.02, toStringWithPrecision(mark, 2));
	}
}

void drawAxes(const float shiftX, const float shiftY)
{
	drawXAxis(shiftX, shiftY);
	drawYAxis(shiftX, shiftY);
}

void resetLineStyle()
{
	glDisable(GL_LINE_STIPPLE);
	glColor3f(1.0, 1.0, 1.0);
}

void enableStrippleLines()
{
	glLineStipple(1, 0xAAAA);
	glEnable(GL_LINE_STIPPLE);
	glColor3f(0.5, 1.0, 0.5);
}

void renderSpectrum()
{
	glClear(GL_COLOR_BUFFER_BIT);

	glLineWidth(1);

	const float shiftX = 0.1;
	const float shiftY = 0.1;

	drawAxes(shiftX, shiftY);
		
	float maxAmplitude = findMaxAmplitude();
		
	int spectrumMapSize = spectrumMap->size();

	enableStrippleLines();

	for (map<int, float>::iterator it = spectrumMap->begin(); it != spectrumMap->end(); ++it)
	{
		glBegin(GL_LINES);
		glVertex2f(toCanvasCoordX(toPixesCoordX(it->first, spectrumMapSize), shiftX), 
			toCanvasCoordY(normalizeAmplitude(it->second, maxAmplitude), shiftY));
		glVertex2f(toCanvasCoordX(toPixesCoordX(it->first, spectrumMapSize), shiftX), toCanvasCoordY(0, shiftY));
		glEnd();
	}
		
	resetLineStyle();

	glFlush();

	glutSwapBuffers();
}

float* getSoundFileSamples(const string soundFilePath, SF_INFO & soundFileInfo)
{
	SNDFILE* rawSoundFile = sf_open(soundFilePath.c_str(), SFM_READ, & soundFileInfo);

	if (rawSoundFile == nullptr)
	{
		throw "Cannot open the file";
	}
	
	if (soundFileInfo.channels > MAX_CHANNELS)
	{
		throw "Cannot process more than 1 channels";
	}

	float* samples = new float[soundFileInfo.frames];

	sf_read_float(rawSoundFile, samples, soundFileInfo.frames);

	sf_close(rawSoundFile);

	return samples;
}

vector<complex<float>> readSoundFileComplexSamples(const string soundFilePath, SF_INFO & soundFileInfo)
{
	float* samples = getSoundFileSamples(soundFilePath, soundFileInfo);

	vector<complex<float>> complexSamples(soundFileInfo.frames);

	for (int i = 0; i < soundFileInfo.frames; ++i)
	{
		complexSamples[i] = complex<float>(samples[i]);
	}

	delete[] samples;

	fixSamplesVectorLength(complexSamples);

	return complexSamples;
}

bool has_suffix(const string & str, const string & suffix)
{
	return str.size() >= suffix.size() &&
		str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool isFilePathValid(const string & soundFilePath)
{
	return has_suffix(soundFilePath, ".wav");
}

void printFileInfo(const SF_INFO soundFileInfo)
{
	cout << "Frames: " << soundFileInfo.frames << endl
		 << "Channels: " << soundFileInfo.channels << endl
		 << "Sample rate: " << soundFileInfo.samplerate << " Hz" << endl;
}

map<int, float>* getSpectrumMap(const string soundFilePath, const int freqStep)
{
	if (!isFilePathValid(soundFilePath))
	{
		throw "Invalid sound file path";
	}

	SF_INFO soundFileInfo;

	vector<complex<float>> complexSamples = readSoundFileComplexSamples(soundFilePath, soundFileInfo);

	printFileInfo(soundFileInfo);

	fft(complexSamples);

	float freqSamplingStep = static_cast<float>(soundFileInfo.samplerate) / complexSamples.size();

	float scaleStep = freqStep / freqSamplingStep;

	int freqLimit = soundFileInfo.samplerate / 2;

	map<int, float>* result = new map<int, float>();

	for (float i = 0.0; i*freqSamplingStep <= freqLimit; i += scaleStep)
	{
		result->insert(pair<int, float>(round(i*freqSamplingStep),
			sqrt(pow(complexSamples[i].real(), 2) + pow(complexSamples[i].imag(), 2))));
	}

	return result;
}

void resize(int width, int height) {
	glutReshapeWindow(WINDOW_WIDTH, WINDOW_HIEGHT);
}

void visualizeSpectrum(string windowTitle, int argc, char* argv[])
{
	glutInit(&argc, argv);

	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
	glutInitWindowPosition(120, 260);
	glutInitWindowSize(WINDOW_WIDTH, WINDOW_HIEGHT);

	glutCreateWindow(windowTitle.c_str());

	glutDisplayFunc(renderSpectrum);
	glutReshapeFunc(resize);
	
	glutMainLoop();
}

void printSpectrumMap(map<int, float>* spectrumMap)
{
	cout << "\nSpectrum map:" << endl;
	cout << "Hz - Amplitude" << endl;

	for (map<int, float>::iterator it = spectrumMap->begin(); it != spectrumMap->end(); ++it)
	{
		cout << it->first << " - " << it->second << endl;
	}
}

int main(int argc, char *argv[])
{
	try {
		if (argc < 2)
		{
			cout << "Usage: SpectrumAnalizer pathToWavFile" << endl;

			return -1;
		}

		string soundFilePath = argv[1];
		
		cout << "\nProcessing...\n" << endl;
		
		spectrumMap = getSpectrumMap(soundFilePath, SPECTRUM_FREQ_STEP);

		printSpectrumMap(spectrumMap);
		
		visualizeSpectrum("Spectrum - " + soundFilePath, argc, argv);

		delete spectrumMap;
	}
	catch(const char* ex)
	{
		cout << ex << endl;
		
		return -1;
	}
					
    return 0;
}