// Example program
#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <sstream>
#include <random>

int main()
{
	int i = 0;
	int j = 0;
	int k = 0;

	double A, B, C, D, a, b, c;
	const int SizeOfArray = 200;
	const double TimeOfsection = 5.0; //microsecond
	const double MeasureTime = .1;
	int MeasureSteps = int(TimeOfsection / MeasureTime);
	double t;
	double PhotonsPerTime[SizeOfArray] = {};
	const int SizeOfOutArray = 10000;
	double Times[SizeOfOutArray] = {};
	double Photons[SizeOfOutArray] = {};
	double Photons_1[SizeOfOutArray] = {};
	double Times_2[SizeOfOutArray + 500] = {};
	double Photons_2[SizeOfOutArray + 500] = {};
	double Efficeny = 0.1;

	std::ifstream infile("ReducedHistoryofVariables.txt");
	std::string line;
	while (std::getline(infile, line)){
		std::istringstream iss(line);
		iss >> A >> B >> C >> D;
		PhotonsPerTime[i] = B;
		i++;
	}
	srand(time(NULL));

	for (i = 0; i < (SizeOfArray + 1); i++){
		if (PhotonsPerTime[i] > 0){
			for (j = 0; j < PhotonsPerTime[i]; j++){
				a = rand();
				b = (TimeOfsection / RAND_MAX);
				c = (i - 1)*TimeOfsection;
				t = (a*b) + c;
				while (t >= (k*MeasureTime)){
					k++;
				}
				Times[k - 1] = Times[k - 1] + 1;
				k = 0;
			}
		}
	}

	for (i = 0; i < SizeOfOutArray; i++){
		for (j = 0; j < Times[i]; j++){
			if (Efficeny >= (rand() % 100) / 100.0){
				Photons[i] = Photons[i] + 1;
			}
		}
	}

	std::default_random_engine generator;


	for (i = 0; i < SizeOfOutArray; i++){
		if (Photons[i] > 0){
			std::normal_distribution<double> distribution(Photons[i], .5*sqrt(Times[i]));
			Photons_1[i] = distribution(generator);
		}
		else{
			Photons_1[i] = 0;
		}
	}

	for (i = 0; i < SizeOfOutArray; i++){
		Photons_2[i + 500] = Photons_1[i];
	}

	std::normal_distribution<double> distribuition_2(0.0, .1);
	for (i = 0; i < (SizeOfOutArray + 500); i++){
		if (i == 500){
			Photons_2[i] = Photons_2[i] + distribuition_2(generator) + 20;
		}
		else{
			Photons_2[i] = Photons_2[i] + distribuition_2(generator);
		}

	}

	std::ofstream yourfile;
	yourfile.open("FinalOutputTest1.txt");
	for (i = 0; i < (SizeOfArray*MeasureSteps); i++){
		yourfile << i*MeasureTime << '\t' << Times[i] << '\t' << Photons[i] << '\t' << Photons_1[i] <<  '\n';
	}
	yourfile.close();

	std::ofstream theirfile;
	theirfile.open("GraphThis1.txt");
	for (i = 0; i < ((SizeOfArray+9)*MeasureSteps); i++){
		theirfile << i*MeasureTime << '\t' << Photons_2[i] << '\n';
	}

	theirfile.close();
}
