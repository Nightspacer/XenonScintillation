// Example program
#include <iostream>
#include <string>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <random>
#include <sstream>


double driftVelocity(double x, double t);

double tFinal = 1 * pow(10, -5);//Total time [s]
double gasPoint = (1/10.0);
const double mTocm = 100;


int main()
{
	double f = 5 * pow(10, 6);//Frequency [1/s]
	double Vmax = 5000;//Voltage [V]
	double d = 1.0;//Separation distance[m]
	double q = -(1.602*pow(10, -19));//Charge of electron [C]
	double m = 9.10938356*pow(10, -31);//Mass of electron [kg]

	const int n = 10000000;
	const int points = 1000;
	double dt = tFinal / n;
	double p[n / points] = {};//position
	double v[n / points] = {};//velocity
	double a[n / points] = {};//acceleration
	double e[n / points] = {};//electric field
	double T[n / points] = {};//Time
	double AP[n / points] = {};//summed photons
	double Po[n / points] = {};//summed photons
	double L[n / points] = {};//total distance
	double K[n / points] = {};//photons per ecm
//	double Ls[n / points] = {};
	double voltageDifference;//voltage at a given time
	double electricField;//electric field at a given time
	double x;//absolute value of electric field
	double dragForce;
	double force;
	double acceleration = 0;
	double velocity = 0;
	double position = 0;
	double t;
	double dx;
	double totaldistance = 0;

	double AttenuationProbability;
	double Probability;
	double MeanFreeLiquid = 1.0;//average absorption distance in liquid [m]
	double MeanFreeGas = MeanFreeLiquid / 5.671 * 3100;//average absorption distance in gas [m]

	int NumParticles = 10;

	double A = .137;
	double b = 177.0;
	double y = 45.7;
	double P = 2.0;
	double PhotonsPer_e_cm = 0;
	double PhotonsTotal = 0;
	double PhotonArray[n / points] = {};
	double PhotonTotalArray[n / points] = {};
	double N[n / points] = {};
	double Photons = 0;

	int i;
	int j;
	int NumParticlesReduction = 0;
	const int z = 3;
	const int SumNum = 50;

	const double PI = 3.141592653589793;

    std::default_random_engine generator;
	srand(time(NULL));

	for (i = 0; i<n; i++){
		t = i*dt;
		voltageDifference = Vmax*cos(2*PI*f*t);
		electricField = voltageDifference / d;
		x = sqrt(pow(electricField,2));
		switch (z){
		case 1:
			dragForce = -(1.602*pow(10, -19))*x*velocity / driftVelocity(x, t);
			force = dragForce + electricField*q;
			acceleration = force / m;
			velocity = velocity + acceleration*dt;
			break;
		case 2:
			dragForce = -(1.602*pow(10, -19))*x*velocity / driftVelocity(x, t);
			force = dragForce + electricField*q;
			acceleration = force / m;
			velocity = velocity + acceleration*dt;
			if (abs(velocity)>driftVelocity(x, t)){
				velocity = -copysign(1.0, electricField)*driftVelocity(x, t);
			}
			break;
		case 3:
			velocity = -copysign(1.0, electricField)*driftVelocity(x, t);
			break;
		}
		position = position + velocity*dt*mTocm;
		totaldistance = totaldistance + sqrt(pow(velocity*dt*mTocm,2));
		dx = sqrt(pow(velocity*dt*mTocm,2));
		if (NumParticles > 0){
			NumParticlesReduction = 0;
			if (t > (tFinal*gasPoint)){
				for (j = 0; j < NumParticles; j++){
					AttenuationProbability = exp(-(dx / MeanFreeGas));
					std::uniform_real_distribution<double> distribution(0,1);
					Probability = distribution(generator);
					//Probability = ((rand()*rand()) % 100000000) / 100000000.0;
					if (AttenuationProbability > Probability){
						PhotonsPer_e_cm = A*x - b*P - y;
						if (PhotonsPer_e_cm >= 0){
							PhotonsTotal = PhotonsTotal + PhotonsPer_e_cm*dx;
							Photons = PhotonsPer_e_cm*dx;
						}
						else{
							Photons = 0;
						}
					}
					else{
						NumParticlesReduction = NumParticlesReduction + 1;
					}
				}
				NumParticles = NumParticles - NumParticlesReduction;
			}
			else{
				for (j = 0; j < NumParticles; j++){
					AttenuationProbability = exp(-dx / MeanFreeLiquid);
					std::uniform_real_distribution<double> distribution(0,1);
					Probability = distribution(generator);
					//Probability = ((rand()*rand()) % 100000000) / 100000000.0;
					if (AttenuationProbability < Probability){
						NumParticlesReduction = NumParticlesReduction + 1;
					}
				}
				NumParticles = NumParticles - NumParticlesReduction;
			}
		}
		else{
			NumParticles = 0;
		}
		if (i % points == 0){
			N[i / points] = NumParticles;
			e[i / points] = electricField;
			a[i / points] = acceleration;
			v[i / points] = velocity;
			p[i / points] = position;
			T[i / points] = t;
			L[i / points] = totaldistance;
			K[i / points] = PhotonsPer_e_cm;
//			Ls[i / points] = dx;
			PhotonArray[i / points] = Photons;
			PhotonTotalArray[i / points] = PhotonsTotal;
			if (i >= SumNum * points){
				AP[i / points] = PhotonTotalArray[i/points]-PhotonTotalArray[i / points-SumNum];
				if ((i % (SumNum * points))==0){
					Po[i / points] = PhotonTotalArray[i / points] - PhotonTotalArray[i / points - SumNum];
				}
			}
		}
	}

	std::ofstream myfile;
	myfile.open("1cm_MFP_liquid_20_sum_500_test_21.txt");
        myfile <<"TIme" << '\t' << "Electric Field" << '\t' << "Acceleration" << '\t' << "Velocity" << '\t' << "Position" << '\t' << "Total Distance" << '\t' << "Photons" << '\t' << "Total Photons" << '\t' << "Electrons" << '\t' << "summed photons in period" << '\t' << "Summed photons in period, with period sample" << '\t' << "Photons per ecm" << '\n';
	for (i = 0; i < (n / points); i++){
		myfile << T[i] << '\t' << e[i] << '\t' << a[i] << '\t' << v[i] << '\t' << p[i] << '\t' << L[i] << '\t' << PhotonArray[i] << '\t' << PhotonTotalArray[i] << '\t' << N[i] << '\t' << AP[i] << '\t' << Po[i] << '\t' << K[i] << '\n';
	}
	myfile.close();

	std::ofstream yourfile;
	yourfile.open("FinalInput_test_10.txt");
	for (i = 0; i < (((n / points) / SumNum)); i++){
		yourfile << T[i * SumNum] << '\t' << Po[i * SumNum] << '\t' << Po[i * SumNum]/1000.0 << '\t' << Po[i*SumNum]/N[i*SumNum] << '\n';
	}
	yourfile.close();
}

double driftVelocity(double x, double t){
	//double tfinal = 5 * pow(10, -5);
	double y = mTocm * x;
	if (t >(tFinal*gasPoint)){
		if (y < 3000){
			return(868.807 / 3000 * y);
		}
		else if (y < 500000){
			return(768.95 + 3.5823*pow(10, -2)*y - 8.8199*pow(10, -7)*pow(y, 2) + 1.2298*pow(10, -11)*pow(y, 3) - 8.1118*pow(10, -17)*pow(y, 4) + 2.9328*pow(10, -22)*pow(y, 5) - 5.9684*pow(10, -28)*pow(y, 6) + 6.41958*pow(10, -34)*pow(y, 7) - 2.8394*pow(10, -40)*pow(y, 8));
		}
		else{
			return(0.0278642*y + 11026.7 - 0.0278642 * 500000);
		}
	}
	else{
		return(2000 + (1 / 9.0) * (x - 1000));
    }
}
