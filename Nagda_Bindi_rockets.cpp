#include <iostream>
#include <fstream>
#include <cmath>
#include <functional>
#include <algorithm>
#include <string>

using namespace std;

// Function declaration
double mult(double [3][3], double [3][1], double [3][1]);   
double cross(double[3][1], double[3][1], double[3][1]);

int main() {
	// Declare variables and constants
	double thetag0 = 1.734560396;
	double pi = 4.0 * atan(1.0);						
	double ae = 1.0;									// Units: DU 
	double re = 6378.145;								// Units: km 
	double fttoDu = 4.77881892 * pow(10.0, -8);			// Conversion factor from feet to DU
	double w[3][1] = { { 0 },{ 0 },{ 0.0588336565 } };  // Units: rad/TU
	double e = 0.08182;									// Ellipsoidal earth
	double kms_to_DUTU = 0.1264963205;					// Conversion factor from km/s to DU/TU

	// Input data
	double l = 19.740 * pi/180.0;						    // Latitude units: rad. North is positive.
	double longitude = -155.124 * pi/180.0;					// Longitude units: rad. East is positive.	
	double Height = 7180.0 * fttoDu;					// Units: DU
	double Time = 0317.02;
	double Day = 219.0;												
	double rho = 404.68/re;								// Rho units: DU
	double rhod = 1.38 * kms_to_DUTU;				    // Units: DU/TU
	double el = 20.7 * pi / 180.0;						// Units: rad		
	double eld = 0.07 * (pi / 180.0) * 806.80415;		// Units: rad/TU
	double az = 42.4 * pi / 180.0;						// Units: rad
	double azd = 0.05 * (pi / 180.0) * 806.80415;		// Units: rad/TU

	// Name for the output file
	string FileName = "Bindi_CaseOne.txt";

	// Retrieving digits from Time input for calculation of time in seconds
	double temp = Time * 1000.0;
	double time1 = int(temp) / 100000;
	double temp2 = int(temp) % 100000;
	double time2 = int(temp2) / 1000;
	double temp3 = int(temp) % 1000;
	double time3 = temp3 / 10;
	double TimeS = (time1 * 3600.0) + (time2 * 60.0) + time3;	// Convert time to seconds
	double Day_Fraction = TimeS / 86400;						// Calculate day fraction
	double D = Day + Day_Fraction;
	double thetag = thetag0 + (1.0027379093 * 2 * pi * D);
	double th = thetag + longitude;

	
	char vectors[3] = { 'I', 'J', 'K' };
	
	// To calculate position (rs) and velocity (vs) of radar site in IJK
	double X = fabs(ae / sqrt(1 - pow(e*sin(l), 2)) + Height) * cos(l);
	double Z = fabs((ae*(1-pow(e,2)) / sqrt(1 - pow(e*sin(l), 2)) + Height)) * sin(l);
	
	// position vector (rs) of radar site
	double rs1 = X * cos(th);
	double rs2 = X * sin(th);
	double rs3 = Z;
	double rs[3][1] = { {rs1},{rs2},{rs3} };

	// print out rs vector
	ofstream myfile;
	myfile.open(FileName);
	myfile << "    Created by Bindi Nagda" << endl;
	char buffer[300];
	myfile << endl;
	myfile << "    rs: {";
	for (int k = 0; k < 3; k++) {
		sprintf_s(buffer,"%12.9f %c", rs[k][0], vectors[k]);
		if (k < 2) {
			myfile << buffer << ", ";
		}
		else { myfile << buffer << "} DU" << endl; }
	}

	// vs = cross product of w and rs
	double vs[3][1];  
	cross(w, rs, vs);

	// print out vs vector
	myfile << endl;
	myfile << "    vs: {";
	for (int k = 0; k < 3; k++) {
		sprintf_s(buffer,"%12.9f %c", vs[k][0], vectors[k]);
		if (k < 2) {
			myfile << buffer << ", ";
		}
		else { myfile << buffer << "} DU/TU" << endl; }
	}

	//To calculate position(r) and velocity(v) of satellite in IJK

	// rho vector given in SEZ coordinate system
	double rh1 = -rho * cos(el)*cos(az);
	double rh2 = rho * cos(el)*sin(az);
	double rh3 = rho * sin(el);
	double rh[3][1] = { { rh1 },{ rh2 },{ rh3 } };  				
	
	// inverse Transformation Matrix
	double T[3][3] = { { sin(l)*cos(th), -sin(th), cos(l)*cos(th) },
					 { sin(l)*sin(th), cos(th), cos(l)*sin(th) },
					 { -cos(l), 0.0, sin(l) } };

	double p[3][1];									// p is product of T and rh 
	mult(T, rh, p);									// Need to convert rh into IJK coordinates using Inverse Transformation Matrix

	// declare column vector for position of satellite in IJK system 
	double ri[3][1];

	// Calculate r in IJK system 
	myfile << endl;
	myfile << "    r:  {";
	for (int b = 0; b < 3; b++) {
		ri[b][0] = rs[b][0] + p[b][0];
	}

	// print out r vector
	for (int k = 0; k < 3; k++) {
		sprintf_s(buffer,"%12.9f %c", ri[k][0], vectors[k]);
		if (k < 2) {
			myfile << buffer << ", ";
		}
		else { myfile << buffer << "} DU" << endl; }
	}

	// rho_dot given in SEZ coordinates
	double rhd1 = -(rhod*cos(el)*cos(az)) + (rho*azd*cos(el)*sin(az)) + (rho*eld*sin(el)*cos(az));
	double rhd2 = (rhod*cos(el)*sin(az)) - (rho*eld*sin(el)*sin(az)) + (rho*azd*cos(el)*cos(az));
	double rhd3 = (rhod * sin(el)) + (rho*eld*cos(el));
	double rhd[3][1] = { { rhd1 },{ rhd2 },{ rhd3 } };					

	// H = w x ri
	double H[3][1];		
	cross(w, ri, H);

	// Need to convert rhd into IJK coordinates
	double rhi[3][1];
	mult(T,rhd,rhi);

	// declare column vector for velocity of satellite in IJK system 
	double v[3][1];    

	// print out v vector
	myfile << endl; 
	myfile << "    v:  {";
	for (int c = 0; c < 3; c++) {
		v[c][0] = rhi[c][0] + H[c][0];						// v = rho_dot + (w x ri)
		sprintf_s(buffer,"%12.9f %c", v[c][0], vectors[c]);
		if (c < 2) {
			myfile << buffer << ", ";
		}
		else { myfile << buffer << "} DU/TU" << endl; }
	}

	cout << endl;
	cout << "Output has been printed to a .txt file located in current directory" << endl;
	return 0;
}

// function to multiply 3x3 matrix with a 3x1 column vector
double mult(double T[3][3], double A[3][1], double B[3][1]) {
	for (int k = 0; k < 3; k++) {

		for (int p = 0; p < 1; p++) {
			double num = 0;

			for (int m = 0; m < 3; m++) {
				num = num + (T[k][m] * A[m][p]);
			}
			B[k][p] = num;     				
		}
	}
return 0;
}

// function to evaluate cross product of 3x1 column vectors
double cross(double A[3][1], double B[3][1], double P[3][1]) {
	{
		P[0][0] = (A[1][0] * B[2][0]) - (A[2][0] * B[1][0]);
		P[1][0] = (A[2][0] * B[0][0]) - (A[0][0] * B[2][0]);
		P[2][0] = (A[0][0] * B[1][0]) - (A[1][0] * B[0][0]);
	}
return 0;
}
