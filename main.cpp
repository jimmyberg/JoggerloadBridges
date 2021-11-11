/**
 * @defgroup   MAIN main
 *
 * @brief      This file implements proper alternative prediction of jogger load on bridges for EN 1991-2.
 *
 * It calculates the maximal vibration amplitude for bridges loaded by a passing groups of joggers.
 * It accounts for admittance due to modeshape, and time of travel of the jogger.
 *
 * This incorporates a model of a single span on a bridge with free rotational supports.
 *
 *
 *   O
 *  / \
 *   |                  'pretty bridge'
 *  / \ -->    ----------------------------------
 *  velocity   -  -  -  -  -  -  -  -  -  -  -  -
 *             ----------------------------------
 *             ^                                ^
 *             |                                |
 *
 *
 *
 *
 *
 *
 * The program will ask for length of suspension of the bridge,
 * eigenfrequency of the bridge and damping ratio.
 *
 * It is meant to compliment the jogger load case described by the HIVOSS and EN 1991-2.
 * This simple model accounts for
 *
 * @author     Jimmy van den Berg - Flow Engineering - The Netherlands
 * @date       2021
 */

// #define _DUBUG_

#include <iostream>
#include <cmath>

using namespace std;


// Handy powers of 2 and 3
double power2(double x){return x*x;}
double power3(double x){return x*x*x;}

/**
 * @brief      Second derivative
 */
double ddy(double t, double a, double T){
	return (power3(a)*T*cos(a*t) - power2(a)*sin(a*t) + a*exp(-t/T)/T) / (power2(a*T) + 1);
}

/**
 * @brief      First derivative
 */
double dy(double t, double a, double T){
	return (power2(a)*T*sin(a*t) + a*cos(a*t) - a*exp(-t/T)) / (power2(a*T) + 1);
}

/**
 * @brief      Function describing the amplitude ratio as function of time.
 *
 * @param[in]  t     time
 * @param[in]  a     angular frequency of modeshape at load
 * @param[in]  T     tau factor being the slowness of the rising amplitude = 1-exp(-t/T)
 *                   at constant load.
 *
 * @return     amplitude at time t
 */
double y(double t, double a, double T){
	return (-a*T*cos(a*t) + a*T*exp(-t/T) + sin(a*t)) / (power2(a*T) + 1);
}

/**
 * @brief      Find t where local maxima is for double y(double t, double a,
 *             double T)
 *
 *             Apply newtons method on the first derivative to find local maxima
 *
 * @return     time t of local maxima
 */
double newtonsMethod(double a, double T){
	const double t_max = M_PI/(a);
	double t = t_max*0.75;
	#ifdef _DUBUG_
	cout << t << endl;
	#endif
	for (unsigned int i = 0; i < 6; ++i){
		const double A = ddy(t, a, T);
		const double B = dy(t, a, T) - A*t;
		t = -B/A;
		if(t>t_max)
			t=t_max;
		if(t<0)
			t=0.45*t_max;
		#ifdef _DUBUG_
		cout << t << ", " << A << ", " << B << ", " << A*t+B << endl;
		#endif
	}
	return t;
}

bool plot = false;
bool overrideVelocity = false;

void processArg(const char* arg){
	if(arg[0] == '-'){
		if(arg[1] == 'p' && arg[2] == 0){
			plot = true;
		}
		else if (arg[1] == 'v' && arg[2] == 0){
			overrideVelocity = true;
		}
	}
}

double joggerLoadFactor(double f){
	if(f <= 1.9 || f >= 3.5)
		return 0;
	else if(f < 2.2)
		return (f-1.9)*(2.2-1.9);
	else if (f <= 2.7)
		return 1;
	else
		return -(f-3.5)*(3.5-2.7);;
}

int main(int argc, char** argv){
	double L = 10;
	double v = 3;
	double f = 2.8;
	double z = 0.006;

	for (int index = 1; index < argc; ++index){
		processArg(argv[index]);
	}

	cout << "Resonance frequency at span [Hz] = " << flush;
	cin >> f;
	const float F = joggerLoadFactor(f)*1250;
	cout << '\n';
	cout << "Jogger load [N]                  = " << joggerLoadFactor(f)*1250 << "  per jogger.\n\n";
	cout << "Length of span [m]               = " << flush;
	cin >> L;
	if(overrideVelocity){
		cout << "Velocity joggers                 = " << flush;
		cin >> v;
	}
	else{
		cout << "Assumed velocity jogger [m/s]    = " << v << "\n";
	}
	cout << "Damping of bridge [-]            = " << flush;
	cin >> z;
	cout << '\n';

	const double a = M_PI * v / L;
	const double T = 1/(2*M_PI*f*z);

	if(plot){
		for (int i = 0; i < 100; ++i)
		{
			float t = i*20.0/100;
			cout << t << ", " << y(t, a, T) << '\n';
		}
	}

	const double t = newtonsMethod(a, T);
	const double walkTime = L/v;
	cout << "t_max                            = " << t << " of " << walkTime << " [s] at " << t*100/walkTime << " %" << '\n';
	float y_max_ratio = y(t, a, T);
	cout << "y_max                            = " << y_max_ratio*100 << " % of maximum.\n\n";

	cout << "Generalized mass [kg]            = " << flush;

	float m;
	cin >> m;

	cout << "Maximal acceleration [m/s^2]     = " << y_max_ratio*1250*joggerLoadFactor(f)/(2*m*z) << " per jogger.\n";

	return 0;
}
