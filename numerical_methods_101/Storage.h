#pragma once
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include "gnuplot-iostream.h"

class Storage
{
private:
	enum class_state{MARIAM, ALEXANDER};
	class_state STATE;
public:
	void set_state();
	///PREDEFINED CONSTANTS
	static constexpr double MARIAM_EPSILON = 0.004;
	static constexpr double ALEXANDER_EPSILON = 0.001;
	static constexpr double PI = 3.141592653589793238462643383279502884;
	//given
	std::vector<double> interval;//x... for error 
	std::vector<double> points;//for eitken
	//mid counts
	unsigned int polynomial_degree; // error function
	std::vector<double> nodes; // chebyshev
	std::vector<double> yi;// can count after chebyshev
	//resulting y
	std::vector<double> result;
	//*********************FUNCTIONS DEFINITION ************************//
	void fill_points(); //points
	void fill_interval();//interval
	void fill_yi();//
	// MARS FUNCTION => 1/(10x^4+7x+55) - TODO LATER
	inline double calculate_function(double x, int mariam_flag)
	{return (1 / (pow(x, 2) + 7 * x + 55));}
	//ALEXANDER FUNCTION => ln(x/6)
	inline double calculate_function(double x)
	{return log(x / 6);}

	//FACTORIAL
	inline long long factorial(unsigned n)
	{return ((n > 0) ? (n * factorial(n - 1)) : 1) ;}
	//DERIVATIVE 
	inline double derivative(const unsigned n, const double node_i)
	{return ( abs(pow((-1), n - 1) * factorial(n - 1)) / (pow(node_i / 6, n)) ) * (1.0/6);}

	double max_derivative(int n);
	void error(double const epsilon);
	void ChebyshevPolynomial();

	void update();

	double CallAitken(
		const double point,
		const int nodes_value,
		const unsigned i,
		const unsigned j);
	void Aitken();

	double divided_differences(
		const unsigned polynomial,
		const unsigned i,
		const unsigned j);
	double coefficient(bool flag, double point, unsigned int n);
	void Newton();
	
	Storage();
	~Storage();
};
