#include "Storage.h"

constexpr double Storage::ALEXANDER_EPSILON;
constexpr double Storage::MARIAM_EPSILON;
constexpr double Storage::PI;

void Storage::fill_interval()
{
	double a, b;
	int num;
	std::cout << "input a from range [a;b]" << std::endl;
	std::cin >> a;
	std::cout << "input b from range [a;b]" << std::endl;
	std::cin >> b;
	std::cout << "enter number of elements: " << std::endl;
	std::cin >> num;
	const auto step = ((b - a) / (num - 1));
	auto temp = a;
	for (auto i = 0; i < num; i++)
	{
		interval.push_back(temp);
		temp += step;
	}
}
void Storage::set_state()
{
	bool state;
	std::cout << "Input 0 to check Mariam`s part or input 1 to check Alexander's...";
	std::cin >> state;
	state == 0 ? STATE = MARIAM : STATE = ALEXANDER;
}
void Storage::fill_points()
{
	int point_number;
	std::cout << "enter the number of point you want to check" << std::endl;
	std::cin >> point_number;
	for (auto i = 1; i <= point_number; i++)
	{
		std::cout << "input point #" << (i) << std::endl;
		double inner_temp;
		std::cin >> inner_temp;
		points.push_back(inner_temp);
	}
}
double Storage::max_derivative(int const n)
{
	//finding max derivative on [a,b]
	double fit_derivative = derivative((n + 1), interval.at(0));
	for (unsigned int i = 0; i < interval.size(); i++)
	{
		if (fit_derivative < derivative(n, interval.at(i)))
		{
			fit_derivative = derivative(n, interval.at(i));
		}
	}
	return  fit_derivative;
}
void Storage::error(double const epsilon)
{
	unsigned int n = 0;

	double function_error = 0;
	do
	{
		function_error = max_derivative(n + 1) / (factorial(n + 1) * pow(2, n));
		++n;
	}
	while (function_error > epsilon);
	polynomial_degree = n + 1;
}
void Storage::ChebyshevPolynomial()
{
	double Xk = 0;
	if (interval.at(0) == -1 && interval.at(interval.size() - 1) == 1)
		for (unsigned int k = 0; k < polynomial_degree; k++)
		{
			Xk = cos(PI / 2 * polynomial_degree + ((PI * k) / polynomial_degree));
			nodes.push_back(Xk);
		}
	else
		for (unsigned int k = 0; k < polynomial_degree; k++)
		{
			Xk =  0.5 * (interval.at(0) + interval.at(interval.size() - 1)) +
				( 0.5 * (interval.at(interval.size() - 1) - interval.at(0)) *
				cos(PI / 2 * polynomial_degree + ((PI * k) / polynomial_degree)));

			nodes.push_back(Xk);
		}
}
void Storage::fill_yi()
{
	for(auto element : nodes)
	{
		yi.push_back(calculate_function(element));
	}
}
void Storage::update()
{
	error(ALEXANDER_EPSILON);
	ChebyshevPolynomial();
	sort(nodes.begin(), nodes.end());
	fill_yi();
}
double Storage::CallAitken(
	const double point,
	const int nodes_value,
	const unsigned i, 
	const unsigned j)
{
	double result;
	if (nodes_value > 2)
	{
		result = ((point - nodes.at(i)) * CallAitken(point, nodes_value - 1, i + 1, j) -
			(point - nodes.at(j)) * CallAitken(point, nodes_value - 1, i, j - 1)) /
			(nodes.at(j) - nodes.at(i));
	}
	else
	{
		result = ((point - nodes.at(i)) * yi.at(j)-
			(point - nodes.at(j)) * yi.at(i)) /
			(nodes.at(j) - nodes.at(i));
	}
	return result;
}
void Storage::Aitken()
{
	for (auto point : points)
	{
		result.push_back(CallAitken(point, polynomial_degree, 0, nodes.size() - 1));
		std::cout << "result" << result.at(0);
	}
}
double Storage::divided_differences(
	const unsigned num_of_args,
	const unsigned i,
	const unsigned j)
{
	double result;
	if(num_of_args > 2)
	{
		result												=
		( divided_differences( num_of_args - 1, i, j + 1 )	-
		 divided_differences( num_of_args - 1, i - 1, j ) )	/
		( nodes.at(i) - nodes.at(j) );
	}
	else 
	{
		result = ( yi.at(i) - yi.at(j) ) / ( nodes.at(i) - nodes.at(j) );
	}
	return result;
}
double Storage::coefficient(bool flag, double point, unsigned int n)
{
	double result = 1;
	while (n > 0)
	{
		flag ? result *= (point - nodes.at(n)) :result *= (point - nodes.at(polynomial_degree - n));
		--n;
	}
	return result;

}
void Storage::Newton()
{
	int i = 0;
	for(auto point : points)
	{
		double push = 0.0;
		bool flag = true;
		if(nodes.at(nodes.size() - 1) - nodes.at(0) > point)
		{
			for (unsigned int n = 1; n <= polynomial_degree; n++)
			{
				if (n == 1)
				{
					push = yi.at(0);
				}
				else 
				{
					push += coefficient(flag, point, n - 1) * divided_differences(n, n - 1, 0);
				}
			}
		}
		else
		{
			flag = false;
			for (unsigned int n = 1; n <= polynomial_degree; n--)
			{
				if (n == 1)
				{
					push = yi.at(polynomial_degree - 1);
				}
				else 
				{
					push += coefficient(flag, point, n - 1) * divided_differences(n, polynomial_degree - 1, polynomial_degree - n);
				}
			}
		}
		result.push_back(push);
		std::cout << "result" << result.at(i);
		++i;
	}
}
Storage::Storage()
{
	set_state();
	fill_interval();
	fill_points();
}
Storage::~Storage()
{
}
