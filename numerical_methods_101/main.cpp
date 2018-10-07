#include "Storage.h"

int main()
{
	Gnuplot gp;


	Storage item;
	item.update();
	std::cout << item.polynomial_degree << "\n";

	std::cout << "nodes given by chebyshev:" << std::endl;
	for(unsigned int i = 0; i < item.nodes.size(); i++)
	{
		
		std::cout << "#" << i << " " << item.nodes.at(i) << std::endl;
	}
	item.Aitken();
	item.Newton();


	std::cin.get(); std::cin.get();
	return 0;
}