
# include <cstdlib>

#include <stdio.h> 

# include <iostream>

# include <iomanip>

# include <fstream>

# include <iomanip>

# include <math.h>

# include <ctime>

# include <cstring>



using namespace std;    

# define POPSIZE 50//种群内个体数量

# define MAXGENS 1000//最大的迭代次数

# define NVARS 3//变量个数，即用以表示基因型的bit数

# define PXOVER 0.8//交换率

# define PMUTATION 0.15//突变率


struct genotype

{

	double gene[NVARS];

	double fitness;

	double upper[NVARS];

	double lower[NVARS];

	double rfitness;

	double cfitness;

};



struct genotype population[POPSIZE + 1];

struct genotype newpopulation[POPSIZE + 1];



int main();

void crossover(int& seed);//交叉操作。selects two parents for the single point crossover

void elitist();//stores the best member of the previous generation

void evaluate();//implements the user-defined valuation function

int i4_uniform_ab(int a, int b, int& seed);//returns a scaled pseudorandom I4 between A and B

int Int_uniform_ab(int a, int b);//my own function to get a value between "a" and "b"

void initialize(string filename, int& seed);//initializes the genes within the variables bounds

void keep_the_best();//keeps track of the best member of the population

void mutate(int& seed);//performs a random uniform mutation

double r8_uniform_ab(double a, double b, int& seed);//returns a scaled pseudorandom R8

double Dou_uniform_ab(double a, double b);//return a double value between "a" and "b"

void report(int generation);//reports progress of the simulation

void selector(int& seed);// is the selection function

double round(double number);//returns the integral value that is nearest to x, with halfway cases rounded away from zero

void timestamp();//prints the current YMDHMS date as a time stamp

void Xover(int one, int two, int& seed);//performs crossover of the two selected parents



//****************************************************************************80



int main()


{

	string filename = "simple_ga_input.txt";

	int generation;

	int i;

	int seed;



	timestamp();

	cout << "\n";

	cout << "SIMPLE_GA:\n";

	cout << "  C++ version\n";

	cout << "  A simple example of a genetic algorithm.\n";



	if (NVARS < 2)

	{

		cout << "\n";

		cout << "  The crossover modification will not be available,\n";

		cout << "  since it requires 2 <= NVARS.\n";

	}



	seed = 123456789;

	srand((unsigned)time(NULL));

	initialize(filename, seed);



	evaluate();



	keep_the_best();



	for (generation = 0; generation < MAXGENS; generation++)

	{

		selector(seed);

		crossover(seed);

		mutate(seed);

		report(generation);

		evaluate();

		elitist();

	}



	cout << "\n";

	cout << "  Best member after " << MAXGENS << " generations:\n";

	cout << "\n";



	for (i = 0; i < NVARS; i++)

	{

		cout << "  var(" << i << ") = " << population[POPSIZE].gene[i] << "\n";

	}



	cout << "\n";

	cout << "  Best fitness = " << population[POPSIZE].fitness << "\n";

	//

	//  Terminate.

	//

	cout << "\n";

	cout << "SIMPLE_GA:\n";

	cout << "  Normal end of execution.\n";

	cout << "\n";

	timestamp();



	return 0;

}

//****************************************************************************80



void crossover(int& seed)
{

	const double a = 0.0;

	const double b = 1.0;

	int mem;

	int one;

	int first = 0;

	double x;



	for (mem = 0; mem < POPSIZE; ++mem)

	{

		//x = r8_uniform_ab ( a, b, seed );

		x = Dou_uniform_ab(a, b);



		if (x < PXOVER)

		{

			++first;



			if (first % 2 == 0)

			{

				Xover(one, mem, seed);

			}

			else

			{

				one = mem;

			}



		}

	}

	return;

}



void elitist()
{

	int i;

	double best;

	int best_mem;

	double worst;

	int worst_mem;



	best = population[0].fitness;

	worst = population[0].fitness;



	for (i = 0; i < POPSIZE - 1; ++i)

	{

		if (population[i + 1].fitness < population[i].fitness)

		{



			if (best <= population[i].fitness)

			{

				best = population[i].fitness;

				best_mem = i;

			}



			if (population[i + 1].fitness <= worst)

			{

				worst = population[i + 1].fitness;

				worst_mem = i + 1;

			}



		}

		else

		{



			if (population[i].fitness <= worst)

			{

				worst = population[i].fitness;

				worst_mem = i;

			}



			if (best <= population[i + 1].fitness)

			{

				best = population[i + 1].fitness;

				best_mem = i + 1;

			}



		}



	}

	if (population[POPSIZE].fitness <= best)

	{

		for (i = 0; i < NVARS; i++)

		{

			population[POPSIZE].gene[i] = population[best_mem].gene[i];

		}

		population[POPSIZE].fitness = population[best_mem].fitness;

	}

	else

	{

		for (i = 0; i < NVARS; i++)

		{

			population[worst_mem].gene[i] = population[POPSIZE].gene[i];

		}

		population[worst_mem].fitness = population[POPSIZE].fitness;

	}



	return;

}



void evaluate()

//    The current function is:  x[1]^2-x[1]*x[2]+x[3]

{

	int member;

	int i;

	double x[NVARS + 1];



	for (member = 0; member < POPSIZE; member++)

	{

		for (i = 0; i < NVARS; i++)

		{

			x[i + 1] = population[member].gene[i];

		}

		population[member].fitness = (x[1] * x[1]) - (x[1] * x[2]) + x[3];

	}

	return;

}



int i4_uniform_ab(int a, int b, int& seed)

{

	int c;

	const int i4_huge = 2147483647;

	int k;

	float r;

	int value;



	if (seed == 0)

	{

		cerr << "\n";

		cerr << "I4_UNIFORM_AB - Fatal error!\n";

		cerr << "  Input value of SEED = 0.\n";

		exit(1);

	}

	//

	//  Guarantee A <= B.

	//

	if (b < a)

	{

		c = a;

		a = b;

		b = c;

	}



	k = seed / 127773;



	seed = 16807 * (seed - k * 127773) - k * 2836;



	if (seed < 0)

	{

		seed = seed + i4_huge;

	}



	r = (float)(seed) * 4.656612875E-10;

	//

	//  Scale R to lie between A-0.5 and B+0.5.

	//

	r = (1.0 - r) * ((float)a - 0.5)

		+ r * ((float)b + 0.5);

	//

	//  Use rounding to convert R to an integer between A and B.

	//Returns the integral value that is nearest to x

	value = round(r);//Vs2008中并没有此函数，需要自己实现。是最近取整

	//

	//  Guarantee A <= VALUE <= B.

	//

	if (value < a)

	{

		value = a;

	}

	if (b < value)

	{

		value = b;

	}



	return value;

}

//my new design function of the distribution

int Int_uniform_ab(int a, int b)

{
	int tmp;

	tmp = (rand() % (b - a + 1)) + a;

	return tmp;



}

//****************************************************************************80



void initialize(string filename, int& seed)
{

	int i;

	ifstream input;

	int j;

	double lbound;

	double ubound;



	input.open(filename.c_str());



	if (!input)

	{

		cerr << "\n";

		cerr << "INITIALIZE - Fatal error!\n";

		cerr << "  Cannot open the input file!\n";

		exit(1);

	}

	// 

	//  Initialize variables within the bounds 

	for (i = 0; i < NVARS; i++)

	{

		input >> lbound >> ubound;



		for (j = 0; j < POPSIZE; j++)

		{

			population[j].fitness = 0;

			population[j].rfitness = 0;

			population[j].cfitness = 0;

			population[j].lower[i] = lbound;

			population[j].upper[i] = ubound;

			//population[j].gene[i] = r8_uniform_ab ( lbound, ubound, seed );

			population[j].gene[i] = Dou_uniform_ab(lbound, ubound);

		}

	}



	input.close();



	return;

}

//****************************************************************************80



void keep_the_best()
{

	int cur_best;

	int mem;

	int i;



	cur_best = 0;



	for (mem = 0; mem < POPSIZE; mem++)

	{

		if (population[POPSIZE].fitness < population[mem].fitness)

		{

			cur_best = mem;

			population[POPSIZE].fitness = population[mem].fitness;

		}

	}

	// 

	//  Once the best member in the population is found, copy the genes.

	//

	for (i = 0; i < NVARS; i++)

	{

		population[POPSIZE].gene[i] = population[cur_best].gene[i];

	}



	return;

}

//****************************************************************************80



void mutate(int& seed)
{

	const double a = 0.0;

	const double b = 1.0;

	int i;

	int j;

	double lbound;

	double ubound;

	double x;



	for (i = 0; i < POPSIZE; i++)

	{

		for (j = 0; j < NVARS; j++)

		{

			//x = r8_uniform_ab ( a, b, seed );

			x = Dou_uniform_ab(a, b);

			if (x < PMUTATION)

			{

				lbound = population[i].lower[j];

				ubound = population[i].upper[j];

				//population[i].gene[j] = r8_uniform_ab ( lbound, ubound, seed );

				population[i].gene[j] = Dou_uniform_ab(lbound, ubound);

			}

		}

	}



	return;

}

//****************************************************************************80



double r8_uniform_ab(double a, double b, int& seed)
{

	int i4_huge = 2147483647;

	int k;

	double value;



	if (seed == 0)

	{

		cerr << "\n";

		cerr << "R8_UNIFORM_AB - Fatal error!\n";

		cerr << "  Input value of SEED = 0.\n";

		exit(1);

	}



	k = seed / 127773;



	seed = 16807 * (seed - k * 127773) - k * 2836;



	if (seed < 0)

	{

		seed = seed + i4_huge;

	}



	value = (double)(seed) * 4.656612875E-10;



	value = a + (b - a) * value;



	return value;

}

// my new function of product a value between  "a" and "b"

double Dou_uniform_ab(double a, double b)

{

	double tmp;

	//rand() / double(RAND_MAX)可以生成0~1之间的浮点数

	tmp = a + static_cast<double>(rand()) / RAND_MAX * (b - a);

	return tmp;

}

void report(int generation)
{

	double avg;

	double best_val;

	int i;

	double square_sum;

	double stddev;

	double sum;

	double sum_square;



	if (generation == 0)

	{

		cout << "\n";

		cout << "  Generation       Best            Average       Standard \n";

		cout << "  number           value           fitness       deviation \n";

		cout << "\n";

	}



	sum = 0.0;

	sum_square = 0.0;



	for (i = 0; i < POPSIZE; i++)

	{

		sum = sum + population[i].fitness;

		sum_square = sum_square + population[i].fitness * population[i].fitness;

	}



	avg = sum / (double)POPSIZE;

	square_sum = avg * avg * POPSIZE;

	stddev = sqrt((sum_square - square_sum) / (POPSIZE - 1));

	best_val = population[POPSIZE].fitness;



	cout << "  " << setw(8) << generation

		<< "  " << setw(14) << best_val

		<< "  " << setw(14) << avg

		<< "  " << setw(14) << stddev << "\n";



	return;

}



void selector(int& seed)

{

	const double a = 0.0;

	const double b = 1.0;

	int i;

	int j;

	int mem;

	double p;

	double sum;

	//

	//  Find the total fitness of the population.

	//

	sum = 0.0;

	for (mem = 0; mem < POPSIZE; mem++)

	{

		sum = sum + population[mem].fitness;

	}

	//

	//  Calculate the relative fitness of each member.

	//

	for (mem = 0; mem < POPSIZE; mem++)

	{

		population[mem].rfitness = population[mem].fitness / sum;

	}

	// 

	//  Calculate the cumulative fitness.

	//

	population[0].cfitness = population[0].rfitness;

	for (mem = 1; mem < POPSIZE; mem++)

	{

		population[mem].cfitness = population[mem - 1].cfitness +

			population[mem].rfitness;

	}

	// 

	//  Select survivors using cumulative fitness. 

	//

	for (i = 0; i < POPSIZE; i++)

	{

		//p = r8_uniform_ab ( a, b, seed );

		p = Dou_uniform_ab(a, b);

		if (p < population[0].cfitness)

		{

			newpopulation[i] = population[0];

		}

		else

		{

			for (j = 0; j < POPSIZE; j++)

			{

				if (population[j].cfitness <= p && p < population[j + 1].cfitness)

				{

					newpopulation[i] = population[j + 1];

				}

			}

		}

	}

	// 

	//  Overwrite the old population with the new one.

	//

	for (i = 0; i < POPSIZE; i++)

	{

		population[i] = newpopulation[i];

	}



	return;

}

//****************************************************************************80

double round(double number)

{

	return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);

}

void timestamp()

//    TIMESTAMP prints the current YMDHMS date as a time stamp.

{

# define TIME_SIZE 40



	static char time_buffer[TIME_SIZE];

	const struct tm* tm;

	size_t len;

	time_t now;



	now = time(NULL);

	tm = localtime(&now);

	//将时间格式化

	len = strftime(time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm);



	cout << time_buffer << "\n";



	return;

# undef TIME_SIZE

}

//****************************************************************************80



void Xover(int one, int two, int& seed)
{

	int i;

	int point;

	double t;

	//  Select the crossover point.

	//point = i4_uniform_ab ( 0, NVARS - 1, seed );

	point = Int_uniform_ab(0, NVARS - 1);

	//

	//  Swap genes in positions 0 through POINT-1.

	//

	for (i = 0; i < point; i++)

	{

		t = population[one].gene[i];

		population[one].gene[i] = population[two].gene[i];

		population[two].gene[i] = t;

	}



	return;

}


