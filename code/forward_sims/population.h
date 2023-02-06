
#ifndef POPULATION_H
#define POPULATION_H
#include <vector>
#include <boost/random/mersenne_twister.hpp>

class population
{
	public:
		// class constructor
		population(int N);
		static void initialize(double sel, double dom);//set selection and dominance coefficients
		void populate_from(population &p, int N=0);//create next generation
		void populate_from(double prob, int N=0);//create next generation
		int choose_allele();//choose an alleleic state (homo beneficial, hetero or homo deleterious) from the population. I think this function isn't actually used anywhere anymore.
		int binom(int n, double p); //binomial distribution
		int boost_binom(int n, double p);
		int alleleholders[3]; //vector of the three classes of the pop: homo beneficial, hetero and homo deleterious

		~population();  // class destructor
		double frequency;
		double freq();
		double prob();//Expected frequency of deleterious allele post selection
		void mutateup(int n);// n deleterious mutations occur
		void mutatedown(int n);// n beneficial mutations occur
		void clear(); //Fix for beneficial
		void fix(); //Fix for deleterious
		int size;
		int allelenum(); //number of copies of  deleterious allel
		population& operator= (const population &Source); //copy one pop into the other. used to create Europeans from Africans
		boost::mt19937 gent;

	private:

            static double s;
            static double h;
            static double hs;


};

#endif // POPULATION_H
