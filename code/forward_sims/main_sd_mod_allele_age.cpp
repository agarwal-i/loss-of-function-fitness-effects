//Modified by Eduardo Amorim (guerraamorim@gmail.com) from the original of Yuval Simons (Simons et al. 2014).
//Modified by Zach Fuller from Eduardo Amorim (Amorim et al. 2017) and Yuval Simons (Simons et al. 2014)

#include <chrono>
#include <random>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "population.h"
#include "BRand.hpp"
#include <stdint.h>
#include <map>
#include <time.h>
#include <sstream>
#include <fstream>
#include <string>
#include <iostream>
  using std::cout;
  using std::endl;
#include <iomanip>
  using std::setprecision;
#include <cstdlib>
  using std::atoi;
  using std::atof;
#include <ctime>
  using std::time;
#include <vector>
  using std::vector;
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/poisson_distribution.hpp>
  using boost::poisson_distribution;
#include <boost/random/variate_generator.hpp>
  using boost::variate_generator;

#define SPLIT 1000
#define Nzero 10000
#define tau 2040
#define taujump 56000
//#define SEED 2710
#define TAU 100
//#define Urate 0.0000000125//0.00000025
#define ratio 1
int Ne[57]= {14448,14068,14068,14464,14464,15208,15208,16256,16256,17618,17618,19347,19347,21534,21534,24236,24236,27367,27367,30416,30416,32060,32060,31284,29404,26686,23261,18990,16490,16490,12958,12958,9827,9827,7477,7477,5791,5791,4670,4670,3841,3841,3372,3372,3287,3359,3570,4095,4713,5661,7540,11375,14310,13292,14522,1000000,5000000};
int T[57]={55940,51395,47457,43984,40877,38067,35501,33141,30956,28922,27018,25231,23545,21951,20439,19000,17628,16318,15063,13859,12702,11590,10517,9482,8483,7516,6580,5672,5520,5156,4817,4500,4203,3922,3656,3404,3165,2936,2718,2509,2308,2116,1930,1752,1579,1413,1252,1096,945,798,656,517,383,252,124,50,0};

boost::mt19937 gent;

double lognormal();
double Uvar;
double ufactor;
int pid_seed;
int nfinal,RUNS,Ntau=14448;
int idx=0;
double sel,DOM,mutU;
int birthdate,deathdate;
using namespace std;

struct freqs
{
  double freq0;
  double freq1;

  freqs(const double a=0,const double b=0) :
    freq0(a),freq1(b) {}
};

class valueComp
{
public:
  bool operator()(const freqs& A,
          const freqs& B)
    const
    { if (A.freq0!=B.freq0)
         return A.freq0<B.freq0;
      else
         return A.freq1<B.freq1; }
};

char* filename(char *);
char* seriesfilename(char * str);
int popsize(int gen, int i);
int poi(double l);
int poicond1(double l);
int boost_poi(double l);
int boost_poicond1(double l);

//double lognormal(double mutU);

//

void Print(const vector<int>& v);

void Print (const vector<int>& v){
  //vector<int> v;
  for (int i=0; i<v.size();i++){
    cout << v[i] << endl;
  }
}

int main(int argc, char *argv[])
{

    //int initgen=150000,stopover;
    int initgen=taujump,stopover;
    double m=0.00015;
    if (argc==8) // User may run the script with "./a.out RUNS sel DOM", where RUNS = number of runs, sel = selective coefficient (in this work always set to 1), and DOM = heterozygote effect (in this work set to 0 or 1%)
       {RUNS=atoi(argv[1]);
       sel=atof(argv[2]);
       DOM=atof(argv[3]);
       mutU=atof(argv[4]);
       ufactor=atof(argv[5]);
       pid_seed=atof(argv[7]);

       }
    else         // User may run the script with "./a.out" and the above mentioned parameters will be asked in the command line prompt.
        {
        cout<<"enter number of runs:";
        cin>>RUNS;
        cout<<"enter selective coefficient:";
        cin>>sel;
        cout<<"enter dominance:";
        cin>>DOM;
        cout<<"enter per gene mutation rate:";
        cin>>mutU;
        cout<<"enter factor for back mutation:";
        cin>>ufactor;
        cout<<"enter process id:";
        cin>>pid_seed;
        }
    stopover=RUNS;
    //cout<<"RUNS:"<<RUNS<<" tau:"<<tau<<" sel:"<<sel<<" initgen:"
    //<<initgen<<" stopover:"<<stopover<<" dominace:"<<DOM<<" m:"<<m<<"\n";
    gent.seed(time(NULL)+pid_seed);
    BRand::Controller.seed(time(NULL)+pid_seed);

    std::uniform_real_distribution<double> mt_rand{0.0, 1.0};

    population::initialize(sel,DOM);
    population* pops[2];
    pops[0]= new population(popsize(0,0));
    pops[1]= new population(popsize(0,0));

    int euroflag=0;
    double halftheta=2*popsize(initgen,0)*mutU;
    double p1=1./(1.+ratio*exp((2*popsize(initgen,0)-1)*sel)*(1-(2*DOM-1)*(2*popsize(initgen,0)-1)*sel*sel/(12*popsize(initgen,0))) );
    //cout<<"P1:"<<halftheta<<",popsize:"<<popsize(initgen,0)<<"\n";
    //map<int,int> count[2],counts[2][2],countd[2],taucount,taucountd;
    //map<freqs,int,valueComp> joint,joints[2],jointd;

    int gen,thispop;
    double baseup=(1./halftheta);
    double basedown=(1./(ratio*halftheta));
    double oneovertwoU=(1./(2*mutU));
    double oneovertwoUratio=(1./(2*mutU*ratio));
    //int initialstate,zeroruns=0,oneruns=0;
    int skip=0,stretch;
    time_t tt;
    struct tm *tim;
    char library[100],comm[200],totalfile[200],totalres[200];
    char outputNameEur[200], outputNameAfr[200], outputNameMut[200];
    double Ep,s1,S2N,Ep2,Ex,Exderived,s12,Ex2,Ex4,Dx,Dx2;
    double europ, euros1, eurox, euroS2N;
    double totfreq;
    int age=0;
    int tauallele,before, after;
    freqs f;

    ofstream myfile,results,afrDistrib,eurDistrib,mutDistribution;

    for(int run=0;run<RUNS;run++)
            {
            gen=initgen;
            idx=0;
            //double Uvar = lognormal(mutU);
            double Uvar = mutU;
            double halfthetaVAR = 2*popsize(initgen,0)*Uvar;
            double baseupVAR = (1./halfthetaVAR);
            double basedownVAR = (1./halfthetaVAR);
            //cout << basedownVAR << " " << baseupVAR << "\n";
            pops[0]->size=popsize(gen,0); //Initializing the population size for the beginning of each run with deleterious allele absent from the population

                pops[0]->clear();
                //initialstate=0;
                //zeroruns++;

            euroflag=0; //A flag indicating if the African-European population split has occured (0 = it hasn't occurred)

            while (gen>0)
                    {
          //if (gen==tau) //tau is the generation of the out of Africa exodus.
            //          {
                        //euroflag=1;
                        //*pops[1]=*pops[0]; //Creating the European population from the African population.
                        //taufreq=pops[0]->freq(); //Keeping track of deleterious allele frequencies at the split for calculation of the change in frequency since the split
              //        }
		      if (gen==920) m=0.000025;
          if (euroflag==1) //Keep the derived allele derived
      {totfreq=0.5*(pops[0]->freq()+pops[1]->freq());
        if (totfreq==1) {pops[0]->clear();pops[1]->clear();totfreq=0.0;}
      }
          else
      {totfreq=pops[0]->freq();
        if (totfreq==1) {pops[0]->clear();totfreq=0.0;}
      }

		      for(int i=0;i<(1);i++)
                    {

//                      printf("%d ", pops[i]->allelenum());
		      
                      if ((pops[i]->allelenum()==0))
                       {
                          deathdate=gen;
                          //printf("%d ",deathdate);
                        }


                      if (pops[i]->allelenum()==0)
                       {
                        pops[i]->clear();
                        //pops[i]->mutateup(1.0);
                        pops[i]->mutateup(std::min(1,boost_poi( Uvar*(2*pops[i]->size) )));
                       }

                      else if (pops[i]->allelenum()>0)
                       {
                        pops[i]->mutateup(0.0);
                       }

		     //if (((pops[i]->allelenum()>0)&&(pops[i]->allelenum()<(2*pops[i]->size)))||(gen<=taujump)) //If the allele is segregating or we are in recent history (last 5920 generations), introduce mutations 
                       // {
                           // if ((pops[i]->allelenum()==1))
			//	{
			  //        birthdate=gen;
                                  //printf("%d ",birthdate);
			//	}
                            //pops[i]->mutateup(std::min(1,boost_poi( Uvar*(2*pops[i]->size-pops[i]->allelenum()) ))) ;
                            //pops[i]->mutatedown(boost_poi( ufactor*Uvar*pops[i]->allelenum())) ;
                            //cout << Uvar << " " << Uvar*(2*pops[i]->size-pops[i]->allelenum()) << " " << Uvar*pops[i]->allelenum() << "\n";
                        //}
		      //else if (pops[i]->allelenum()==0) //If the ancestral allele is fixed calculate when will the next deleterious allele appear and advance population state to that point in time (unless that point in time is in recent history; then just set time to initial history).
                        //{
                             //pops[i]->clear();
                             //cout << BRand::Controller.nextOpened() << endl;
                             //cout<< "mt:" << mt_rand(gent) << endl;
                             //gen-=int(-log( mt_rand(gent))*baseupVAR);
                             //if (gen<=taujump)
                             //   gen=taujump+1;
                             //else
                                // pops[i]->size=popsize(gen,i);
                                 //pops[i]->clear(); 
                                 //pops[i]->mutateup(100.0);
                                // pops[i]->mutateup(std::min(1,poicond1(2*Uvar*pops[i]->size))); 

                        //}
                        //else if (pops[i]->allelenum()==(2*pops[i]->size))//If the deleterious allele is fixed calculate when will the next deleterious allele appear and advance population state to that point in time (unless that point in time is in recent history and then just set time to initial history). Note that with s = 1 this cannot happen
                        //{
                             //pops[i]->clear();
                            // gen-=int(-log(mt_rand(gent))*basedownVAR);
                             //if (gen<=taujump)
                               // gen=taujump+1;
                             //else
                               //  {pops[i]->size=popsize(gen,i);
                               //  pops[i]->fix();
                               //  pops[i]->mutatedown(poicond1(2*ufactor*Uvar*pops[i]->size)); }

                        //}

		      if (euroflag==0) //If a split between Af and Eur populations has occured apply migration
                           pops[i]->populate_from(pops[i]->prob(),popsize(gen-1,i));
                        else
                           pops[i]->populate_from((1-m)*(pops[i]->prob())+m*(pops[1-i]->prob()),popsize(gen-1,i));
                    }
          //std::stringstream result;
          //std::copy(pops[0]->alleleholders.begin(), pops[0]->alleleholders.end(), std::ostream_iterator<int>(result, " "));

		      gen--; //Next generation please

}

cout << gen << " ";
for (int i = 3 - 1; i >= 0; i--)
    cout << pops[0]->alleleholders[i]<< " ";
for (int i = 3 - 1; i >= 0; i--)
    cout << pops[1]->alleleholders[i] << " ";
cout << Uvar;
cout << " " << sel << " " << DOM << " " << deathdate << " " << birthdate;
cout << "\n";

    }
}


int popsize(int gen, int i) //Calculates the size of population i (0=African, 1=European) at generation gen according to Schiffles & Durbin's model
{
  //cout << gen << " " << idx << "\n" ;
  if (gen<=T[idx]) idx+=1;
  //return 150000;
  return Ne[idx];



}

double lognormal(double mutU)
{

  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  gent.seed(seed);
  //std::default_random_engine generator (seed);

  // The following lines have to be uncommented or commented prior to compilation according to the type of mutation that is going to be simulated.
  // There are 4 types (CpGti, CpGtv, nonCpGti and nonCpGtv) and the genomic average, based on Kong et al. (2012) - See Methods for detail
  // We also apply the correction from Harpak et al. (2016) - See Methods for detail.

  //std::normal_distribution<double> distribution(-7.324836926,0.57); //CpGti
  //std::normal_distribution<double> distribution(-8.392236341,0.57); //CpGtv
  //std::normal_distribution<double> distribution(-8.583066473,0.57); //nonCpGti
  //std::normal_distribution<double> distribution(-8.798867103,0.57); //nonCpGtv
  //std::normal_distribution<double> distribution(-5.906195,0.57); //PRDM9
  std::normal_distribution<double> distribution(mutU,0.57);

  double u = distribution(gent);
  return pow(10,u);
}

int poi(double l) //Regular Poisson random variate

{double p=1,expl=exp(-l);
int k=0;
while (p>=expl)
{k++;
p=p*BRand::Controller.nextClosed();
}
return k-1;

}

int boost_poi(double l) //Regular Poisson random variate

{
    //cout << l << endl;
    if (l==0.0){
      return 0.0;
    }

    else{
      boost::random::poisson_distribution<> dist(l);
      return dist(gent);
    }

}

int boost_poicond1(double l) //Regular Poisson random variate

{
    //cout << l << endl;
    if (l==0.0){
      return 1.0;
    }

    else{
      int res=0;
      while (res<1){
          boost::random::poisson_distribution<> dist(l);
          res = dist(gent);
      }
      return(res);
    }

}


int poicond1(double l) //Poisson random variate conditional on the result being at least 1
{

double expl=exp(-l);
double p=expl+BRand::Controller.nextClosed()*(1-expl);
int k=1;
while (p>=expl)
{k++;
p=p*BRand::Controller.nextClosed();
}
return k-1;

}
