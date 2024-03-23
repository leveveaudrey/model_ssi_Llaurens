#ifndef _SSI
#define _SSI

#include <iostream>
#include <stdlib.h>
#include "MersenneTwister.h"
#include <vector>
#include <fstream>
#include<cmath>
#include <algorithm>
#include <bitset>

#include "randomc.h"                   // define classes for random number generators
#include "stocc.h"                     // define random library classes


using namespace std;
class genotypes;
class user_parameters;

/****************************************************************************/
 /******************     ALLELIC STRUCTURE BY DEME    ************************/
 /****************************************************************************/

 class alleles
 {
    vector<double > A;

    public:
    alleles(){};
    alleles(long na){A.resize(na,0.0);}
    alleles(genotypes);
    alleles(vector<alleles >);
    double& operator[](int x){return A[x];}
    double operator-(alleles);
    ~alleles(){};

    void display(){for (long i=0;i<A.size();i++) cout << A[i]<<'\t'; cout << endl;}
    void file_out(long);
    void diversity_file_out(long);
    void final_file_out();
    
    double diversity();
    long size(){return A.size();}
};



/****************************************************************************/
/*********************** DOMINANCE : the fast way****************************/
/****************************************************************************/

class dominance_parameters{

int dom_type; // 0 : dom-dom ; 1 : dom-cod ; 2 : cod - dom ;
                //  3 : cod-cod ; 999 : neutral (dans le sens pollen-stigma)
long class_allelic_frequency_number; // number of classes for the computation of allelic frequencies distribution
vector<long > dom;
vector <long> ANb;//contains the number of allele by class for which its frequency is different from 0
vector<vector<long> > ADistrib;
vector<long > nclass; // contains the class to which allele i belongs
vector<vector<long> > AFDistrib;//contains the allelic frequency distribution
vector<vector <vector <long > > > GFDistrib;//contains the genotypic frequeny distribution

public:

dominance_parameters(){}
dominance_parameters(const char*);
dominance_parameters(const char*, long);
dominance_parameters(int, vector<long >);
dominance_parameters(int, vector<long >,long);
dominance_parameters(long);
dominance_parameters(long, long);
~dominance_parameters(){}

long& operator[](int x){return dom[x];}
void display();
long get_nclass(long i){if (i>nclass.size()-1) cout << "Ach probleme !"; return nclass[i];}
long get_size(){return dom.size();}
long get_allele_number(){return nclass.size();}
double dom_stigma(long, long);
double dom_pollen(long, long);

long check_allele_number(alleles, double);
void check_allele_distrib();
void check_allelic_frequency_distrib(alleles,long );
void check_genotypic_frequency_distrib(genotypes,long );

void file_out();
void ADistrib_file_out();
void AFDistrib_file_out();
void GFDistrib_file_out();
void ANb_file_out();
};



/****************************************************************************/
/********** DOMINANCE : the fast way, several at a time *********************/
/****************************************************************************/

class dp_table{

vector<dominance_parameters > dp;

public:

dp_table(){}
dp_table(const char*);
dp_table(const char*, long);
dp_table(string);
dp_table(string, long);
~dp_table(){}

dominance_parameters& operator[](int x){return dp[x];}
void file_out(){for (long i=0;i<dp.size();i++) dp[i].file_out();}
long size(){return dp.size();}

    };
    


/****************************************************************************/
/***************** DOMINANCE RELATIONSHIP CALCULATIONS***********************/
/****************************************************************************/

class dominance{
vector <vector <double > > D_pollen, D_stigma; // relation de dominance entre allèles dans pollen et stigmate
vector <vector <vector <double > > > A, P; //

char file_D_pollen, file_D_stigma;
vector <vector <int > > D;
double threshold;
vector <vector <vector < vector <double > > > > p; // takes value 1 if cross ij x kl is compatible else 0

public:

dominance();
dominance(dominance_parameters, double);
dominance(const char*, const char*, double);
~dominance(){};
long get_allele_number(){if (D_pollen.size()==D_stigma.size()) return D_pollen.size(); else return 0;}
vector <vector < vector <double > > >& operator[](int x){return p[x];}

void display();
    };


/****************************************************************************/
/********************     SHELTERED LOAD CLASS    ***************************/
/****************************************************************************/

class sheltered_load_loci{
     
  vector<vector<bitset<100> > > vb; // contains structure of sheltered load loci associated to each S-alleles
  vector<vector<long > > vbf;    // contains the frequencies of each sheltered load allele inside each S-allele.
 long loci_number;
 double sel;// coefficient of selection by locus;
 vector <vector<long> > distrib;
 vector <long> dist_ho; //distribution des homozygotes;
 vector <vector<long> > F; 
 
  public:
      
sheltered_load_loci(long, long, double, alleles, long);
sheltered_load_loci(){};     
~sheltered_load_loci(){}; 

bool selection_ID_multiplicative(sheltered_load_loci, long, long);
bool selection_ID_dominance(sheltered_load_loci, long, long);
void  fixation_number(long, bool);

void clean();
void initialize();

void display(long, long);
void mutation(alleles,double, double,long);
long count();
long freq_count();
long distrib_all();
void ADistrib_file_out()   ;
void ADistrib_file_out(long)   ;
void distribution_homozygotes_ID (long);
      
      };



/****************************************************************************/
/*****************     GENOTYPIC STRUCTURE BY DEME    ***********************/
/****************************************************************************/

class genotypes{
vector<vector <double > > G;
//vector <vector <vector < vector <double > > > > P;//contient les proba conditionnelles de fecondation (dans le cas du modèle de Wright)
bool extinction;// indicates if population goes extinct because of compatible genotypes loss
public:

genotypes(){extinction=false;};
genotypes(long na){G.resize(na); for (long i=0; i<na; i++) G[i].resize(na);extinction=false;}
genotypes(vector<genotypes >, long, double); // permet de fabriquer le nuage pollinique qui arrive par migration sur chacun des demes

vector<double >& operator[](int x){return G[x];}
genotypes(const char*);
~genotypes(){};

void initialisation();
void initialisation(double f){for (long i=0;i<G.size(); i++) for (long j=0;j<G.size();j++) G[i][j]=f;}
void initialisation(const char*);
void initialisation(int, genotypes); // permet l'initialisation des fréquences dans le cas de l'introduction d'un allèle en faible fréquence en se basant sur les fréquences à l'équilibre d'un cas avec un allèle ne moins

void display(){for(long i=0;i<G.size(); i++){ for (long j=0;j<G.size();j++) cout << G[i][j] << '\t'; cout << endl;}}
long size(){return G.size();}

// determine la fréquence des génotypes des graines produites dans ce dème après fécondation en fonction du nuage pollinique  qui arrive
void fecondation_FS(genotypes, dominance);//avec modèle FS;
void fecondation_wFS(genotypes, dominance); // avec modèle type Wright (sans FS);
void fecondation_FS(genotypes, dominance, dominance_parameters);//avec modèle FS, ne sert que dans le cas où on entre les relations de dominance de manière simplifiée;
void fecondation_wFS(genotypes, dominance, dominance_parameters); // avec modèle type Wright (sans FS);
void fecondation_FS(genotypes, dominance, alleles);//avec modèle FS, sert lorsque dérive essentiellement : ne calcule pas la fréquence dun génotype si un des allèles a une fréquence nulle
void fecondation_wFS(genotypes, dominance, alleles); // avec modèle type Wright (sans FS);

void derive(long); // Argument is the number of individuals N (2 N = number of chromosomes in the population)
sheltered_load_loci derive_ID(long, sheltered_load_loci); // Argument is the number of individuals N (2 N = number of chromosomes in the population)
sheltered_load_loci derive_ID_dominance(long, sheltered_load_loci);

void mutation_KAM(double, long);
void mutation_KAM_forcee(int, long, long);

void allele_number_estimation_paxman(long, long);

void file_out(long);
void final_file_out();
void final_file_out(dominance_parameters);

bool _is_extinct(){return extinction;}

void sampling_without_replacement(long, long);
    };



/****************************************************************************/
/***********************     METAPOPULATION    ******************************/
/****************************************************************************/

class metapopulation{

vector<genotypes > mp;
//genotypes pol_tot;
double m, mu; //migration rate and mutation rate
long Ni, Nd; // individuals and demes number

public:

metapopulation(){};
metapopulation(genotypes , double, double, long, long);
~metapopulation(){};

void mutation_KAM(){for (long i=0;i<mp.size();i++) mp[i].mutation_KAM(mu, Ni);}
void derive(){for (long i=0;i<mp.size();i++) mp[i].derive(Ni);}
void fecondation_FS(dominance);
void fecondation_wFS(dominance);

genotypes islands_migration_model(long);
void display(){for (long i=0;i<mp.size();i++) {cout << "POP" << i << endl; mp[i].display();}}
void display_alleles(){for (long i=0;i<mp.size();i++) {cout << "POP" << i << endl; (alleles(mp[i])).display(); cout << endl;}}

double global_Fst();
alleles mean_allelic_frequency();

void file_out(long);
void final_file_out();

    };

/************************************************************************************/
/************************************************************************************/
/***********************************Dialog Class*************************************/
/************************************************************************************/
/************************************************************************************/

class user_parameters{
      
      public:
    
    string parameters_file, //name of the file containing the "simple" cases (dom, domcod, or cod, for any number of alleles by class)
           pollen_dominance, //name of the file containing the dominance relationships in pollen
           pistil_dominance, //name of the file containing the dominance relationships in pistil
           initial_genotypic_frequencies, // name of the file containing the initial genotypic frequencies (used only for frequencies change in a generation)
           user_file;     //name of the file containing all informations of a given repetition;
    long gen_max,
         Nindividuals, //number of diploid individuals in the population
         oubli,        // number of "forget" generation, after which various statistics computations are performed
         repet_max,    // number of independent repetitions
         class_allelic_frequency_number,// number of classes used for the frequencies distribution estimation
         Nsample,                          // sample size (number of diploid individuals sampled used)
         generations_between_sampling,     // number of generations between two samples in a single repetition
         loci_number_ID; //Number of loci involved in sheltered load
    double stop,                           // threshold value defining equilibrum for deterministic equilibrium computations
    threshold_ANb,                         //minimal allelic frequency for an allele to be taken into account in statistics estimation
    mu,                                    //mutation rate
    selection_coefficient,                  //selection coefficient by locus at the sheltered load loci
    mu_ID,                                  //mutation rate by locus by chromosome by generation at the sheltered load loci in the sense 0 to 1
    eta_ID;                                 //mutation rate by locus by chromosome by generation at the sheltered load loci in the sense 1 to 0
    
    bool _derive,               // drift or not drift ?
         _mutation,             // mutation or not ?
         _is_model_FS,          // Selection model through male AND female ways ?
         _is_model_wFS,       // Selection model through male way ONLY ?
         _is_sampling,        // If true, only a sample is used for statistics computations
         _deterministe,         // Deterministic equilibrium computations ?
         _simple,               // Use of simple dominance relationships (dom, domcod or cod) ? 
         _allelic_frequency_distribution, // Computation of allelic frequencies distributions ?
         _genotypic_frequency_distribution,// Computation of genotypic frequencies distributions ?
         _one_generation_change, //Computations of allelic and genotypic frequencies change in a generation ?
         _mutation_ID;
  public:
         user_parameters();
         ~user_parameters(){};
      };

#endif
