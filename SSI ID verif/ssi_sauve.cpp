
#include <iostream>
#include <stdlib.h>
#include "MersenneTwister.h"
#include <vector>
#include <fstream>
#include<cmath>
#include <algorithm>
#include <time.h>
#include "ssi.h"


#include "randomc.h"                   // define classes for random number generators
#include "stocc.h"                     // define random library classes


using namespace std;

extern MTRand al;
extern StochasticLib1 sto;
class genotypes;

long factorial(long n)
{
  long fact = 1;
  for (long i=1; i<=n; fact *= i++);
  return fact;
    }

long combinaison(long k, long n)
{
        if (k==0) return 1; else if (n==k) return 1;
         else return combinaison(k,n-1)+combinaison(k-1,n-1);
    }

double binomial(double p, long k, long n)
{
    return combinaison(k,n)*std::pow((double) p,(double) k)*std::pow((double) (1.0-p), (double)n-k);
    }
    
double paxman_estim(long m, long a)
{
    long it_max=20, it=0;;
    double stop=1e-5;
    double est=(double) a, est_mem=0;
    
    while (it<it_max||abs(est-est_mem)>stop)
    {
        est_mem=est;
        est-=(-a+est-est*std::pow((double) (est-2)/(double)est,(double) m))
                /double(1-std::pow((double) (est-2),(double) m-1)
                *std::pow((double) (est),(double)-m)*(est-2+2*m));
        it++;
        }
    return est;
    }

/****************************************************************************/
 /******************     ALLELIC STRUCTURE BY DEME    ************************/
 /****************************************************************************/

void alleles::file_out(long gen )
{
    ofstream myFile("Alleles.txt",ios::app);
    if (! myFile) cout << "Impossible to open output file Alleles.txt (it must be closed)";
    myFile << gen << '\t';
    for (long i=0;i<A.size();i++) myFile << A[i] << '\t';
    myFile << endl; myFile.close();

    }
    
void alleles::diversity_file_out(long gen )
{
    double ibs=0.0;
    ofstream myFile("Alleles.txt",ios::app);
    if (! myFile) cout << "Impossible to open output file Alleles.txt (it must be closed)";
myFile << gen << '\t';
    for (long i=0;i<A.size();i++) ibs+=A[i]*A[i];
    myFile << 1-ibs << endl; myFile.close();

    }
    
double alleles::diversity()
{
    double ibs=0.0;
    for (long i=0;i<A.size();i++) ibs+=A[i]*A[i];
    return 1-ibs;
}
    
void alleles::final_file_out()
{
    ofstream myFile("F_Alleles.txt",ios::app);
    if (! myFile) cout << "Impossible to open output file F_Alleles.txt (it must be closed)";
    for (long i=0;i<A.size();i++) myFile << A[i] << '\t';
    myFile << endl;
    myFile.close();
    }


double alleles::operator-(alleles a)
{
    double max_diff=0;

for (long i=0;i<A.size();i++) if (max_diff<abs(a.A[i]-A[i])) max_diff=abs(a.A[i]-A[i]);
return max_diff;
    }


alleles::alleles(vector <alleles > va)
{
    A.resize(va[0].size(),0.0);

    for (long i=0;i<va.size();i++)
        {
            for (long j=0;j<=va.size();j++) A[i]+=va[j][i];
            
            A[i]/=va.size();
        }

    }

/****************************************************************************/
/*********************** DOMINANCE : the fast way****************************/
/****************************************************************************/


void dominance_parameters::ANb_file_out()
{
     ofstream myFile("Alleles_Number.txt",ios::app);
     if (! myFile) cout << "Impossible to open output file Alleles.txt (it must be closed)";
     
     for (long i=0;i<ANb.size();i++)
         {
                myFile<< ANb[i]<<'\t';
         }
     myFile<<endl;
     myFile.close();
     }


void dominance_parameters::ADistrib_file_out()
{
     ofstream myFile("Alleles_Number.txt",ios::app);
     if (! myFile) cout << "Impossible to open output file Alleles.txt (it must be closed)";
     
     for (long i=0;i<ADistrib.size();i++)
         {
               for (long j=0;j<ADistrib[i].size();j++) myFile<< ADistrib[i][j]<<'\t';
               myFile<<endl;
         }
     
     myFile.close();
     }
     
void dominance_parameters::AFDistrib_file_out()
{
     ofstream myFile("Allelic_Frequency_Distribution.txt",ios::app);
     if (! myFile) cout << "Impossible to open output file Allelic_Frequency_Distribution.txt (it must be closed)";
     
     double inc= 1.0/(double) (AFDistrib[0].size()-2);
     
     myFile << 0 << '\t';
     for (int i=1;i<=AFDistrib[0].size()-2;i++) myFile<< "]" << (i-1)*inc << ", " << i*inc << "]" << '\t';
     myFile << 1 << endl << endl;
     
     for (long i=0;i<AFDistrib.size();i++)
         {
               for (long j=0;j<AFDistrib[i].size();j++) myFile<< AFDistrib[i][j]<<'\t';
               myFile<<endl;
         }
     myFile<<endl;
     myFile.close();
     }

void dominance_parameters::GFDistrib_file_out()
{
     ofstream myFile("Genotypic_Frequency_Distribution.txt",ios::app);
   if (! myFile) cout << "Impossible to open output file Genotypic_Frequency_Distribution.txt (it must be closed)";
   
   double inc= 1.0/(double) (GFDistrib[0][0].size()-2);
     
     myFile << 0 << '\t';
     for (int i=1;i<=GFDistrib[0][0].size()-2;i++) myFile<< "]" << (i-1)*inc << ", " << i*inc << "]" << '\t';
     myFile << 1 << endl << endl;
   
      for (long i=0;i<GFDistrib.size();i++)
         for (long j=i;j<GFDistrib.size();j++)
                 {
                  for (long k=0;k<GFDistrib[i][j].size();k++) myFile<< GFDistrib[i][j][k]<<'\t';
                  myFile<<endl;
                  }
              
      myFile<<endl;
      myFile.close();
     }


void dominance_parameters::check_allele_distrib()
{
     for (long i=0;i<ANb.size();i++) ADistrib[i][ANb[i]]++;
     }

long dominance_parameters::check_allele_number(alleles a, double t) // t=threshold pour tenir compte d'un allèle
{
            long tot=0;                                           
     for (long i=0;i<ANb.size();i++) ANb[i]=0;
     for (long i=0;i<a.size();i++) if (a[i]>t) ANb[nclass[i]]++;
     for (long i=0;i<ANb.size();i++) tot+=ANb[i];
     return tot;
     }
     
void dominance_parameters::file_out()
{
    ofstream myFile1("Genotypes.txt",ios::app), myFile2("Alleles.txt",ios::app),
                             myFile3("F_Alleles.txt",ios::app),
                             myFile4("F_Genotypes.txt",ios::app),
                             myFile5("Allelic_Frequency_Distribution.txt",ios::app),
                             myFile7("Genotypic_Frequency_Distribution.txt",ios::app),
                             myFile6("Alleles_Number.txt",ios::app);
    if (! myFile1 || !myFile2 || !myFile3|| !myFile4 || !myFile5|| !myFile6) cout << "Impossible to open output files (they must be closed).";

    myFile1<< dom_type << endl; myFile2<< dom_type << endl; myFile3<< dom_type << endl;
        myFile4<< dom_type << endl; myFile5<< dom_type << endl;myFile6<< dom_type << endl;
    for (long i=0;i<dom.size();i++)
        {
            myFile1 << dom[i] << '\t';
            myFile2 << dom[i] << '\t';
            myFile3 << dom[i] << '\t';
            myFile4 << dom[i] << '\t';
            myFile5 << dom[i] << '\t';
            myFile6 << dom[i] << '\t';
            myFile7 << dom[i] << '\t';
        }
    myFile1 << endl<<endl; myFile1.close();
    myFile2 << endl<<endl; myFile2.close();
    myFile3 << endl<<endl; myFile3.close();
    myFile4 << endl<<endl; myFile4.close();
    myFile5 << endl<<endl; myFile5.close();
    myFile6 << endl<<endl; myFile6.close();
    myFile7 << endl<<endl; myFile7.close();
    }

void dominance_parameters::display()
 {
    cout << "Dominance type : " << dom_type << endl << "Allele number by class : ";
    for (long i=0;i<dom.size();i++) cout << dom[i] << '\t';
    cout << endl << endl << "Dominance class by allele : " ;
    for (long i=0;i<nclass.size();i++) cout << nclass[i] << '\t';
}

double dominance_parameters::dom_pollen(long i, long j)
{
 int dec; // takes value 0 if both allele are in different class 1 else
 double res;
 if (nclass[i]==nclass[j]) dec=1; else dec =0;
 
 switch(dom_type) // 0 : dom-dom ; 1 : dom-cod ; 2 : cod - dom ;
                    // 3 : cod-cod  ; 999 : neutral (dans le sens pollen-stigma)
 {
    case 0 : res= dec*0.5 + (1-dec)*1; break;
    case 1 : res=  dec*0.5 + (1-dec)*1; break;
    case 2 : res=  0.5; break;
    case 3 : res=  0.5; break;
    case 999 : res = 0.0; break;
    default:break;
        }
 return res;
    }
    
double dominance_parameters::dom_stigma(long i, long j)
{
int dec; // takes value 0 if both allele are in different class 1 else
 double res;
 if (nclass[i]==nclass[j]) dec=1; else dec =0;

 switch(dom_type) // 0 : dom-dom ; 1 : dom-cod ; 2 : cod - dom ;
                    //  3 : cod-cod ; 999 : neutral (dans le sens pollen-stigma)
 {
    case 0 : res= dec*0.5 + (1-dec)*1; break;
    case 1 : res=  0.5; break;
    case 2 : res=  dec*0.5 + (1-dec)*1; break;
    case 3 : res=  0.5; break;
    case 999 : res = 0.0; break;
    default:break;
        }
 return res;
    }
    
    
dominance_parameters::dominance_parameters(int dt, vector<long > vl)
{
long temp, c, imax=0, i=0;

dom_type=dt;
dom=vl;

temp=0;c=0; imax=dom[0];
for (long k=0;k<dom.size();k++) temp+=dom[k];
nclass.resize(temp);
while(i<temp)
    {
        if(i<imax)
        {
            nclass[i]=c;
            i++;
            }
            else {c++;imax+=dom[c];}
        }
    ANb=dom;
    ADistrib.resize(dom.size());
    for (long j=0;j<ADistrib.size();j++) ADistrib[j].resize(dom[j]+1,0);
    
    }
    
dominance_parameters::dominance_parameters(int dt, vector<long > vl, long nc)
{
long temp, c, imax=0, i=0;
class_allelic_frequency_number=nc;
dom_type=dt;
dom=vl;

temp=0;c=0; imax=dom[0];
for (long k=0;k<dom.size();k++) temp+=dom[k];
nclass.resize(temp);
while(i<temp)
    {
        if(i<imax)
        {
            nclass[i]=c;
            i++;
            }
            else {c++;imax+=dom[c];}
        }
    ANb=dom;
    ADistrib.resize(dom.size());
    AFDistrib.resize(temp); 
    for (long j=0;j<temp;j++) AFDistrib[j].resize(class_allelic_frequency_number);
    for (long j=0;j<ADistrib.size();j++) ADistrib[j].resize(dom[j]+1,0);
    
    }
    
dominance_parameters::dominance_parameters(const char *fp)
{
char ch;
 long temp, c, imax=0, i=0;
 ifstream fparam(fp);
 
 if (! fparam){cout << "Error opening parameter file (for simple cases computation)" << endl; cin >> ch; exit(1);}
 
 fparam >> dom_type;

 while( fparam >> temp )  dom.push_back(temp);
fparam.close();

temp=0;c=0; imax=dom[0];
for (long k=0;k<dom.size();k++) temp+=dom[k];
nclass.resize(temp);
while(i<temp)
    {
        if(i<imax)
        {
            nclass[i]=c;
            i++;
            }
            else {c++;imax+=dom[c];}
        }
        
    ANb=dom;
    ADistrib.resize(dom.size());
    for (long i=0;i<ADistrib.size();i++) ADistrib[i].resize(dom[i]+1,0);
}

dominance_parameters::dominance_parameters(const char *fp, long nc)
{
char ch;
 long temp, c, imax=0, i=0;
 class_allelic_frequency_number=nc;
 ifstream fparam(fp);
  if (! fparam){cout << "Error opening parameter file (for simple cases computation)" << endl; cin >> ch;exit(1);}
 
 fparam >> dom_type;

 while( fparam >> temp )  dom.push_back(temp);
fparam.close();

temp=0;c=0; imax=dom[0];
for (long k=0;k<dom.size();k++) temp+=dom[k];
nclass.resize(temp);
while(i<temp)
    {
        if(i<imax)
        {
            nclass[i]=c;
            i++;
            }
            else {c++;imax+=dom[c];}
        }
        
    ANb=dom;
    ADistrib.resize(dom.size());
    for (long i=0;i<ADistrib.size();i++) ADistrib[i].resize(dom[i]+1,0);
}

dominance_parameters::dominance_parameters(long nb)
{
 ANb.resize(nb);
 ADistrib.resize(nb);
 for (long i=0;i<nb;i++) {nclass.push_back(i);ADistrib[i].resize(2);}                                               
                                                }
                                                
dominance_parameters::dominance_parameters(long nb, long nc)
{
class_allelic_frequency_number=nc;
 ANb.resize(nb);
 ADistrib.resize(nb);
 AFDistrib.resize(nb); 
 for (long i=0;i<nb;i++) {nclass.push_back(i);ADistrib[i].resize(2);AFDistrib[i].resize(nc);}
}                                                

void dominance_parameters::check_allelic_frequency_distrib(alleles a, long t)
{
     long c;
     double wd=1.0/(double) t;
     
     if (AFDistrib.empty()) {AFDistrib.resize(a.size()); for (int i=0;i<a.size();i++) AFDistrib[i].resize(t+2);}
     
     
     for (long i=0;i<a.size();i++)
         {
          c=0;
          while(a[i]>c*wd) c++;
          AFDistrib[i][c]++;
          }
     }
     
void dominance_parameters::check_genotypic_frequency_distrib(genotypes g, long t)
{
     long c;
     double wd=1.0/(double) t;
     
     if (GFDistrib.empty()) 
        {
         GFDistrib.resize(g.size()); 
             for (int i=0;i<g.size();i++) 
                 {
                      GFDistrib[i].resize(g.size());
                      for (int j=0;j<g.size();j++)
                         GFDistrib[i][j].resize(t+2);
                      
                  }
         }
     
      for (int i=0;i<g.size();i++)                  
          for (int j=i;j<g.size();j++)
         {
          c=0;
          while(g[i][j]>c*wd) c++;
          GFDistrib[i][j][c]++;
          if (i!=j) GFDistrib[j][i][c]++;
          }
     }

/****************************************************************************/
/********** DOMINANCE : the fast way, several at a time *********************/
/****************************************************************************/
   
dp_table::dp_table(const char *fp)
{
    ifstream fparam(fp);
    vector<long > vt;
    long dom_type, temp;
    bool test=true;
    char ch;
    
   if (! fparam){cout << "Error opening parameter file (for simple cases computation)" << endl; cin >> ch;exit(1);}
            
while( fparam >> temp) {
 if (test) {test=false; dom_type=temp;}
        else { vt.push_back(temp);
         if (fparam.peek()=='\n') {
            test=true;
            dp.push_back(dominance_parameters(dom_type,vt));
            vt.resize(0);
            }}
            }
dp.push_back(dominance_parameters(dom_type,vt));
fparam.close();
}

dp_table::dp_table(const char *fp, long nc)
{
    ifstream fparam(fp);
    vector<long > vt;
    long dom_type, temp;
    bool test=true;
    char ch;
    
   if (! fparam){cout << "Error opening parameter file (for simple cases computation)" << endl; cin >> ch;exit(1);}
            
while( fparam >> temp) {
 if (test) {test=false; dom_type=temp;}
        else { vt.push_back(temp);
         if (fparam.peek()=='\n') {
            test=true;
            dp.push_back(dominance_parameters(dom_type,vt, nc));
            vt.resize(0);
            }}
            }
dp.push_back(dominance_parameters(dom_type,vt));
fparam.close();
}

dp_table::dp_table(string fp)
{

    ifstream fparam(fp.c_str());
    vector<long > vt;
    long dom_type, temp;
    bool test=true;
    char ch;
    
    
    
   if (! fparam){cout << "Error opening parameter file (for simple cases computation)" << endl; cin >> ch;exit(1);}
        
while( fparam >> temp) { 
 if (test) {test=false; dom_type=temp;   }
        else { vt.push_back(temp);
         if (fparam.peek()=='\n') {
            test=true;
            dp.push_back(dominance_parameters(dom_type,vt));
            vt.resize(0);
            }}
            }            
dp.push_back(dominance_parameters(dom_type,vt));


fparam.close();
}

dp_table::dp_table(string fp, long nc)
{

    ifstream fparam(fp.c_str());
    vector<long > vt;
    long dom_type, temp;
    bool test=true;
    char ch;
    
    
    
  if (! fparam){cout << "Error opening parameter file (for simple cases computation)" << endl; cin >> ch;exit(1);}
        
while( fparam >> temp) { 
 if (test) {test=false; dom_type=temp;   }
        else { vt.push_back(temp);
         if (fparam.peek()=='\n') {
            test=true;
            dp.push_back(dominance_parameters(dom_type,vt, nc));
            vt.resize(0);
            }}
            }            
dp.push_back(dominance_parameters(dom_type,vt));


fparam.close();
}

/****************************************************************************/
/***************** DOMINANCE RELATIONSHIP CALCULATIONS***********************/
/****************************************************************************/

void dominance::display()
{
    cout << " Dpollen" << endl;
    for (long i=0;i<D_pollen.size();i++)
        {
            for (long j=0;j<D_pollen.size();j++)  cout << D_pollen[i][j] << '\t';
            cout << endl;
        }
    
    cout << endl << endl << " Dstigma" << endl;
     for (long i=0;i<D_stigma.size();i++)
        {
            for (long j=0;j<D_stigma.size();j++)  cout << D_stigma[i][j] << '\t';
            cout << endl;
        }
    }
    
dominance::dominance()
{
D_pollen, D_stigma;D;
threshold=0.0;
    }
    
dominance::dominance(const char *fp, const char *fs, double t)
{
char ch;
threshold=t;
 double temp, tems;
 double calc;
 long nb_alleles =0;
 int i=0, j=0;
 ifstream fparam(fp), finp(fp), fins(fs);
 
 if (! fparam || !finp || !fins){cout << "Error opening parameter file (for specific dominance relationships computation)" << endl; /*cin >> ch;*/abort();}
 
 cin >> ch;
 while( fparam >> temp )  nb_alleles++;
 nb_alleles = (long) sqrt((double) nb_alleles);
fparam.close();
 
 D_pollen.resize(nb_alleles); D_stigma.resize(nb_alleles);
 A.resize(nb_alleles);  P.resize(nb_alleles);
 p.resize(nb_alleles);
 
 for (long k=0;k<nb_alleles;k++)
    {
        D_pollen[k].resize(nb_alleles);
        D_stigma[k].resize(nb_alleles);
        A[k].resize(nb_alleles);
        P[k].resize(nb_alleles);
        p[k].resize(nb_alleles);
        for (long l=0;l<nb_alleles;l++)
            {
                A[k][l].resize(nb_alleles);
                P[k][l].resize(nb_alleles);
                p[k][l].resize(nb_alleles);
                 for (long m=0;m<nb_alleles;m++)
                    p[k][l][m].resize(nb_alleles);
                }
     }
     
while( finp >> temp )
    {
        fins >> tems;
        if (j==nb_alleles) {j=0;i++;}
        D_pollen[i][j]=(double) temp;D_stigma[i][j]=(double) tems;
        j++;
    }
finp.close(); fins.close();

for (long i=0;i<nb_alleles;i++)
    for(long j=0;j<nb_alleles;j++)
        for(long k=0;k<nb_alleles;k++)
    {
        if (k==i) {A[i][j][k]= D_pollen[i][j]; P[i][j][k]= D_stigma[i][j]; }
            else if (k==j){A[i][j][k]= D_pollen[j][i]; P[i][j][k]= D_stigma[j][i]; }
                else {A[i][j][k]=P[i][j][k]= 0; }
    }
for (long i=0;i<nb_alleles;i++)
    for(long j=0;j<nb_alleles;j++)
        for(long k=0;k<nb_alleles;k++)
            for(long l=0;l<nb_alleles;l++)
                {   calc=0;
                    for (long m=0;m<nb_alleles;m++)
                        calc+=A[i][j][m]*P[k][l][m];
                    if (calc<= threshold) p[i][j][k][l]=1; else p[i][j][k][l]=0;
                    }
    }
    
dominance::dominance(dominance_parameters d, double t)
{
    threshold=t;
    long nb_alleles=0;
    double calc;
    char ch;
 for (long k=0;k<d.get_size();k++) nb_alleles+=d[k];
 
  D_pollen.resize(nb_alleles); D_stigma.resize(nb_alleles);
 A.resize(nb_alleles);  P.resize(nb_alleles);
 p.resize(nb_alleles);

 for (long k=0;k<nb_alleles;k++)
    {
        D_pollen[k].resize(nb_alleles);
        D_stigma[k].resize(nb_alleles);
        A[k].resize(nb_alleles);
        P[k].resize(nb_alleles);
        p[k].resize(nb_alleles);
        for (long l=0;l<nb_alleles;l++)
            {
                A[k][l].resize(nb_alleles);
                P[k][l].resize(nb_alleles);
                p[k][l].resize(nb_alleles);
                 for (long m=0;m<nb_alleles;m++)
                    p[k][l][m].resize(nb_alleles);
                }
     }

for (long i=0;i<nb_alleles;i++)
    for(long j=i;j<nb_alleles;j++)
    {
        if (i==j) D_pollen[i][i]=D_stigma[i][i]=1;
            else {
                    D_pollen[i][j]= d.dom_pollen(i,j);
                    D_stigma[i][j]=d.dom_stigma(i,j);
                    D_pollen[j][i]=1-D_pollen[i][j];
                    D_stigma[j][i]=1-D_stigma[i][j];
                }
        }


for (long i=0;i<nb_alleles;i++)
    for(long j=0;j<nb_alleles;j++)
        for(long k=0;k<nb_alleles;k++)
    {
        if (k==i) {A[i][j][k]= D_pollen[i][j]; P[i][j][k]= D_stigma[i][j]; }
            else if (k==j){A[i][j][k]= D_pollen[j][i]; P[i][j][k]= D_stigma[j][i]; }
                else {A[i][j][k]=P[i][j][k]= 0; }
    }
for (long i=0;i<nb_alleles;i++)
    for(long j=0;j<nb_alleles;j++)
        for(long k=0;k<nb_alleles;k++)
            for(long l=0;l<nb_alleles;l++)
                {   calc=0;
                    for (long m=0;m<nb_alleles;m++)
                        calc+=A[i][j][m]*P[k][l][m];
                    if (calc<= threshold) p[i][j][k][l]=1; else p[i][j][k][l]=0;
                    }
    }


    
/****************************************************************************/
/********************     SHELTERED LOAD CLASS    ***************************/
/****************************************************************************/
    
sheltered_load_loci::sheltered_load_loci(long n, long nS, double s, alleles a, long N)
{
  loci_number= n;
  sel=s;
  
  vb.resize(nS);
  vbf.resize(nS);
  
  for (long i=0;i<nS;i++) 
  {
      vb[i].resize(1);                                            
      vbf[i].resize(1,lround(2*N*a[i]));  
  }
}
bool sheltered_load_loci::selection_ID_multiplicative(sheltered_load_loci sll_mem, long a1, long a2)
{
     bitset<100> b;
     long inc1=0, inc2=0;
     long r1,r2;
     
     r1=al.randInt(sll_mem.vbf[a1].back()-1)+1;
     r2=al.randInt(sll_mem.vbf[a2].back()-1)+1;
     
     while(r1>sll_mem.vbf[a1][inc1])inc1++;
     while(r2>sll_mem.vbf[a2][inc2])inc2++;
     
     b=sll_mem.vb[a1][inc1];
     b&=sll_mem.vb[a2][inc2];
     
     if (al() < pow((double) (1.0-sel),(int) b.count()) )
        {
             vbf[a1][inc1]++;
             vbf[a2][inc2]++;             
             return true;
             }     
     else return false;
     }
 
 void sheltered_load_loci::initialize()
 {
 for (int i=0;i<vbf.size();i++) 
     for (int j=0;j<vbf[i].size();j++)
      vbf[i][j]=0;
     }
      
 void sheltered_load_loci::clean()
 {
 vector<vector<bitset<100> > > vbtemp; 
  vector<vector<long > > vbftemp; 
   
   vbtemp.resize(vb.size());
   vbftemp.resize(vbf.size());
   
   for (int i=0;i<vbf.size();i++) 
     for (int j=0;j<vbf[i].size();j++)
         
         if(vbf[i][j]!=0)
         {
          vbftemp[i].push_back(vbf[i][j]);  
          vbtemp[i].push_back(vb[i][j]);             
                         }
   for (int i=0;i<vbf.size();i++) 
     for (int j=1;j<vbftemp[i].size();j++)
         vbftemp[i][j]+=vbftemp[i][j-1];
   
vb=vbtemp;
 vbf=vbftemp;    

      } 
      
      
 void sheltered_load_loci::display()
 {
 
 for (int i=0; i<vbf.size();i++) 
 {
 for (int j=0;j<vbf[i].size();j++){for (int kk=0;kk<loci_number;kk++) cout << vb[i][j][kk]; 
      cout<< '\t' <<vbf[i][j] << '\t'; cout << endl;}
 cout << endl<< endl;
 }
}

void sheltered_load_loci::mutation(alleles a, double mu,long N )
 { 
int i=0,j=0,k,nb_mut,l,ii;  
long inc=0,r;
string ch;
double rand_temp;

bitset<100> btemp;

if (2*N*loci_number*mu>=1.0) nb_mut= lround(2*N*loci_number*mu);
    else nb_mut=sto.Poisson(2*N*loci_number*mu);
                                       //cout<<nb_mut<<endl;  
alleles atemp=a;
atemp[1]=atemp[0]+a[1];
atemp[2]=atemp[1]+a[2];

for (l=0;l<nb_mut;l++) 
 { 
 i=0;
 j=0;
 rand_temp=al();          
     
 while (rand_temp>atemp[i])i++; 
                                            //cout<<"i : "<<i<<endl;
 r=al.randInt(vbf[i].back()-1)+1;
                                        //cout<<"r : "<<r<<endl;
 while (r>vbf[i][j])j++; 
                                        //cout<<"j : "<<j<<endl;
                                        
                                        //cout << vb[i][j] << endl ;
 k=al.randInt(loci_number-1);           //cout << vb[i][j] << endl;
  
btemp=vb[i][j];
btemp[k]= !btemp[k];  

ii=0;

while(btemp!=vb[i][ii] && ii<vb[i].size())ii++;

if (ii<vb[i].size()) for (int kk=j;kk<=ii-1;kk++) vbf[i][kk]--;
  
     else {for (int kk=j;kk<vbf[i].size();kk++) vbf[i][kk]--; vbf[i].push_back(vbf[i].back()+1); vb[i].push_back(btemp);}                        
 }
}
       
/****************************************************************************/
/*****************     GENOTYPIC STRUCTURE BY DEME    ***********************/
/****************************************************************************/

genotypes::genotypes(const char *fp)
{
    extinction=false;                           

    ifstream fparam(fp);
    
    vector<double > vt;
    
    double temp;
    bool test=false;
    char ch;
    
   if (! fparam){cout << "Error opening parameter file (for initial genotypic frequencies)" << endl; cin >> ch;exit(1);}
            
while( fparam >> temp) {
         vt.push_back(temp);
         if (fparam.peek()=='\n' ||fparam.eof() ) 
            {
            G.push_back(vt);
            vt.resize(0);
            }
         }
fparam.close();

for (long i=0;i<G.size();i++) if (G.size()!=G[i].size()) test=true;

if (test)
   {
   cout << "Error in Initial Genotypic Frequencies Data : missing data." ;
   cin >> ch;
    }
for (long i=0;i<G.size();i++) for (long j=i+1;j<G.size();j++) G[j][i]=G[i][j];

}

void genotypes::allele_number_estimation_paxman(long it, long sample_size)
{
    vector<double > estim, sorted_estim;
    vector<vector <double > > Gtemp;
    vector<bool > sampled_alleles;
    vector<long >distrib, sorted_distrib;
    long tot, jit, m, n;
    double r;
    double estn, mean=0, var=0;
    char ch;

Gtemp.resize(G.size());
sampled_alleles.resize(G.size());
estim.resize(1,0.0); distrib.resize(1,0);

 for(long i=0;i<Gtemp.size();i++)
    {
        Gtemp[i].resize(G.size());
        if (i==0) Gtemp[0][0]=G[0][0];
        for (long j=i;j<G.size();j++)
            {
                if (j!=i && j!=0) Gtemp[i][j]=Gtemp[i][j-1]+G[i][j];
                else if (j!=0) Gtemp[i][j]=Gtemp[i-1][G.size()-1]+G[i][j];
            }
        }
    
for (long i=0;i<it;i++)
{
    tot=0;
    for (long j=0;j<sampled_alleles.size();j++) sampled_alleles[j]=false;
    for(long ni=0;ni<sample_size;ni++)
            {
                r=al();
                m=n=0;
                while (r>Gtemp[m][n]) {if (n<Gtemp.size()-1) n++; else {m++;n=m;}}
               if (!sampled_alleles[m]) {sampled_alleles[m]=true; tot++;}
               if (!sampled_alleles[n]) {sampled_alleles[n]=true; tot++;}
                }
        estn=paxman_estim(sample_size,tot);
        jit=0;
        if (estim[0]==0.0) estim[0]=estn;

        while(estim[jit]!=estn && jit<estim.size()-1)jit++;
        if (estim[jit]!=estn) {estim.push_back(estn);distrib.push_back(1);}
            else distrib[jit]++;
}

sorted_estim=estim;
sort(sorted_estim.begin(), sorted_estim.end());
sorted_distrib.resize(distrib.size());

for (long i=0;i<distrib.size();i++)
{
    jit=0;
    while(sorted_estim[i]!=estim[jit]) jit++;
    sorted_distrib[i]=distrib[jit];
    }
    
 for (long i=0;i<distrib.size();i++)
    {
        mean+=estim[i]*(double) distrib[i]/(double) it;
        var+=estim[i]*estim[i]*(double) distrib[i]/(double) it;
        }
cout << endl;
for (long i=0;i<distrib.size();i++) {cout << sorted_distrib[i] << '\t' << sorted_estim[i]<< endl;}
cout << "Mean : " << mean << '\t' << "Variance : " << var-mean*mean << endl;
cout << endl;
}

genotypes :: genotypes(vector<genotypes > v, long d_focal, double mig)
{

    
   G.resize(v[0].size());
    for (long i=0;i<v[0].size();i++) G[i].resize(v[0].size(),0.0);
    
        for (long i=0; i<G.size(); i++)
            for (long j=i;j<G.size(); j++)
                {
                    for (long k=0;k<v.size();k++)
                        if (k!= d_focal) G[i][j]+=mig*v[k][i][j];

                        G[i][j]/(double) (v.size()-1);
                        G[i][j]+=(1.0-mig)*v[d_focal][i][j];

                    G[j][i]=G[i][j];
                }
                
    extinction=false;
}

void genotypes::file_out(long gen)
{
    ofstream myFile1("Genotypes.txt",ios::app);

     if (! myFile1) cout << "Impossible to open output file Genotypes.txt (it must be closed)";
myFile1 << gen << '\t';

    for (long i=0;i<G.size();i++)
        for (long j=i;j<G.size();j++) myFile1 << G[i][j] << '\t';
    myFile1 << endl; myFile1.close();
    
    }

void genotypes::final_file_out()
{
    ofstream myFile1("F_Genotypes.txt",ios::app);

    if (! myFile1) cout << "Impossible to open output file F_Genotypes.txt (it must be closed)";

    for (long i=0;i<G.size();i++){
        for (long j=0;j<G.size();j++) myFile1 << G[i][j] << '\t'; myFile1  << endl;}
    myFile1 << endl; myFile1.close();
    
    }

void genotypes::final_file_out(dominance_parameters dp)
{
    ofstream myFile1("F_Genotypes.txt",ios::app), myFile2("Heterozygotes.txt",ios::app);

    if (! myFile1 || !myFile2) cout << "Impossible to open output file F_Genotypes.txt and/or Heterozygotes.txt (they must be closed)";

    for (long i=0;i<G.size();i++)
        {
        for (long j=0;j<G.size();j++) myFile1 << G[i][j] << '\t';
           myFile1 << endl; }
           
/*vector<vector <double > >  Ghet;

Ghet.resize(dp.get_size());
for (long i=0;i<dp.get_size();i++) Ghet[i].resize(dp.get_size());
       
           for (long m=0;m<G.size();m++)
                for (long n=m+1;n<G.size();n++)
                     Ghet[dp.get_nclass(m)][dp.get_nclass(n)]=G[m][n];
           
    for (long i=0;i<Ghet.size();i++)
        {
        for (long j=0;j<Ghet.size();j++) myFile2 << Ghet[i][j] << '\t';
           myFile2 << endl; }
     */      
    myFile1 << endl; myFile1.close();
     myFile2 << endl; myFile2.close();

    }



void genotypes::mutation_KAM(double mu, long N)
{
long nb_mut,m,n,x;
double r=al();
vector<vector <double > > Gtemp;

if (2*N*mu>=1.0) nb_mut= lround(2*N*mu);
    else nb_mut=sto.Poisson(2*N*mu);


if (nb_mut!=0)
{
    Gtemp.resize(G.size());

   for(long i=0;i<Gtemp.size();i++)
    {
        Gtemp[i].resize(G.size());
        if (i==0) Gtemp[0][0]=G[0][0];
        for (long j=i;j<G.size();j++)
            {
                if (j!=i && j!=0) Gtemp[i][j]=Gtemp[i][j-1]+G[i][j];
                else if (j!=0) Gtemp[i][j]=Gtemp[i-1][G.size()-1]+G[i][j];
            }
        }
}
    
    for(long i=0;i<nb_mut;i++)
            {
                        r=al();
                        m=n=0;
                        while (r>Gtemp[m][n])
                             {if (n<Gtemp.size()-1) n++; else {m++;n=m;}}
                        
                G[m][n]-=1/(double) N;
                
                G[n][m]=G[m][n];
                x=n;n=al.randInt(G.size()-1);
                
                while(n==x)n=al.randInt(G.size()-1);
                G[m][n]+=1/(double) N;
                G[n][m]=G[m][n];

               }
    }

void genotypes::mutation_KAM_forcee(int rare_allele, long N, long nmut)
{
long m,n,x;
double r;
vector<vector <double > > Gtemp;



    Gtemp.resize(G.size());

   for(long i=0;i<Gtemp.size();i++)
    {
        Gtemp[i].resize(G.size());
        if (i==0) Gtemp[0][0]=G[0][0];
        for (long j=i;j<G.size();j++)
            {
                if (j!=i && j!=0) Gtemp[i][j]=Gtemp[i][j-1]+G[i][j];
                else if (j!=0) Gtemp[i][j]=Gtemp[i-1][G.size()-1]+G[i][j];
            }
        }
                for (int k=0;k<nmut;k++)
                {
                        
                        r=al();
                        m=n=0;
                        while (r>Gtemp[m][n])
                             {if (n<Gtemp.size()-1) n++; else {m++;n=m;}}
                 if (G[m][n]>=1/(double) N && G[m][rare_allele]<1-1/(double) N)
                    {       
                            G[m][n]-=1/(double) N;
                            G[n][m]=G[m][n];
                            if (al()<0.5) m=n;
                            G[m][rare_allele]+=1/(double) N;
                            G[rare_allele][m]=G[m][rare_allele];
                    }
                 else k--; 
                }
    }

void genotypes::derive(long N)
{
    vector<vector <double > > Gtemp;
    double r;
    long m,n;
    string ch;
    Gtemp.resize(G.size());
   for(long i=0;i<Gtemp.size();i++)
    {
        Gtemp[i].resize(G.size());
        if (i==0) Gtemp[0][0]=G[0][0];
        for (long j=i;j<G.size();j++)
            {
                if (j!=i && j!=0) Gtemp[i][j]=Gtemp[i][j-1]+G[i][j];
                
                else if (j!=0) Gtemp[i][j]=Gtemp[i-1][G.size()-1]+G[i][j];
                G[i][j]=G[j][i]=0.0;
            }
        }
    if (fabs((Gtemp.back()).back()-1.0)<10e-10) (Gtemp.back()).back()=1.0; else {cout << "Rounding problem"; cin >> ch;}
  
    for(long ni=0;ni<N;ni++)
            {
                r=al();
                m=n=0;
                while (r>Gtemp[m][n]) {if (n<Gtemp.size()-1) n++; else {m++;n=m;}if(m>=Gtemp.size()) {cout << m << '\t'<< r << '\t' << (Gtemp.back()).back();cin >> ch;}}
                 
                G[m][n]+=1/(double) (N);
                G[n][m]=G[m][n];
                }
}

sheltered_load_loci genotypes::derive_ID(long N, sheltered_load_loci sll)
{
      vector<vector <double > > Gtemp;
    double r, survival;
    long m,n, ni=0;
    string ch;
    
    sheltered_load_loci sll_mem=sll;
    
    sll.initialize();
    
    
    Gtemp.resize(G.size());
   for(long i=0;i<Gtemp.size();i++)
    {
        Gtemp[i].resize(G.size());
        if (i==0) Gtemp[0][0]=G[0][0];
        for (long j=i;j<G.size();j++)
            {
                if (j!=i && j!=0) Gtemp[i][j]=Gtemp[i][j-1]+G[i][j];
                
                else if (j!=0) Gtemp[i][j]=Gtemp[i-1][G.size()-1]+G[i][j];
                G[i][j]=G[j][i]=0.0;
            }
        }
    if (fabs((Gtemp.back()).back()-1.0)<10e-10) (Gtemp.back()).back()=1.0; else {cout << "Rounding problem"; cin >> ch;}
  
    while(ni<N)
            {
                r=al();
                m=n=0;
                while (r>Gtemp[m][n]) {if (n<Gtemp.size()-1) n++; else {m++;n=m;}if(m>=Gtemp.size()) {cout << m << '\t'<< r << '\t' << (Gtemp.back()).back();cin >> ch;}}
                    
                if (sll.selection_ID_multiplicative(sll_mem,m,n)) {G[m][n]+=1/(double) (N); G[n][m]=G[m][n];ni++;}
                }
    
    /*cout << "Avant clean :" << endl;
    sll.display();
    //cin >> ch;*/
                
   sll.clean(); 
   return sll;            
}

void genotypes::initialisation(int k, genotypes g)
{
     G.clear();
     G.resize(g.size()+1); 
     for (long i=0;i<G.size();i++) G[i].resize(g.size()+1,0.0);
          
     for (long i=0;i<G.size();i++)
         for (long j=0;j<G.size();j++)
         {
            if (i<k && j<k) G[i][j]=g[i][j];
               else if (i<k && j>k) G[i][j]=g[i][j-1];
               else if (i>k && j<k) G[i][j]=g[i-1][j];
               else if (i>k && j>k) G[i][j]=g[i-1][j-1];
             
         }     
     }

void genotypes::initialisation()
{
    double f=1/(double) (G.size()*(G.size()+1)/ (double) 2);
    for (long i=0;i<G.size(); i++)
        for (long j=0;j<G.size();j++)
            G[i][j]=f;
    }
    
void genotypes::initialisation(const char *file)
{
ifstream f(file);
double temp;
long i=0,j=0;

while( f >> temp )
    {
        if (j==G.size()) {j=0;i++;}
        G[i][j]=temp;
        j++;
    }
}

void genotypes::fecondation_FS(genotypes pollen, dominance dom)
{
vector<vector <double > > Gtemp=G;
double tot=0;

for (long m=0;m<G.size();m++)
    for (long n=m;n<G.size();n++)
        {
            G[m][n]=0.0;
            if (m==n)
                for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=m && k!= j&&j!=m)
                                G[m][m]+=0.125*(dom[m][j][m][k]
                                *pollen[m][j]*Gtemp[m][k]+dom[m][k][m][j]
                                *pollen[m][k]*Gtemp[m][j]);

                        G[m][m]+=0.5*(dom[m][m][m][j]*pollen[m][m]*Gtemp[m][j]
                                +dom[m][j][m][m]*pollen[m][j]*Gtemp[m][m]);
                        if(j!=m) G[m][m]+=0.25*dom[m][j][m][j]*pollen[m][j]*Gtemp[m][j];
                    }
            else
                {
                    for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=n && j!= m)
                                G[m][n]+=0.25*(dom[m][j][n][k]
                                *pollen[m][j]*Gtemp[n][k]+dom[n][k][m][j]
                                *pollen[n][k]*Gtemp[m][j]);

                        G[m][n]+=0.5*(dom[m][m][n][j]*pollen[m][m]*Gtemp[n][j]
                                +dom[n][j][m][m]*pollen[n][j]*Gtemp[m][m])
                                +0.5*(dom[m][j][n][n]*pollen[m][j]*Gtemp[n][n]
                                +dom[n][n][m][j]*pollen[n][n]*Gtemp[m][j]);
                    }
                }
               tot+=G[m][n];
            }

   for (long m=0;m<G.size();m++)  for (long n=m;n<G.size();n++) G[m][n]/=(double) tot;

   for (long m=0;m<G.size();m++) for (long n=0;n<m;n++) G[m][n]=G[n][m];
}

 
void genotypes::fecondation_FS(genotypes pollen, dominance dom, dominance_parameters dp)
{
vector<vector <double > > Gtemp=G, EQhet;//EQhet contains genotypic frequencies at equilibrium for alleles in class {i,j} for i!=j (heterozygotes)
vector <double >EQhom;//EQhom contains genotypic frequencies at equilibrium for homozygotes
double tot=0;

EQhet.resize(dp.get_size()); EQhom.resize(dp.get_size());

for (long i=0;i<dp.get_size();i++)
    {EQhet[i].resize(dp.get_size());EQhom[i]=-9.0;
       for (long j=0;j<dp.get_size();j++) EQhet[i][j]=-9.0; }

for (long m=0;m<G.size();m++)
    for (long n=m;n<G.size();n++)
        {
            G[m][n]=0.0;
            if (m==n)
                {
                    if (EQhom[dp.get_nclass(m)]==-9.0)
                    {
                        for (long j=0;j<G.size();j++)
                         {
                                for (long k=0;k<G.size();k++)
                            if (k!=m && k!= j&&j!=m)
                                      G[m][m]+=0.125*(dom[m][j][m][k]
                                    *pollen[m][j]*Gtemp[m][k]+dom[m][k][m][j]
                                     *pollen[m][k]*Gtemp[m][j]);

                               G[m][m]+=0.5*(dom[m][m][m][j]*pollen[m][m]*Gtemp[m][j]
                                       +dom[m][j][m][m]*pollen[m][j]*Gtemp[m][m]);
                                if(j!=m) G[m][m]+=0.25*dom[m][j][m][j]*pollen[m][j]*Gtemp[m][j];
                         }
                       EQhom[dp.get_nclass(m)]=G[m][m];
                    }
                    else G[m][m]=EQhom[dp.get_nclass(m)];
                }
            else
               if (EQhet[dp.get_nclass(m)][dp.get_nclass(n)]==-9.0)
                {
                    for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=n && j!= m)
                                G[m][n]+=0.25*(dom[m][j][n][k]
                                *pollen[m][j]*Gtemp[n][k]+dom[n][k][m][j]
                                *pollen[n][k]*Gtemp[m][j]);

                        G[m][n]+=0.5*(dom[m][m][n][j]*pollen[m][m]*Gtemp[n][j]
                                +dom[n][j][m][m]*pollen[n][j]*Gtemp[m][m])
                                +0.5*(dom[m][j][n][n]*pollen[m][j]*Gtemp[n][n]
                                +dom[n][n][m][j]*pollen[n][n]*Gtemp[m][j]);
                    }
                    EQhet[dp.get_nclass(m)][dp.get_nclass(n)]=G[m][n];
                }
                else G[m][n]=EQhet[dp.get_nclass(m)][dp.get_nclass(n)];
               tot+=G[m][n];
            }

   for (long m=0;m<G.size();m++)  for (long n=m;n<G.size();n++) G[m][n]/=(double) tot;

   for (long m=0;m<G.size();m++) for (long n=0;n<m;n++) G[m][n]=G[n][m];
}
 
 
 void genotypes::fecondation_FS(genotypes pollen, dominance dom, alleles a)
{
vector<vector <double > > Gtemp=G;
double tot=0;

for (long m=0;m<G.size();m++)
    for (long n=m;n<G.size();n++)
        {
            G[m][n]=0.0;
            if (a[m]!=0.0 && a[n]!=0.0)
            {
            if (m==n)
                for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=m && k!= j&&j!=m)
                                G[m][m]+=0.125*(dom[m][j][m][k]
                                *pollen[m][j]*Gtemp[m][k]+dom[m][k][m][j]
                                *pollen[m][k]*Gtemp[m][j]);

                        G[m][m]+=0.5*(dom[m][m][m][j]*pollen[m][m]*Gtemp[m][j]
                                +dom[m][j][m][m]*pollen[m][j]*Gtemp[m][m]);
                        if(j!=m) G[m][m]+=0.25*dom[m][j][m][j]*pollen[m][j]*Gtemp[m][j];
                    }
            else
                {
                    for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=n && j!= m)
                                G[m][n]+=0.25*(dom[m][j][n][k]
                                *pollen[m][j]*Gtemp[n][k]+dom[n][k][m][j]
                                *pollen[n][k]*Gtemp[m][j]);

                        G[m][n]+=0.5*(dom[m][m][n][j]*pollen[m][m]*Gtemp[n][j]
                                +dom[n][j][m][m]*pollen[n][j]*Gtemp[m][m])
                                +0.5*(dom[m][j][n][n]*pollen[m][j]*Gtemp[n][n]
                                +dom[n][n][m][j]*pollen[n][n]*Gtemp[m][j]);
                    }
                }
               tot+=G[m][n];
            }
}
if (tot!=0)   for (long m=0;m<G.size();m++)  for (long n=m;n<G.size();n++) { G[m][n]/=(double) tot;}
   else extinction=true;
   for (long m=0;m<G.size();m++) for (long n=0;n<m;n++) G[m][n]=G[n][m];
}



 
void genotypes::fecondation_wFS(genotypes pollen, dominance dom)
  {
    vector <vector <vector < vector <double > > > > P;
    vector<vector <double > > Gtemp=G;

            P.resize(G.size());
            for (long i=0;i<G.size();i++)
            {
                    P[i].resize(G.size());
                for (long j=0;j<G.size();j++)
                {
                    P[i][j].resize(G.size());
                    for (long k=0;k<G.size();k++) P[i][j][k].resize(G.size());
                }
            }

 
           for (long i=0;i<G.size();i++)
               for (long j=i;j<G.size();j++)
                  for (long k=0;k<G.size();k++)
                       for (long l=k;l<G.size();l++)
                             {
                                P[i][j][k][l]=0.0;
                                if (dom[k][l][i][j]!=0.0)
                                    {
                                for (long u=0;u<G.size();u++)
                                    for (long v=u;v<G.size();v++)
                                        P[i][j][k][l]+=dom[u][v][i][j]*pollen[u][v];

                               P[i][j][k][l]=dom[k][l][i][j]*pollen[k][l]
                                                *G[i][j]/(double) P[i][j][k][l];
                                    }
                               P[j][i][l][k]=P[j][i][k][l]=P[i][j][l][k]=P[i][j][k][l];
                             }

for (long m=0;m<G.size();m++)
    for (long n=m;n<G.size();n++)
        {
            G[m][n]=0.0;
            if (m==n)
                for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=m && k!= j&&j!=m)
                                G[m][m]+=0.125*(P[m][j][m][k]+P[m][k][m][j]);

                        G[m][m]+=0.5*(P[m][m][m][j]+P[m][j][m][m]);
                        if(j!=m) G[m][m]+=0.25*P[m][j][m][j];
                    }
            else
               {
                    for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=n && j!= m)
                                G[m][n]+=0.25*(P[m][j][n][k]+P[n][k][m][j]);

                        G[m][n]+=0.5*(P[m][m][n][j]+P[n][j][m][m])
                                +0.5*(P[m][j][n][n]+P[n][n][m][j]);
                    }
                }
            }
   for (long m=0;m<G.size();m++) for (long n=0;n<m;n++) G[m][n]=G[n][m];

}



void genotypes::fecondation_wFS(genotypes pollen, dominance dom, dominance_parameters dp)
  {
        long count=0, nci,ncj;
vector<vector <double > > Gtemp=G, EQhet, temphet;//EQhet contains genotypic frequencies at equilibrium for alleles in class {i,j} for i!=j (heterozygotes)
vector <double >EQhom, temphom;//EQhom contains genotypic frequencies at equilibrium for homozygotes
vector <vector <vector < vector <double > > > > P;
char ch;
double verif=0.0;

EQhet.resize(dp.get_size()); EQhom.resize(dp.get_size()); 
for (long i=0;i<dp.get_size();i++)
    {EQhet[i].resize(dp.get_size());EQhom[i]=-9.0;
       for (long j=0;j<dp.get_size();j++) EQhet[i][j]=-9.0; }
       
temphet=EQhet;
temphom=EQhom;

            P.resize(G.size());
            for (long i=0;i<G.size();i++)
            {
                    P[i].resize(G.size());
                for (long j=0;j<G.size();j++)
                {
                    P[i][j].resize(G.size());
                    for (long k=0;k<G.size();k++) P[i][j][k].resize(G.size());
                }
            }
        
for (long i=0;i<G.size();i++)
    for (long j=i;j<G.size();j++)
        for (long k=0;k<G.size();k++)
            for (long l=k;l<G.size();l++)
                {
                    nci=dp.get_nclass(i);ncj=dp.get_nclass(j);
                    
                    if (dom[k][l][i][j]!=0.0)
                        {
                            if(i==j && temphom[nci]==-9.0)
                           {
                                temphom[nci]=0.0;
                            for (long u=0;u<G.size();u++)
                             for (long v=u;v<G.size();v++)
                                temphom[nci]+=dom[u][v][i][j]*pollen[u][v];
                                }
                          
                            if(i!=j && temphet[nci][ncj]==-9.0)
                           {
                               temphet[nci][ncj]=0.0;
                            for (long u=0;u<G.size();u++)
                             for (long v=u;v<G.size();v++)
                                temphet[nci][ncj]+=dom[u][v][i][j]*pollen[u][v];
                                temphet[ncj][nci]=temphet[nci][ncj];
                                }
                         
                          if (i==j) P[i][j][k][l]=dom[k][l][i][j]*pollen[k][l]
                                                *G[i][j]/(double) temphom[nci];
                                else P[i][j][k][l]=dom[k][l][i][j]*pollen[k][l]
                                                *G[i][j]/(double) temphet[nci][ncj];
                        }
                        else P[i][j][k][l]=0.0;

                    P[j][i][l][k]=P[j][i][k][l]=P[i][j][l][k]=P[i][j][k][l];
                }
                
/*for (long i=0;i<G.size();i++)
    for (long j=i;j<G.size();j++)
        for (long k=0;k<G.size();k++)
            for (long l=k;l<G.size();l++)
                {
                    P[i][j][k][l]=0.0;
                    if (dom[k][l][i][j]!=0.0)
                        {
                           for (long u=0;u<G.size();u++)
                             for (long v=u;v<G.size();v++)
                                P[i][j][k][l]+=dom[u][v][i][j]*pollen[u][v];

                          P[i][j][k][l]=dom[k][l][i][j]*pollen[k][l]
                                                *G[i][j]/(double) P[i][j][k][l];
                        }
                       
                    P[j][i][l][k]=P[j][i][k][l]=P[i][j][l][k]=P[i][j][k][l];
                }*/


for (long m=0;m<G.size();m++)
    for (long n=m;n<G.size();n++)
        {
            G[m][n]=0.0;
            if (m==n)
            {
                if (EQhom[dp.get_nclass(m)]==-9.0)
                    {
                         for (long j=0;j<G.size();j++)
                            {
                                 for (long k=0;k<G.size();k++)
                                    if (k!=m && k!= j&&j!=m)
                                        G[m][m]+=0.125*(P[m][j][m][k]+P[m][k][m][j]);

                                        G[m][m]+=0.5*(P[m][m][m][j]+P[m][j][m][m]);
                                    if(j!=m) G[m][m]+=0.25*P[m][j][m][j];
                             }
                     EQhom[dp.get_nclass(m)]=G[m][m];
                    }
                else G[m][m]=EQhom[dp.get_nclass(m)];
            }
            else
                if (EQhet[dp.get_nclass(m)][dp.get_nclass(n)]==-9.0)
                {
                    for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=n && j!= m)
                                G[m][n]+=0.25*(P[m][j][n][k]+P[n][k][m][j]);

                        G[m][n]+=0.5*(P[m][m][n][j]+P[n][j][m][m])
                                +0.5*(P[m][j][n][n]+P[n][n][m][j]);
                    }
                    EQhet[dp.get_nclass(m)][dp.get_nclass(n)]=G[m][n];
                }
                else G[m][n]=EQhet[dp.get_nclass(m)][dp.get_nclass(n)];
                
                G[n][m]=G[m][n];
        }
}

/*void genotypes::fecondation_wFS(genotypes pollen, dominance dom, alleles a)
  {
char ch;
    vector<vector <double > > Gtemp=G;

        if (P.size()!=G.size())
        {
            P.resize(G.size());
            for (long i=0;i<G.size();i++)
            {
                    P[i].resize(G.size());
                for (long j=0;j<G.size();j++)
                {
                    P[i][j].resize(G.size());
                    for (long k=0;k<G.size();k++) P[i][j][k].resize(G.size());
                }
            }
        }


           for (long i=0;i<G.size();i++)
               for (long j=i;j<G.size();j++)
                  for (long k=0;k<G.size();k++)
                       for (long l=k;l<G.size();l++)
                             {
                                P[i][j][k][l]=0.0;
                            if (dom[k][l][i][j]!=0.0 && dom[k][l][i][j]!=0.0 && pollen[k][l]!=0.0&&G[i][j]!=0.0)
                                {
                                for (long u=0;u<G.size();u++)
                                    for (long v=u;v<G.size();v++)
                                        if (dom[u][v][i][j]!=0.0 && pollen[u][v] !=0.0)
                                            P[i][j][k][l]+=dom[u][v][i][j]*pollen[u][v];

                                P[i][j][k][l]=dom[k][l][i][j]*pollen[k][l]
                                                *G[i][j]/(double) P[i][j][k][l];
                                }
                               P[j][i][l][k]=P[j][i][k][l]=P[i][j][l][k]=P[i][j][k][l];
                             }

for (long m=0;m<G.size();m++)
    for (long n=m;n<G.size();n++)
        {
            G[m][n]=0.0;
            if (a[m]!=0.0 && a[n]!=0.0)
            {
            if (m==n)
                for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=m && k!= j&&j!=m)
                                G[m][m]+=0.125*(P[m][j][m][k]+P[m][k][m][j]);

                        G[m][m]+=0.5*(P[m][m][m][j]+P[m][j][m][m]);
                        if(j!=m) G[m][m]+=0.25*P[m][j][m][j];
                    }
            else
               {
                    for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=n && j!= m)
                                G[m][n]+=0.25*(P[m][j][n][k]+P[n][k][m][j]);

                        G[m][n]+=0.5*(P[m][m][n][j]+P[n][j][m][m])
                                +0.5*(P[m][j][n][n]+P[n][n][m][j]);
                    }
                }
            }
}
   for (long m=0;m<G.size();m++) for (long n=0;n<m;n++) G[m][n]=G[n][m];
}*/

void genotypes::fecondation_wFS(genotypes pollen, dominance dom,alleles a)
  {
    vector<vector <double > > Gtemp=G;
    vector <vector <vector < vector <double > > > > P;
    double temp;
    bool verif=false;
    char ch;

            P.resize(G.size());
            for (long i=0;i<G.size();i++)
            {
                    P[i].resize(G.size());
                for (long j=0;j<G.size();j++)
                {
                    P[i][j].resize(G.size());
                    for (long k=0;k<G.size();k++) P[i][j][k].resize(G.size());
                }
            }

for (long i=0;i<G.size();i++)
    for (long j=i;j<G.size();j++)
        {
                temp=0.0;
                for (long u=0;u<G.size();u++)
                    for (long v=u;v<G.size();v++)
                                temp+=dom[u][v][i][j]*pollen[u][v];
                if (temp==0.0)  verif=true;
            for (long k=0;k<G.size();k++)
            for (long l=k;l<G.size();l++)
                {
                    if (dom[k][l][i][j]!=0.0 && temp!=0.0)
                          P[i][j][k][l]=dom[k][l][i][j]*pollen[k][l]
                                                *G[i][j]/(double) temp;
                                                
                        else P[i][j][k][l]=0.0;
                    P[j][i][l][k]=P[j][i][k][l]=P[i][j][l][k]=P[i][j][k][l];
                }
}

for (long m=0;m<G.size();m++)
    for (long n=m;n<G.size();n++)
        {
            G[m][n]=0.0;
            if (a[m]!=0.0 && a[n]!=0.0)
            {
            if (m==n)
                for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=m && k!= j&&j!=m)
                                G[m][m]+=0.125*(P[m][j][m][k]+P[m][k][m][j]);

                        G[m][m]+=0.5*(P[m][m][m][j]+P[m][j][m][m]);
                        if(j!=m) G[m][m]+=0.25*P[m][j][m][j];
                    }
            else
               {
                    for (long j=0;j<G.size();j++)
                    {
                        for (long k=0;k<G.size();k++)
                        if (k!=n && j!= m)
                                G[m][n]+=0.25*(P[m][j][n][k]+P[n][k][m][j]);

                        G[m][n]+=0.5*(P[m][m][n][j]+P[n][j][m][m])
                                +0.5*(P[m][j][n][n]+P[n][n][m][j]);
                    }
                }
            }
}
   for (long m=0;m<G.size();m++) for (long n=0;n<m;n++) G[m][n]=G[n][m];
   temp=0.0;
   if (verif) for (long m=0;m<G.size();m++) for (long n=0;n<=m;n++) temp+=G[m][n]; 
   if (verif && temp==0.0) extinction=true;
}

alleles::alleles(genotypes g)
{
    A.resize(g[0].size());
    for(long i=0;i<A.size();i++)
        {A[i]=0.0;
        for(long j=0;j<A.size();j++)
            if (i!=j)A[i]+=0.5*g[i][j];
                else A[i]+=g[i][j];
        }
    }

/****************************************************************************/
/***********************     METAPOPULATION    ******************************/
/****************************************************************************/

metapopulation::metapopulation(genotypes g, double mig, double mut, long N, long d)
{
    mp.resize(d,g);
    m=mig; mu=mut;Ni=N; Nd=d;
    //pol_tot=g;
    }
 
void metapopulation::fecondation_FS(dominance matrixD)
{
    genotypes pol;
    alleles a;
    vector<genotypes > mptemp=mp;
      //char ch;

    for (long i=0;i<mp.size();i++)
        {
            pol=islands_migration_model(i);
            a=alleles(pol);
            
        /*    cout << "après migration" << endl;
             pol.display();
            cout << endl;
            cin >> ch;*/

            mptemp[i].fecondation_FS(pol,matrixD,a);
        }

        mp=mptemp;
}

void metapopulation::fecondation_wFS(dominance matrixD)
{
    genotypes pol;
    alleles a;
    vector<genotypes > mptemp=mp;
  

    for (long i=0;i<mp.size();i++)
        {
            pol=islands_migration_model(i);
            a=alleles(pol);

            mptemp[i].fecondation_wFS(pol,matrixD,a);
        }
        
        mp=mptemp;
}
 
genotypes metapopulation:: islands_migration_model(long i)
{
    return genotypes(mp, i, m);
    }
    
double metapopulation::global_Fst()
{
vector<alleles > va;
alleles meanf((alleles(mp[0])).size());
va.resize(mp.size());
double Q3=0.0, Q2=0.0;

    for (long i=0;i<mp.size();i++)
    {
        va[i]=alleles(mp[i]);
        for (long k=0;k<meanf.size();k++)
            {
                meanf[k]+=va[i][k];
                Q2+=va[i][k]*va[i][k];
            }
    }

 for (long k=0;k<meanf.size();k++) {meanf[k]/=mp.size(); Q3+=meanf[k]*meanf[k];}
 Q2/=mp.size();
 
 return (Q2-Q3)/(1.0-Q3);

    }
    
alleles metapopulation::mean_allelic_frequency()
{
        vector <alleles > va(mp.size());
        
        for (long i=0;i<mp.size();i++) va[i]=alleles(mp[i]);
        
        return alleles(va);
}
    
void metapopulation::file_out(long gen )
{
    mean_allelic_frequency().file_out(gen);
    /*ofstream myFile("Fst",ios::app);
    if (! myFile) cout << "Problèmes à l'ouverture du fichier de sortie Fst";
myFile << gen << '\t';
    //for (long i=0;i<A.size();i++) myFile << A[i] << '\t';
    myFile << endl; myFile.close();*/

    }

void metapopulation::final_file_out()
{
         for (long i=0;i<mp.size();i++) (alleles(mp[i])).final_file_out();
    }

/************************************************************************************/
/************************************************************************************/
/***********************************Dialog Class*************************************/
/************************************************************************************/
/************************************************************************************/

 user_parameters::user_parameters()
{
     int choice_type=0, choice=0,est=0;
     _simple=true;
     _is_model_FS=false; _is_model_wFS=true;
     _is_sampling=false;
     _allelic_frequency_distribution=false;
     _genotypic_frequency_distribution=false;
     class_allelic_frequency_number=-1;
     generations_between_sampling=-1;
     stop=0.0001;
     _one_generation_change=false;
     char ch;
     Nsample=-1;
     string line;
  time_t rawtime;
  struct tm * timeinfo;
  
  long graine1;
     
     int count_line=0;
     
     string temp_string;
     
     
     user_file="user_parameters.txt";
     
     
     
     
     cout << "                      *************************"<<endl<<endl
          << "    Welcome in NESSI:" << endl << "    Numerical Estimations for Sporophytic Self-Incompatibility." << endl << endl
          << '\t' <<  "program written by Sylvain Billiard." << endl
          << '\t' << "see Billiard S, V. Castric and X. Vekemans 2007 Genetics  for details."
          << endl <<endl << "                      *************************";
     cout <<endl << endl << endl << endl;
     
cout << "Avant de commencer, donnez moi une graine pour le generateur de nombre aleatoire svp (entier quelconque). C'est pour aider au débugage, merci." << endl << endl
<< "Graine : "; cin >> graine1;     
     al.seed(graine1);

cout <<endl << endl << endl << endl;     
     
cout << " *** Would you like to compute: "<< endl
     << "    (1) Deterministic genotypic and allelic equilibrium frequencies? " << endl
     << " or (2) Distributions in finite populations? " << '\t' << endl
     << " or (3) Genotypic and Allelic frequencies change distributions in a generation ? ";

cout << endl << endl;

while (choice_type!=1 && choice_type !=2 && choice_type !=3 && choice_type!=9 ) cin >> choice_type; 

if (choice_type ==1) {  
     _derive=false; _mutation=false; _deterministe=true;
     }

if (choice_type==2) {
_derive=true; _mutation=true; _deterministe=false;_allelic_frequency_distribution=true; _genotypic_frequency_distribution=true;
                    }

if (choice_type==3) {
_derive=true; _mutation=true; _deterministe=false; _one_generation_change=true;_allelic_frequency_distribution=true;_genotypic_frequency_distribution=true;
gen_max=1;
generations_between_sampling=0;
}

               
if (choice_type==9)
{
_derive=true; _mutation=true; _deterministe=false;_one_generation_change=false;
}               
  
if (choice_type !=9){               

cout << " *** Dominance relationships : "<< endl
     << "    (1) Simple: dom or domcod models?" << endl
     << " or (2) Specific?" << '\t';

cout << endl << endl;

while (choice!=1 && choice !=2) cin >> choice; 
if (choice ==2) _simple=false; 

cout << endl << endl;

choice=0;

cout << " *** Frequency-dependent selection model through: "<< endl
     << "    (1) Male way only (Wright's model)?" << endl
     << " or (2) Male and female ways (fecundity selection)?" << '\t';
while (choice!=1 && choice !=2) cin >> choice; 

if (choice==2) {_is_model_FS=true;_is_model_wFS=false;}

}

cout << endl << endl;  

while (temp_string != "y") {cout << "Is the input parameters file (user_parameters.txt) ready for reading? (y/n)  "  ; cin >> temp_string;}

std::ifstream FileIn(user_file.c_str());
     if (!FileIn){cout << "Error opening user_parameters.txt file" << endl; cin >> ch;exit(1);}
     

while(!getline(FileIn, line).eof())
{
  getline(FileIn, temp_string); 
  
  switch (count_line)
{
case 0: if (strtod(temp_string.c_str(),NULL)>10e-15) stop = strtod(temp_string.c_str(),NULL);break;
case 1: class_allelic_frequency_number = strtol(temp_string.c_str(),NULL,0);break;
case 2: Nindividuals = strtol(temp_string.c_str(),NULL,0);break;
case 3: mu = strtod(temp_string.c_str(),NULL);break;
case 4: threshold_ANb = strtod(temp_string.c_str(),NULL);break;
case 5: oubli = strtol(temp_string.c_str(),NULL,0);break;
case 6: gen_max = strtol(temp_string.c_str(),NULL,0);break;
case 7: repet_max = strtol(temp_string.c_str(),NULL,0);break;
case 8: Nsample = strtol(temp_string.c_str(),NULL,0);break;
case 9: generations_between_sampling = strtol(temp_string.c_str(),NULL,0);break;
case 10: initial_genotypic_frequencies=temp_string;break;
case 11: pollen_dominance=temp_string;break;
case 12: pistil_dominance=temp_string;break;
case 13: parameters_file=temp_string;break;
default: break;
}

count_line++;
     
}

FileIn.close();

if (class_allelic_frequency_number==0){ _allelic_frequency_distribution=false;_genotypic_frequency_distribution=false;}
if (Nsample!=0 && Nsample!=Nindividuals) _is_sampling=true;

 time ( &rawtime );
 timeinfo = localtime ( &rawtime );
  
    
ofstream FileOut(user_file.c_str());
if (!FileOut){cout << "Error opening user_parameters.txt file (it must be closed) " << endl; cin >> ch;}
   FileOut << "//*** Genotypic frequencies equilibrium criteria: 10e-x, x = ? (0<x<9)"<<endl;
           if (choice_type==1) FileOut << stop;  
   FileOut << endl << "//*** Number of frequency classes for frequencies distribution estimation? (0=no distribution)"<<endl;
           if (choice_type==2||choice_type==3) FileOut << class_allelic_frequency_number;
   FileOut << endl << "//*** Number of diploid individuals in the populations?"<<endl;
           if (choice_type==2||choice_type==3) FileOut << Nindividuals;
   FileOut << endl << "//*** Mutation rate (per individual per generation)?"<<endl ;
           if (choice_type==2||choice_type==3) FileOut << mu;
   FileOut << endl << "//*** Minimal frequency of an allele to be counted?"<<endl;
           if (choice_type==2||choice_type==3) FileOut << threshold_ANb;
   FileOut << endl << "//*** Number of forget generations?"<<endl;
           if (choice_type==2||choice_type==3) FileOut << oubli;
   FileOut << endl << "//*** Number of generations in a single repetition (after forget generations)?"<<endl;
           if (choice_type==2||choice_type==3) FileOut << gen_max;
   FileOut << endl << "//*** Number of repetitions?"<<endl;
           if (choice_type==2||choice_type==3) FileOut << repet_max;
   FileOut << endl << "//*** Sample size? "<<endl;
           if (choice_type==2||choice_type==3) FileOut << Nsample;
   FileOut << endl << "//*** Number of generations between sampling?"<<endl;
           if (choice_type==2||choice_type==3) FileOut << generations_between_sampling;
           
   FileOut << endl << "//*** File name for initial (observed) genotypic frequencies? (with extension)"<<endl;
           if (choice_type==3) FileOut << initial_genotypic_frequencies;
   FileOut << endl << "//*** File name for dominance relationships in pollen? (with extension)"<<endl;
           if (!_simple) FileOut << pollen_dominance;
   FileOut << endl << "//*** File name for dominance relationships in pistil? (with extension)"<<endl;
           if (!_simple) FileOut << pistil_dominance;
   FileOut << endl << "//*** File name for simple dominance relationships? (with extension)"<<endl;
           if (_simple) FileOut << parameters_file;
   FileOut << endl << endl << "//*** First seed for the pseudo-random number generator: "<< graine1 <<endl;
           //<< "//*** Second seed for the pseudo-random number generator: " << graine2 <<endl;

FileOut <<endl<<endl
<< "// Last Analysis: " <<  asctime (timeinfo)<< endl<< endl
<< "//      Computations: ";
if (choice_type==1) FileOut << "Deterministic genotypic and allelic equilibrium frequencies." << endl;
if (choice_type==2) FileOut << "Distributions in finite populations." << endl;
if (choice_type==3) FileOut << "Genotypic and Allelic frequencies change distributions in a generation." <<endl;
if (choice_type==9) FileOut << "It was a TEST." << endl;

FileOut << "//      Type of frequency dependent selection: ";
if (_is_model_FS) FileOut << "Male and female ways (fecundity selection)." << endl;
 else FileOut << "Male way only (Wright's classical model)."<< endl;
 
FileOut << "//      Dominance relationships: ";
if (_simple) FileOut << "simple (dom, domcod or cod)"<< endl;
else FileOut << "Specific (dominance matrices)"<< endl;
     FileOut.close();
     
     }

