#include <iostream>
#include <stdlib.h>
#include "MersenneTwister.h"
#include <vector>
#include<cmath>
#include<string.h>
#include "ssi.h"
#include "randomc.h"                   // define classes for random number generators
#include "stocc.h"
#include <bitset>

using namespace std;
MTRand al(212);

//int32 graine2 =time(0);                // random seed
int32 graine2=3211;
StochasticLib1 sto(graine2);     // defines a random generator to use for specific distributions

int main(int argc, char *argv[])
{
 alleles a, at;
 sheltered_load_loci sll;
vector <sheltered_load_loci> tab_distrib;
     string ch;   
    dominance_parameters dom;
    user_parameters up;
    dominance matrixD ;
      
    dp_table dp;
    genotypes pol, g, geq, gsample;
    long gen, step=1000, step_count=1,step_count_sampling=1;
               
    int cases_number;
 //   string ch;
    double threshold=0.0;
    double disc,ret;

    ofstream myFile1("Genotypes.txt"), myFile2("Alleles.txt"),myFile3("F_Alleles.txt")
                ,myFile4("F_Genotypes.txt"), myFile5("Alleles_Number.txt"),myFile6("Allelic_Frequency_Distribution.txt"),
                 myFile7("Genotypic_Frequency_Distribution.txt"), myFile8("Allelic_Frequency_Distribution_ID.txt"), myFile9("load_ID.txt")
                 , myFile10("fixation_number_ID.txt"),myFile11("distribution_homozygotes_ID.txt");
    if (! myFile1 || ! myFile2 || !myFile3 || !myFile4 ||  !myFile5||  !myFile6 || !myFile7|| !myFile8|| !myFile9|| !myFile10|| !myFile11){cout << "Error opening output file" << endl; cin >> ch;return -1;}
    myFile1.close(); myFile2.close();myFile3.close();myFile4.close();myFile5.close();myFile6.close();myFile7.close();myFile8.close();myFile9.close();myFile10.close();myFile11.close();
  
    if (up._simple) {  dp=dp_table(up.parameters_file);
cout << "Nombre de combinaisons à effectuer : "<< dp.size() << endl<< endl;
cases_number=dp.size();}
  
else {matrixD=dominance(up.pollen_dominance.c_str(),up.pistil_dominance.c_str(),0);cases_number=1; dom=dominance_parameters(matrixD.get_allele_number());}

 

for (long k=0;k<cases_number;k++)
{                                   //1
gen=0;
step_count=1;

if (up._simple)
   {               
   dom=dp[k];
   matrixD=dominance(dom,threshold);
   dom.file_out();
   cout << dp.size()-k << '\t' << endl;  
   }

else  matrixD.display();
        
if (!up._one_generation_change)  
        {g=genotypes(matrixD.get_allele_number());
        g.initialisation();}
else {geq=genotypes(up.initial_genotypic_frequencies.c_str()); g=geq;g.display();}

a=alleles(g);

if (up._one_generation_change)  a.file_out(1);

at=a;
  
/********************************************************************/
/*******************Deterministic Computations***********************/
/********************************************************************/

if (!up._one_generation_change)
   while(gen==0||(a-at>up.stop)) 
   {   pol=g;

    if (gen==0) gen++;
    if (up._is_model_FS) {if (up._simple) g.fecondation_FS(pol, matrixD, dom); else g.fecondation_FS(pol, matrixD);}
     else {if (up._simple) g.fecondation_wFS(pol, matrixD,dom); else g.fecondation_wFS(pol, matrixD);}
    at=a;
    a=alleles(g);
    geq=g;
    g.file_out(gen);
    a.file_out(gen);
    }
    
if(!up._derive) {a.final_file_out(); if (up._simple) g.final_file_out(dom); else g.final_file_out();}

/********************************************************************/
/*********************Stochastic Simulations*************************/
/********************************************************************/
if (up._derive) 
   {  
   
             for (long ii=1;ii<=up.repet_max;ii++)
   {    
       bool t=true;
        step_count=1;   step_count_sampling=0;
        gen=0;
        g=geq;
         sll =sheltered_load_loci(up.loci_number_ID,geq.size(),up.selection_coefficient, alleles(g),up.Nindividuals );
        cout << "Repetition number "<< ii << endl;
   while(gen==0|| gen<up.oubli+up.gen_max+1) 
      
      {     
            
             
           /*  if (!g._is_extinct()) sll=g.derive_ID(up.Nindividuals,sll);
                 sll.distrib_all(); */
            
            if (gen>up.oubli+step_count_sampling*up.generations_between_sampling) {
                              if (up._is_sampling) {gsample=g; gsample.derive(up.Nsample);}
                                 else gsample=g;
                                 
                               a=alleles(gsample);                                
                               dom.check_allele_number(a, up.threshold_ANb);
                               dom.check_allele_distrib();
                               if (up._allelic_frequency_distribution) dom.check_allelic_frequency_distrib(a, up.class_allelic_frequency_number) ;
                               if (up._genotypic_frequency_distribution) dom.check_genotypic_frequency_distrib(gsample, up.class_allelic_frequency_number) ;
                               step_count_sampling++;
                              }

            
           if (!g._is_extinct()) 
            { if(up._mutation) g.mutation_KAM(up.mu,up.Nindividuals);
            
              
            a=alleles(g);
            if(up._mutation_ID) sll.mutation(a,up.mu_ID,up.eta_ID,up.Nindividuals);
        
            pol=g;
            
            if (up._is_model_FS) g.fecondation_FS(pol, matrixD, a);
              else g.fecondation_wFS(pol, matrixD, a);
                  
            //a=alleles(g);
            
            
            if (!g._is_extinct()) sll=g.derive_ID(up.Nindividuals,sll);
                 sll.distrib_all();
                 
            a=alleles(g);     
            }
            if (gen==step_count*step){ cout << gen << endl;  step_count++; /*sll.fixation_number(ii,t );tab_distrib.push_back(sll);t=false;sll.distribution_homozygotes_ID(ii);*/sll.display(ii,gen);}
                        
                                
                        gen++;
                        
                       
        } 
       
 
        
  }

dom.ADistrib_file_out();
    if (up._allelic_frequency_distribution) dom.AFDistrib_file_out();
    if (up._genotypic_frequency_distribution) dom.GFDistrib_file_out();
 
  
    
for (int k=0;k<a.size();k++)
{
    ofstream myFile("Allelic_Frequency_Distribution_ID.txt",ios::app);
    
    myFile << "Allele :" << k <<endl;
    
    for (int kk=0;kk<tab_distrib.size();kk++) 
    tab_distrib[kk].ADistrib_file_out(k);
    
    myFile << endl;
    myFile.close();}
}  }        


  system("PAUSE");	
  return 0;
}

