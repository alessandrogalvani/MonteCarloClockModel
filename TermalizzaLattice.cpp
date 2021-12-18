// NON FA MISURE MA EVOLVE SISTEMA PER UN PO' E POI QUANDO SI SPERA SIA TERMALIZZATO STAMPA SU FILE LA MATRICE.


#include <iostream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <fstream>
#include <cmath>
#include <sys/resource.h>
#include <vector>
using namespace std;

/* come ising 2d ma adesso tipo blume capel,punto di reticolo vale 0 se non c'è spin
si calcolano una volta all'inizio possibili valori di coseno di differenza angoli, così dopo basta prendere elemento che 
si vuole di questo vettore. bondenergy dà solo energia tra i e j dovuta a coseno */


// PARAMETRI DA VARIARE
long int rng=-time(NULL);
const int p=8; // numero angoli
const int Lx=160+1;
const int Ly=(Lx-1)*6;
const int nmontecarlo=1;
const double D= 1.02;// piùmeno 60 è abbastanza per rendere tutto vuoto, cioè 0, (se D negativo) 
double beta=.5637963; // beta critico per D=1.02 è beta=0.5637963; beta=5 è praticamente T=0, andare oltre è peggio
// per ising D critico e' -0.655, beta crit=
double Dfrattobeta=D/beta;
int fr=1;

// PARAMETRI CHE SI CALCOLA LUI
const int Lz=Ly;        // SE SI CAMBIA QUESTO BISOGNA CAMBIARE NEXT-PREVIOUS E PURE GENERATORE RANDOM DI SITO DA 1 A LY
const int N=Lx*Ly*Lz;   const int sitivariabili=Lz*Ly*(Lx-2);
const double pi= 4.*atan(1);


// PESI CHE VENGONO SETTATI DA SETTAPESI

double expD=exp(D);
double expmD=1./expD;
double v[p];
double vseni[p];
double pesi[p];
double pesim[p];// pesi di boltzmann per vari possibili valori di angolo e loro inversi, per non ricalcolarli
double pbond[2*p][2*p];

static bool cluster[Lx][Ly][Lz];
static int m[Lx][Ly][Lz];

int megaspin=1; // da 1 a p;

#include "funzioniclock.h"

void randfill()
{
  for(int k=0;k<Ly;k++)
    {
      for(int j=0;j<Ly;j++)
	   {
	     m[0][j][k]=megaspin; m[Lx-1][j][k]=megaspin;    
	     for(int i=1;i<Lx-1;i++)
	     { m[i][j][k]=megaspin;
	        //if(randombit()==1) m[i][j][k]=0; // E' PER NON AVERE ZERI, CORREGGERE 
	        //else m[i][j][k]=randomspin();
	     }
	    }
    }
}

void  tutto ()
{
  srand (time(NULL));


  double mediafinale=0;   // magnetizzazione complessiva
  double magnxfinale[Lx]={0};   // magnetizz in funzione della distanza dal bordo
  double meansq[Lx]={0};           // mediaquadra di diversi giri MC per avere varianza
  double s4finale=0;
  double binder=0;
   
 
            
      randfill();
      int i=1; int j=0; int k=0;

      //ofstream nomequalsiasi; nomequalsiasi.open("dati.txt");
      for(int mosse=0;mosse<nmontecarlo; mosse++) //for(int mosse=0;mosse<nmontecarlo;mosse++)
	{

	  for(int sweep=0;sweep<sitivariabili;sweep++)
	    {
	      spinflip2(i,j,k);
	      sitodopo(i,j,k);
	    }
      
	  for(int rr=0;rr<6;rr++)
	    {
	      for(int sweep=0;sweep<sitivariabili;sweep++)
      		{
      		  spinflip1(i,j,k);
      		  sitodopo(i,j,k);
      		}
       //printpiano(1);
        wolff();
       // cout<<megaspin<<m[0][1][2]<<m[Lx-1][4][6]<<endl;
        //printcluster(1);
	    }

	  
	      //cout<<mean()<<endl;

                     
            
               
	  // conta[m[2][3][1]]++; }
               
	}      
      //
         /*   ofstream nomequalsiasi;
    nomequalsiasi.open("ReticoloTermalizzato.txt");
    for(int i=0;i<Lx;i++)
       for(int j=0;j<Ly;j++)
        for(int k=0;k<Lz;k++)
        nomequalsiasi << m[i][j][k]<<" ";
        
    nomequalsiasi.close();
*/


}
    


int main()
{ 

    int tiniz=time(NULL);
    
    
    settapesi(beta);
     cout<< "dimensione del sistema e'" <<Lx<<"x"<<Ly<<"x"<<Lz<<endl
        <<"numero di sweep e' "<<nmontecarlo<<endl;
       
    
    tutto();
    int tfin=time(NULL);
    cout<<endl<< "tempo trascorso e'"<<tfin-tiniz<<endl;

  /* for(int tt=0;tt<20;tt++)
        {
            double bet=.28+3*double(tt)/100;
            settapesi(bet);
           
	    cout<<endl<<"{ "<<bet<<", "; tutto();cout<< "},";
            
        }
  */

  

  /*    ofstream nomequalsiasi;
    nomequalsiasi.open("dati.txt");
    for(double i=1;i<8;i++)
        { beta+=.04; settapesi();
        nomequalsiasi <<"{"<<beta <<" , "<< tutto()<<" },"<< endl;
        }
    nomequalsiasi.close();
  */
   

  return 0;
}
