// Author: Wes Kendall
// Copyright 2011 www.mpitutorial.com
// This code is provided freely with the tutorials on mpitutorial.com. Feel
// free to modify it for your own use. Any distribution of the code must
// either provide a link to www.mpitutorial.com or keep this header intact.
//
// MPI_Send, MPI_Recv example. Communicates the number -1 from process 0
// to process 1.
//
#include <mpi.h>
#include <iostream>
#include <stdio.h>      /* printf, scanf, puts, NULL */
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <fstream>
#include <cmath>
#include <sys/resource.h>
#include <vector>
#include<string>
using namespace std;

/* come ising 2d ma adesso tipo blume capel,punto di reticolo vale 0 se non c'è spin
si calcolano una volta all'inizio possibili valori di coseno di differenza angoli, così dopo basta prendere elemento che 
si vuole di questo vettore. bondenergy dà solo energia tra i e j dovuta a coseno */


// PARAMETRI DA VARIARE
const int p=8; // numero angoli
const int Lx=10+1;
const int Ly=(Lx-1)*4;

int tempo=36;
int tempotermal=tempo/2;
//const double numgiri=1;
const double D= 1.02;// piùmeno 60 è abbastanza per rendere tutto vuoto, cioè 0, (se D negativo) 
double beta=.5637963; // beta critico per D=1.02 è beta=0.5637963; beta=5 è praticamente T=0, andare oltre è peggio
// per ising D critico e' -0.655, beta crit=
double Dfrattobeta=D/beta;

// PARAMETRI CHE SI CALCOLA LUI
const int Lz=Ly;        // SE SI CAMBIA QUESTO BISOGNA CAMBIARE NEXT-PREVIOUS E PURE GENERATORE RANDOM DI SITO DA 1 A LY
const int N=Lx*Ly*Lz;   const int sitivariabili=Lz*Ly*(Lx-2);
const double pi= 4.*atan(1);
int fr=1;   // serve per reticolo di corr2
int frSUX=(Lx-1)/16;
string nomem="correlazionimclockELLEICS.txt";
string nomee="correlazioniVICINECENTRALIeclockELLEICS.txt";
// PESI CHE VENGONO SETTATI DA SETTAPESI
long int rng;
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


#include "funzioniclockold.h"

void randfill()
{
  for(int k=0;k<Ly;k++)
    {
      for(int j=0;j<Ly;j++)
     {
       m[0][j][k]=megaspin; m[Lx-1][j][k]=megaspin;    
       for(int i=1;i<Lx-1;i++)
       {
          if(randombit()==1) m[i][j][k]=0; 
          else m[i][j][k]=randomspin();
       }
      }
    }
}

void filldafile()
{ 
  ifstream nomequalsiasi;
  nomequalsiasi.open("aReticoloTermalizzatoELLEICS.txt");
  for(int i=0;i<Lx;i++)
    {
      for(int j=0;j<Ly;j++)
     {
       for(int k=0;k<Lz;k++)
       {
        nomequalsiasi>>m[i][j][k];
       }
      }
    }
  nomequalsiasi.close();
  megaspin=m[0][0][0];
}


int  main ()
{

    int tiniz=time(NULL);
    
    
    settapesi(beta);
     

  double mediafinale=0;   // magnetizzazione complessiva
  double magnxfinale[Lx/2+1]={0};   // magnetizz in funzione della distanza dal bordo
  //  double corr2finale[16][16][16][16]={0};
  // double meansqcorr2[16][16][16][16]={0};
  double meansq[Lx/2+1]={0};           // mediaquadra di diversi giri MC per avere varianza
  double exfinale[Lx/2+1]={0};
  double meansqE[Lx/2+1]={0};

  
        rng=-time(NULL);     // settare entrambi in modo che dipenda da che core si sta usando
        srand (time(NULL) );

    int stepdopotermalizz=1;
        // correlazioni dell'energia hanno variabili con stessi nomi di corr magn ma con E davanti
        
       double mediatemp=0;
       double magnxtemp[Lx/2+1]={0};
       double extemp[Lx]={0};

       //   double corr2temp[16][16][16][16]={0};


       
       randfill(); 
       megaspin=m[0][0][0];
       // randfill();
          
          //ofstream nomequalsiasi; nomequalsiasi.open("dati.txt");
	 cout<<"inizia"<<time(NULL)-tiniz<<endl;
	  while(time(NULL)-tiniz<tempotermal)
	    {
	      SweepEWolff();
	    }
	  while(time(NULL)-tiniz>=tempotermal&&time(NULL)-tiniz<tempo) //for(int mosse=0;mosse<nmontecarlo; mosse++) //for(int mosse=0;mosse<nmontecarlo;mosse++)
	    {
	      SweepEWolff();
         
           for(int ii=0;ii<Lx/2+1;ii++) 
              {
                extemp[ii]+=(energiax(ii)+energiax(Lx-1-ii))/2;
                magnxtemp[ii]+=(magnx(ii)+magnx(Lx-1-ii))/2;
              }
	   //   corr(corr2temp);            
	     stepdopotermalizz++;
	    }
	  cout<<"giri post termal "<<stepdopotermalizz<<endl;
             // finito un giro di montecarlo
         for(int ii=0;ii<Lx/2+1;ii++) 
              {
                magnxtemp[ii]/=(stepdopotermalizz-1);
                extemp[ii]/=(stepdopotermalizz-1);
              }        
          mediatemp/=(stepdopotermalizz-1);

      for(int ii=0;ii<Lx/2+1;ii++) 
          {
            cout<<magnxtemp[ii]<<" , ";;
          }
          cout<<endl<<time(NULL)-tiniz;
  return 0;
}
    


