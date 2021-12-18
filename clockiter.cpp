// DOVREBBE ESSERE QUELLO DEFINITIVO

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
const int Lx=128+1;
const int Ly=(Lx-1)*6;
const int nmontecarlo=3;
const int steptermalizz=1;
const double numgiri=1;
const double D= 1.02;// piùmeno 60 è abbastanza per rendere tutto vuoto, cioè 0, (se D negativo) 
double beta=.5637963; // beta critico per D=1.02 è beta=0.5637963; beta=5 è praticamente T=0, andare oltre è peggio
// per ising D critico e' -0.655, beta crit=
double Dfrattobeta=D/beta;


// PARAMETRI CHE SI CALCOLA LUI
const int Lz=Ly;        // SE SI CAMBIA QUESTO BISOGNA CAMBIARE NEXT-PREVIOUS E PURE GENERATORE RANDOM DI SITO DA 1 A LY
const int N=Lx*Ly*Lz;   const int sitivariabili=Lz*Ly*(Lx-2);
const double pi= 4.*atan(1);
const int stepdopotermalizz=nmontecarlo-steptermalizz;
int fr=(Lx-1)/16;   // serve per reticolo di corr2
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
  nomequalsiasi.open("ReticoloTermalizzato.txt");
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
int tempotemp;
void  tutto ()
{


  double mediafinale=0;   // magnetizzazione complessiva
  double magnxfinale[Lx]={0};   // magnetizz in funzione della distanza dal bordo
  double meansq[Lx]={0};           // mediaquadra di diversi giri MC per avere varianza
  double s4finale=0;
  double binder=0;
  //double correlazioni[16][16][16][16]={0};
  double exfinale[Lx]={0};
  for(int iciclogrosso=0;iciclogrosso<numgiri;iciclogrosso++)
    {//cout <<"iniziato il giro "<<1+ iciclogrosso<< endl;
      double mediatemp=0;        
      double magnxtemp[Lx]={0};
      double extemp[Lx]={0};

      //filldafile();
      randfill();
      
      
      //ofstream nomequalsiasi; nomequalsiasi.open("dati.txt");
      for(int mosse=0;mosse<nmontecarlo; mosse++) //for(int mosse=0;mosse<nmontecarlo;mosse++)
  {
           SweepEWolff();
           cout<<mosse<<" "<<mean()<<endl;
    

     if(mosse>steptermalizz) 
      { 
        //corr(correlazioni);

        for(int ii=1;ii<Lx/2+1;ii++) 
        {
          magnxtemp[ii]+=(magnx(ii)+magnx(Lx-1-ii))/2;
          
        }
        for(int ii=0;ii<Lx/2+1;ii++) 
        {
          extemp[ii]+=(energiax(ii)+energiax(Lx-1-ii))/2;
          
        }
                             
        // mediatemp+= mean(m);
        //  s4temp+=s4bulk(m);
      }
                
    // conta[m[2][3][1]]++; }
               
  }       // finito un giro di montecarlo
     
            
      //mediatemp/=(stepdopotermalizz-1);
      // s4temp/=(stepdopotermalizz-1);
      //mediafinale+=mediatemp;
      //binder+=s4temp/pow(mediatemp,2);
                
      for(int u=1;u<Lx/2+1;u++) 
  {
    magnxfinale[u]+=magnxtemp[u]/(stepdopotermalizz-1);
   // meansq[u]+=pow(magnxtemp[u]/(stepdopotermalizz-1),2);
  }
         for(int u=0;u<Lx/2+1;u++) 
  {
    exfinale[u]+=extemp[u]/(stepdopotermalizz-1);
   // meansq[u]+=pow(magnxtemp[u]/(stepdopotermalizz-1),2);
  }
            

      //cout <<binder<<endl;
      //mediafinale+=binder;
    }
    
  /*ofstream nomequalsiasi;
  nomequalsiasi.open("dati24con1E4step2giri.txt"); 
  for(int u=1;u<Lx/2+1;u++)
    {
      nomequalsiasi<<magnxfinale[u]/numgiri<<" , "
       <<sqrt(meansq[u]/numgiri-pow(magnxfinale[u]/numgiri,2))<<endl; 
    }
    
           
  nomequalsiasi.close();
   */
    for(int u=1;u<Lx/2+1;u++)
    {
      cout<<magnxfinale[u]/numgiri<<" , ";
       
    }cout<<endl<<"energie sono:";
    for(int u=0;u<Lx/2+1;u++)
    {
      cout<<exfinale[u]/numgiri<<" , ";
       
    }cout<<endl;
    /*cout<<endl<<"le deviazioni standard sono"<<endl;
     for(int u=1;u<Lx/2+1;u++)
    {
      cout<<sqrt(meansq[u]/numgiri-pow(magnxfinale[u]/numgiri,2))<<" , "; 
    }*/
    //tabella(correlazioni,(stepdopotermalizz-1)*numgiri);
}
    


int main()
{

        int tiniz=time(NULL);
    srand (time(NULL));
    
    settapesi(beta);
     cout<< "dimensione del sistema e'" <<Lx<<"x"<<Ly<<"x"<<Lz<<endl
        <<"numero di sweep e' "<<nmontecarlo<<endl
        <<"numero di giri e' "<<numgiri<<endl;
    
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
