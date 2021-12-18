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
const int Lx=32+1;
const int Ly=(Lx-1)*6;
const int nmontecarlo=50;
const int steptermalizz=33;
//const double numgiri=1;
const double D= 1.02;// piùmeno 60 è abbastanza per rendere tutto vuoto, cioè 0, (se D negativo) 
double beta=.5637963; // beta critico per D=1.02 è beta=0.5637963; beta=5 è praticamente T=0, andare oltre è peggio
// per ising D critico e' -0.655, beta crit=
double Dfrattobeta=D/beta;

// PARAMETRI CHE SI CALCOLA LUI
const int Lz=Ly;        // SE SI CAMBIA QUESTO BISOGNA CAMBIARE NEXT-PREVIOUS E PURE GENERATORE RANDOM DI SITO DA 1 A LY
const int N=Lx*Ly*Lz;   const int sitivariabili=Lz*Ly*(Lx-2);
const int Ly2=Ly*Ly;
const double pi= 4.*atan(1);
const int stepdopotermalizz=nmontecarlo-steptermalizz;
int fr=(Lx-1)/16;   // serve per reticolo di corr2

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


int  main ()
{
    MPI_Init(NULL, NULL);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  if (world_size < 2) 
  {
    cout<<"non ci sono abbastanza core";
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

    int tiniz=time(NULL);
    
    
    settapesi(beta);
     

  double mediafinale=0;   // magnetizzazione complessiva
  double magnxfinale[Lx/2+1]={0};   // magnetizz in funzione della distanza dal bordo
  double corr2finale[16][16][16][16]={0};
  double meansqcorr2[16][16][16][16]={0};
  double Ecorr2finale[16][16][16][16]={0};
  double meansqEcorr2[16][16][16][16]={0};
  double meansq[Lx/2+1]={0};           // mediaquadra di diversi giri MC per avere varianza
  double exfinale[Lx/2+1]={0};
  double meansqE[Lx/2+1]={0};

  

  if (world_rank < world_size-1) 
  {
        rng=-time(NULL)*world_rank;     // settare entrambi in modo che dipenda da che core si sta usando
        srand (time(NULL)*world_rank);

        // correlazioni dell'energia hanno variabili con stessi nomi di corr magn ma con E davanti
        
       double mediatemp=0;
       double magnxtemp[Lx/2+1]={0};
       double extemp[Lx]={0};

       double corr2temp[16][16][16][16]={0};
       double Ecorr2temp[16][16][16][16]={0};

       
          //filldafile(); 
          randfill();
          
          //ofstream nomequalsiasi; nomequalsiasi.open("dati.txt");
          for(int mosse=0;mosse<nmontecarlo; mosse++) //for(int mosse=0;mosse<nmontecarlo;mosse++)
      {
        SweepEWolff();
         if(mosse>steptermalizz) 
          {
           for(int ii=0;ii<Lx/2+1;ii++) 
              {
                extemp[ii]+=(energiax(ii)+energiax(Lx-1-ii))/2;
                magnxtemp[ii]+=(magnx(ii)+magnx(Lx-1-ii))/2;
              }
             corr(corr2temp);
             Ecorr(Ecorr2temp); 
             mediatemp+= mean();                 
          }
                   
      }       // finito un giro di montecarlo
         for(int ii=0;ii<Lx/2+1;ii++) 
              {
                magnxtemp[ii]/=(stepdopotermalizz-1);
                extemp[ii]/=(stepdopotermalizz-1);
              }
          mediatemp/=(stepdopotermalizz-1);
              // manda vettore magnetizzazione, vettore energie, matriciona correlazioni, usando tag diversi
        MPI_Send(
              /* data         = */ magnxtemp,//&mediatemp, 
              /* count        = */ Lx/2+1, 
              /* datatype     = */ MPI_DOUBLE, 
              /* destination  = */ world_size-1, 
              /* tag          = */ world_rank, 
              /* communicator = */ MPI_COMM_WORLD
                );
          MPI_Send(
              /* data         = */ extemp,//&mediatemp, 
              /* count        = */ Lx/2+1, 
              /* datatype     = */ MPI_DOUBLE, 
              /* destination  = */ world_size-1, 
              /* tag          = */ 5000+world_rank,   // questi numeri non importano, basta che siano gli stessi in MPI_Receive, e non si usi due volte lo stesso per due send diversi
              /* communicator = */ MPI_COMM_WORLD
                );
        MPI_Send(
              /* data         = */ corr2temp,//&mediatemp, 
              /* count        = */ 65536,    // e' 16^4 
              /* datatype     = */ MPI_DOUBLE, 
              /* destination  = */ world_size-1, 
              /* tag          = */ 1000+world_rank, 
              /* communicator = */ MPI_COMM_WORLD
                );
        MPI_Send(
              /* data         = */ Ecorr2temp,//&mediatemp, 
              /* count        = */ 65536,    // e' 16^4 
              /* datatype     = */ MPI_DOUBLE, 
              /* destination  = */ world_size-1, 
              /* tag          = */ 9999+world_rank, 
              /* communicator = */ MPI_COMM_WORLD
                );
  }
    else if (world_rank == world_size-1) 
    { cout<< "dimensione del sistema e'" <<Lx<<"x"<<Ly<<"x"<<Lz<<endl
        <<"numero di sweep e' "<<nmontecarlo<<endl
        <<"numero di core e'"<<world_size<<endl;

      double ricevi[Lx/2+1];
      double ricevi2[Lx/2+1];
     
     //double cosotemp=0;
    for(int i=0;i<world_size-1;i++)             
    { 
      MPI_Recv(                     // ultimo core riceve magnetizzazioni
      /* data         = */ ricevi,//&cosotemp, 
      /* count        = */ Lx/2+1, 
      /* datatype     = */ MPI_DOUBLE, 
      /* source       = */ i, 
      /* tag          = */ i,           
      /* communicator = */ MPI_COMM_WORLD, 
      /* status       = */ MPI_STATUS_IGNORE);
      MPI_Recv(                     // ultimo core riceve energie
      /* data         = */ ricevi2,//&cosotemp, 
      /* count        = */ Lx/2+1, 
      /* datatype     = */ MPI_DOUBLE, 
      /* source       = */ i, 
      /* tag          = */ 5000+i, 
      /* communicator = */ MPI_COMM_WORLD, 
      /* status       = */ MPI_STATUS_IGNORE);
      
      for(int si=0;si<Lx/2+1;si++)                // ultimo core aggiunge magnetizzazioni ed energie alla media totale
      {
        magnxfinale[si]+=ricevi[si];
        exfinale[si]+=ricevi2[si];
        meansq[si]+=pow(ricevi[si],2);
        meansqE[si]+=pow(ricevi2[si],2);
      }
      
      double ricevicorr2[16][16][16][16]={0};
      double riceviEcorr2[16][16][16][16]={0};
      
      MPI_Recv(
      /* data         = */ ricevicorr2,//&cosotemp, 
      /* count        = */ 65536, 
      /* datatype     = */ MPI_DOUBLE, 
      /* source       = */ i, 
      /* tag          = */ 1000+i, 
      /* communicator = */ MPI_COMM_WORLD, 
      /* status       = */ MPI_STATUS_IGNORE);
      MPI_Recv(
      /* data         = */ riceviEcorr2,//&cosotemp, 
      /* count        = */ 65536, 
      /* datatype     = */ MPI_DOUBLE, 
      /* source       = */ i, 
      /* tag          = */ 9999+i, 
      /* communicator = */ MPI_COMM_WORLD, 
      /* status       = */ MPI_STATUS_IGNORE);
      double unacorr2divisa;
      double unaEcorr2divisa;
      for (int x = 1; x <= 15; x++) // ricevicorr contiene somme, per avere le medie bisogna dividere per numero giri dopo termalizz e anche per ampiezza di intervalli su cui si variano y0 e z0, che e' 5*(Lx-1) per ognuno
      {
        for (int x1 =x; x1 <= 16-x; x1++)
        {
          for (int deltay = 0; deltay < 15; deltay++)
          {
            for (int deltaz = deltay; deltaz < 15; deltaz++)
            { 
              unacorr2divisa=ricevicorr2[x][x1][deltay][deltaz]/( (stepdopotermalizz-1)*25*(Lx-1)*(Lx-1) );
              corr2finale[x][x1][deltay][deltaz]+=unacorr2divisa; 
              meansqcorr2[x][x1][deltay][deltaz]+=pow(unacorr2divisa,2);//

              unaEcorr2divisa=riceviEcorr2[x][x1][deltay][deltaz]/( (stepdopotermalizz-1)*25*(Lx-1)*(Lx-1) );
              Ecorr2finale[x][x1][deltay][deltaz]+=unaEcorr2divisa; 
              meansqEcorr2[x][x1][deltay][deltaz]+=pow(unaEcorr2divisa,2);//


            }
          }
        }
      }
    
    
    } //fine ciclo su vari core che mandano roba all'ultimo
      for(int ii=0;ii<Lx/2+1;ii++) 
          {
            cout<<magnxfinale[ii]/(world_size-1)<<" , ";;
          }
      cout<<endl<<"errori sulle medie"<<endl;
     for(int u=0;u<Lx/2+1;u++)
    {
      cout<<sqrt(meansq[u]/(world_size-1)-pow(magnxfinale[u]/(world_size-1),2))/sqrt(world_size-1)<<" , "; //divisione ulteriore per worldsize-1 è perché errore sulla media è varianza/sqrt n
    }
    cout<<"energie: "<<endl;
    for(int ii=0;ii<Lx/2+1;ii++) 
              {
                cout<<exfinale[ii]/(world_size-1)<<" , ";
              }
          cout<<endl<<"errori sulle medie delle energie"<<endl;
         for(int u=0;u<Lx/2+1;u++)
        {
          cout<<sqrt(meansqE[u]/(world_size-1)-pow(exfinale[u]/(world_size-1),2))/sqrt(world_size-1)<<" , "; 
        }
        
    tabella(corr2finale,meansqcorr2,(world_size-1),world_size,0);         //tabella prende corr finali e le salva in file con quel nome
    tabella(Ecorr2finale,meansqEcorr2,(world_size-1),world_size,1); 
    
    cout<<endl<< "tempo trascorso e'"<<time(NULL)-tiniz<<endl;

  }

  MPI_Finalize();

      //mediafinale+=mediatemp;
      
      
    /*for(int u=1;u<Lx/2+1;u++)
    {
      cout<<magnxfinale[u]<<" , ";
       //<<sqrt(meansq[u]/numgiri-pow(magnxfinale[u]/numgiri,2))<<endl; 
    }*/
    /*cout<<endl<<"le deviazioni standard sono"<<endl;
     for(int u=1;u<Lx/2+1;u++)
    {
      cout<<sqrt(meansq[u]/numgiri-pow(magnxfinale[u]/numgiri,2))<<" , "; 
    }*/
    
    
  /*ofstream nomequalsiasi;
  nomequalsiasi.open("dati24con1E4step2giri.txt"); 
  for(int u=1;u<Lx/2+1;u++)
    {
      nomequalsiasi<<magnxfinale[u]/numgiri<<" , "
       <<sqrt(meansq[u]/numgiri-pow(magnxfinale[u]/numgiri,2))<<endl; 
    }
    
           
  nomequalsiasi.close();
   */
  return 0;
}
    


