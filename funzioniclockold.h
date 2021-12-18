
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
double randomreal(long *idum) // LA SUA VARIABILE DEV'ESSERE & coso, dove coso definito all'inizio coso=-time(NULL) SERVE IL MENO
{ 
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;
  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX; 
  else {return temp;}
}


int randint(int max){return 1+(int) (double(max)*rand()/(RAND_MAX+1.0)); }
int randombit(){return randint(2);}
int randomspin(){return randint(p);}
int randomx(){return randint(Lx-2);}
int randomy(){return randint(Ly)-1;}




//std::uniform_int_distribution<int> cento(1,100); 

int next(int i)         //se si vuole avere pbc anche su x serve Lx=Ly oppure definire nuovi next/previous
{ 
  if (i==Ly-1) return 0;
  else return i+1;
}
int previous(int i)
{
  if (i==0) return Ly-1;
  else return i-1;
}
int ndaijk(int i, int j, int k){return k+j*Ly+i*Ly*Ly;}
int idan(int n){return n/(Ly*Ly);}
int jdan(int n){return n%(Ly*Ly)/Ly;}
int kdan(int n){return n%Ly;}

void sitodopo(int &i,int &j, int &k)

{
  if (i<Lx-2) 
            
    i++;
        
  else 
    {
      i=1;
      if(j<Ly-1)
	j++;
      else
	{   j=0;
	  if(k<Lz-1) k++;
	  else k=0;
	}
    }
        
}// si usa per far crescere indici ijk in spinflip1 e spinflip2

void printpiano( int k)
    {
        for(int i=0;i<Lx;i++)
        {
            for(int j=0;j<Ly;j++)
            {
                cout << m[i][j][k];
            }
            cout << endl ;
        }
        cout << endl;
    }
void printcluster(int k)
    {
        for(int i=0;i<Lx;i++)
        {
            for(int j=0;j<Ly;j++)
            {
                cout <<" "<< cluster[i][j][k];
            }
            cout << endl ;
        }
        cout << endl;
    }//printa un piano di cluster di wolff



double mean() 
{
  double mediax=0; double mediay=0;

  for(int i=0;i<Lx;i++)
    {
      for(int j=0;j<Ly;j++)
      {
        for(int k=0;k<Lz;k++)
          {
            if(m[i][j][k]!=0) {mediax+=v[abs(m[i][j][k]-1)]; mediay+=vseni[abs(m[i][j][k]-1)];}
          }
      }
    }
  mediax/=N; mediay/=N;

  return mediax*v[megaspin-1]+mediay*vseni[megaspin-1];    
}
double energialocale(int i, int j, int k)
{
  double E=0;
  int sito=m[i][j][k];
  if(sito!=0)
  { E+=Dfrattobeta;
    if(i>0)if(m[i-1][j][k]!=0)    E+=v[ abs(sito-m[i-1][j][k]) ]/2. ;
    if(i<Lx-1)if(m[i+1][j][k]!=0) E+=v[ abs(sito-m[i+1][j][k]) ]/2 ;
    if(m[i][next(j)][k]!=0)       E+=v[ abs(sito-m[i][next(j)][k]) ]/2 ;
    if(m[i][previous(j)][k]!=0)   E+=v[ abs(sito-m[i][previous(j)][k]) ]/2 ;
    if(m[i][j][next(k)]!=0)       E+=v[ abs(sito-m[i][j][next(k)]) ]/2 ;
    if(m[i][j][previous(k)]!=0)   E+=v[ abs(sito-m[i][j][previous(k)]) ]/2 ;
    
  }
  return E;
}
double magnx(int i)    
{
  double mediax=0; double mediay=0;

       
  for(int j=0;j<Ly;j++)
    {
      for(int k=0;k<Lz;k++)
  {
    if(m[i][j][k]!=0) mediax+=v[m[i][j][k]-1];
    if(m[i][j][k]!=0) mediay+=vseni[m[i][j][k]-1];
  }
    }
  mediax/=(Ly*Lz); mediay/=(Ly*Lz);
  return mediax*v[megaspin-1]+mediay*vseni[megaspin-1];

}


double energiax(int i)    
{
  double media=0;
       
  for(int j=0;j<Ly;j++)
      for(int k=0;k<Lz;k++)
      media+=energialocale(i,j,k);

  return media/(Ly*Lz);

}

void spinflip1(int i, int j, int k)      // svuota o riempie sito
{
  int angolonuovo;
  double expmenodeltaE;

  if(m[i][j][k]!=0)
    {
      angolonuovo= 0;

      expmenodeltaE=expmD;
      if(i<Lx-1)  if (m[i+1][j][k]!=0)    expmenodeltaE*= pesim[ abs(m[i+1][j][k]-m[i][j][k]) ];
      if(i>0)     if(m[i-1][j][k]!=0)     expmenodeltaE*= pesim[ abs(m[i-1][j][k]-m[i][j][k]) ];
      if(m[i][next(j)][k]!=0) expmenodeltaE*= pesim[ abs(m[i][next(j)][k]-m[i][j][k]) ]; 
      if(m[i][previous(j)][k]!=0) expmenodeltaE*= pesim[ abs(m[i][previous(j)][k]-m[i][j][k]) ];
      if(m[i][j][next(k)]!=0) expmenodeltaE*= pesim[ abs(m[i][j][next(k)]-m[i][j][k]) ]; 
      if(m[i][j][previous(k)]!=0) expmenodeltaE*= pesim[ abs(m[i][j][previous(k)]-m[i][j][k]) ];
            
    }

  if(m[i][j][k]==0)
    {
      angolonuovo= randomspin();

      expmenodeltaE=expD;             
      if(i<Lx-1)  if(m[i+1][j][k]!=0) expmenodeltaE*= pesi[ abs(m[i+1][j][k]-angolonuovo)];       
      if(i>0)     if(m[i-1][j][k]!=0) expmenodeltaE*= pesi[ abs(m[i-1][j][k]-angolonuovo)];
      if(m[i][next(j)][k]!=0) expmenodeltaE*= pesi[ abs(m[i][next(j)][k]-angolonuovo)];
      if(m[i][previous(j)][k]!=0) expmenodeltaE*= pesi[ abs(m[i][previous(j)][k]-angolonuovo)];
      if(m[i][j][next(k)]!=0) expmenodeltaE*= pesi[ abs(m[i][j][next(k)]-angolonuovo)];
      if(m[i][j][previous(k)]!=0) expmenodeltaE*= pesi[ abs(m[i][j][previous(k)]-angolonuovo)];
    }

  if(randomreal(& rng)<expmenodeltaE) {m[i][j][k]=angolonuovo; }
  
}
void spinflip2(int i, int j, int k)      // ergodico
{
  int angolonuovo;
  double expmenodeltaE=1;
        
  if(randombit()==1) angolonuovo=0;    // SE SI FA ISING SI DEVE PROPORRE 0 CON PROBABILITA' 1/3, NON 1/2
  else angolonuovo=randomspin();

            
  if(angolonuovo!=0)
    {       // copincollato da spinflip1, potrebbe essere sbagliato
      expmenodeltaE*=expD;  
      if(i<Lx-1)  if (m[i+1][j][k]!=0) expmenodeltaE*= pesi[ abs(m[i+1][j][k]-angolonuovo)];      
      if(i>0)     if(m[i-1][j][k]!=0)  expmenodeltaE*= pesi[ abs(m[i-1][j][k]-angolonuovo)];
      if(m[i][next(j)][k]!=0) expmenodeltaE*= pesi[ abs(m[i][next(j)][k]-angolonuovo)];
      if(m[i][previous(j)][k]!=0) expmenodeltaE*= pesi[ abs(m[i][previous(j)][k]-angolonuovo)];
      if(m[i][j][next(k)]!=0) expmenodeltaE*= pesi[ abs(m[i][j][next(k)]-angolonuovo)];
      if(m[i][j][previous(k)]!=0) expmenodeltaE*= pesi[ abs(m[i][j][previous(k)]-angolonuovo)];
               
    }

  if(m[i][j][k]!=0)
    {   
      expmenodeltaE*=expmD; 
      if(i<Lx-1)  if (m[i+1][j][k]!=0) expmenodeltaE*= pesim[ abs(m[i+1][j][k]-m[i][j][k]) ];
      if(i>0)     if(m[i-1][j][k]!=0)  expmenodeltaE*= pesim[ abs(m[i-1][j][k]-m[i][j][k]) ];
      if(m[i][next(j)][k]!=0) expmenodeltaE*= pesim[ abs(m[i][next(j)][k]-m[i][j][k]) ]; 
      if(m[i][previous(j)][k]!=0) expmenodeltaE*= pesim[ abs(m[i][previous(j)][k]-m[i][j][k]) ];
      if(m[i][j][next(k)]!=0) expmenodeltaE*= pesim[ abs(m[i][j][next(k)]-m[i][j][k]) ]; 
      if(m[i][j][previous(k)]!=0) expmenodeltaE*= pesim[ abs(m[i][j][previous(k)]-m[i][j][k]) ];
           
    }




      
  if(randomreal(& rng)<expmenodeltaE) {m[i][j][k]=angolonuovo;}
  // printpiano(m,k); }//cout << "spin flippato in "<< angolonuovo << endl;}
            
            
  //    cout << "si è flippato" <<peso<< endl;
}

int flipangolo(int a,int r)    // prende numero da 1 a p e dà numero da 1 a p
{
 return ( (r-a+1+p/2+p)%p ) +1;
}

void settapesi(double be)
{
  // vettore con possibili valori coseni
  for(int jj=0;jj<p;jj++)
    {
      v[jj]=cos(2.*pi*double(jj)/p); if( abs(v[jj])<1E-7) v[jj]=0;
      vseni[jj]=sin(2.*pi*double(jj)/p);
    }

  // vettore con possibili pesi boltzmann

  for(int jj=0;jj<p;jj++)
    {
      pesi[jj]=exp(be*v[jj]);
      pesim[jj]=1/pesi[jj];     
    }

  for(int jj=0;jj<2*p;jj++)
    {
      for(int kk=0;kk<2*p;kk++)
        {
    pbond[jj][kk]=1-exp(-2* be *cos(pi*double(jj)/p)* cos(pi*double(kk)/p) );  
        if(pbond[jj][kk]<0) pbond[jj][kk]=0;
        }
    }


}




void wolff() // matrice Lx*Ly*Lz cluster serve per controllare che sito non sia gia' stato preso. vettore cluster serve per flippare tutti assieme alla fine. vettore sitidafare e' quello che si svuota gradualmente
{ 
  for (int i = 0; i < Lx; i++)
    for (int j = 0; j < Ly; j++)
      for (int k = 0; k < Lz; k++)
      cluster[i][j][k] = 0;


  int i0 = randomx();
  int j0 = randomy();
  int k0 = randomy();

  int r= randomspin()-1;
 // per evitare di fare roba ricorsiva si fa un while su siti da fare. ogni terna che indica un sito viene convertita in un numero solo, che è scritto in una base in cui una cifra è in base Lx e le altre due in base Ly. poi le tre funzioni ida,jdan,kdan riconvertono questo numero in i,j,k


  vector<int> sitidafare;
  vector<int> siticluster;

  if(m[i0][j0][k0]!=0)
    {
     cluster[i0][j0][k0]=1;
     sitidafare.push_back(ndaijk(i0, j0,k0));
     siticluster.push_back(ndaijk(i0, j0,k0));
    }
  int aa=0;
  while(!sitidafare.empty())
  {   //prima si determina chi va in cluster, poi si flippa tutti
    int i=idan(sitidafare[0]);
    int j=jdan(sitidafare[0]);
    int k=kdan(sitidafare[0]);
    if( (i==0||i==Lx-1)&&aa==0 )
      { aa=1;
        megaspin=flipangolo(megaspin,r);
        for (int jj = 0; jj < Ly; jj++)
        {
          for (int kk = 0; kk < Lz; kk++)
          {
            if(!cluster[0][jj][kk])//if(i!=0||jj!=j||kk!=k)
            {
              cluster[0][jj][kk]=1;
              sitidafare.push_back(ndaijk(0, jj,kk)); 
              siticluster.push_back(ndaijk(0, jj,kk));
            }
            if(!cluster[Lx-1][jj][kk])//if(i!=Lx-1||jj!=j||kk!=k)
            {
              cluster[Lx-1][jj][kk]=1;
              sitidafare.push_back(ndaijk(Lx-1, jj,kk)); 
              siticluster.push_back(ndaijk(Lx-1, jj,kk));
            }
          }
         }
        
      }
      

    
      
        if(i>0&&!cluster[i-1][j][k]&&m[i-1][j][k]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i-1][j][k]-1)-r) ] )
          { 
          cluster[i-1][j][k]=1;
          sitidafare.push_back(ndaijk(i-1, j,k));
          siticluster.push_back(ndaijk(i-1, j,k));
          }
        if(i<Lx-1&&!cluster[i+1][j][k]&&m[i+1][j][k]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i+1][j][k]-1)-r) ] )
          {
          cluster[i+1][j][k]=1;
          //if(i!=Lx-1||cluster[0][j][k]==0)  // questo e' per evitare che si raggiunga bordo contemporaneamente da due
          sitidafare.push_back(ndaijk(i+1, j,k));
          siticluster.push_back(ndaijk(i+1, j,k));
          }
        if(i>0&&i<Lx-1&&!cluster[i][previous(j)][k]&&m[i][previous(j)][k]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][previous(j)][k]-1)-r) ] )
          {
          cluster[i][previous(j)][k]=1;
          sitidafare.push_back(ndaijk(i, previous(j),k)); 
          siticluster.push_back(ndaijk(i, previous(j),k));
          }

        if(i>0&&i<Lx-1&&!cluster[i][next(j)][k]&&m[i][next(j)][k]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][next(j)][k]-1)-r) ] )
          {
          cluster[i][next(j)][k]=1;
          sitidafare.push_back(ndaijk(i, next(j),k)); 
          siticluster.push_back(ndaijk(i, next(j),k));
          }        

        if(i>0&&i<Lx-1&&!cluster[i][j][next(k)]&&m[i][j][next(k)]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][j][next(k)]-1)-r) ] )
          {
          cluster[i][j][next(k)]=1;
          sitidafare.push_back(ndaijk(i, j , next(k) ) ); 
          siticluster.push_back(ndaijk(i, j , next(k) ) );
          }  

        if(i>0&&i<Lx-1&&!cluster[i][j][previous(k)]&&m[i][j][previous(k)]!=0)
          if (randomreal(& rng) < pbond[abs(2*(m[i][j][k]-1)-r)][abs(2*(m[i][j][previous(k)]-1)-r) ] )
          {
          cluster[i][j][previous(k)]=1;
          sitidafare.push_back(ndaijk(i, j,previous(k) )); 
          siticluster.push_back(ndaijk(i, j,previous(k) ));
          }
          
    
      
    sitidafare.erase(sitidafare.begin());
  }// fine del while
//cout<<siticluster.size()<<endl;
  for(vector<int>::iterator q=siticluster.begin();q!=siticluster.end();q++)
  { 
    int is=idan(*q);
    int js=jdan(*q);
    int ks=kdan(*q);
   //cout<<is<<" "<<js<<" "<<ks<<endl;
    m[is][js][ks]=flipangolo(m[is][js][ks],r);
  }//cout<<"taglia "<<siticluster.size()<<endl;
  
}



void SweepEWolff()
{
  int i=1; int j=0; int k=0;
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
       // wolff();
       // cout<<megaspin<<m[0][1][2]<<m[Lx-1][4][6]<<endl;
        //printcluster(1);
      }
}

//int rifletti(int x){return 16-x;}

double c2(int x, int y, int z, int x1, int y1, int z1)
{
  if( m[x][y][z]==0 || m[x1][y1][z1]==0 ) return 0;
  else  return v[ abs( m[x][y][z]-m[x1][y1][z1] ) ] ;
}

double corr(double correlazioni[16][16][16][16])    // ci sono 7672 correlazioni indi, per ora non si trasla parallelamente su y e z ma poi va fatto
{
  //int y0=0; int z0=0; // bisogna fare due for pure su questi
  
 //v[ abs( m[i][j][k]-m[][][]) ];
  for(int j0=0;j0<5*(Lx-1);j0++)    // j e k sono punti veri, x y e z sono da moltiplicare perche' vanno fino a 16
    {
      for(int k0=0;k0<5*(Lx-1);k0++)
        {
           for(int x=1;x<=15;x++)
            {
             for(int x1=x;x1<=16-x;x1++) // per avere piramide, x' ha due valori possibili in meno ogni volta che x cresce DA X A 16-X
              {
                for (int deltay = 0; deltay < 15; deltay++)
                {
                  for (int deltaz = deltay; deltaz < 15; deltaz++)
                  {
                    correlazioni[x][x1][deltay][deltaz]+=( c2(fr*x,j0,k0,fr*x1,j0+fr*deltay,k0+fr*deltaz) + c2(fr*(16-x),j0,k0,fr*(16-x1),j0+fr*deltay,k0+fr*deltaz)+
                      c2(fr*x,j0,k0,fr*x1,j0+fr*deltaz,k0+fr*deltay) + c2( fr*(16-x),j0,k0,fr*(16-x1),j0+fr*deltaz,k0+fr*deltay) )/4;
                  }
                }
              }
            }
        }
    }
}


double Ecorr(double Ecorrelazioni[16][16][16][16])    // ci sono 7672 correlazioni indi, per ora non si trasla parallelamente su y e z ma poi va fatto
{
  //int y0=0; int z0=0; // bisogna fare due for pure su questi
  
 //v[ abs( m[i][j][k]-m[][][]) ];
  for(int j0=0;j0<5*(Lx-1);j0++)    // j e k sono punti veri, x y e z sono da moltiplicare perche' vanno fino a 16
    {
      for(int k0=0;k0<5*(Lx-1);k0++)
        {
           for(int x=1;x<=15;x++)
            {
             for(int x1=x;x1<=16-x;x1++) // per avere piramide, x' ha due valori possibili in meno ogni volta che x cresce DA X A 16-X
              {
                for (int deltay = 0; deltay < 15; deltay++)
                {
                  for (int deltaz = deltay; deltaz < 15; deltaz++)
                  {
                    Ecorrelazioni[x][x1][deltay][deltaz]+=( energialocale(fr*x,j0,k0)*energialocale(fr*x1,j0+fr*deltay,k0+fr*deltaz) + 
                      energialocale(fr*(16-x),j0,k0)*energialocale(fr*(16-x1),j0+fr*deltay,k0+fr*deltaz)+
                      energialocale(fr*x,j0,k0)*energialocale(fr*x1,j0+fr*deltaz,k0+fr*deltay) + 
                      energialocale( fr*(16-x),j0,k0)*energialocale(fr*(16-x1),j0+fr*deltaz,k0+fr*deltay) )/4;
                  }
                }
              }
            }
        }
    }
}


void tabella(double correlazioni[16][16][16][16], double meansqcorr[16][16][16][16], int N, int worldsize, bool indice)
{cout<<endl;
  double corrp;
  /*string nomefile="corr2";
  int tempo=time(NULL);
  string txt=".txt";
  string str= to_string(tempo);
  string lunghezza= to_string(Lx-1);
  string ngiri=to_string(nmontecarlo);
  string nomecompleto=lunghezza+","+ngiri+nomefile+tempo+txt;*/
  ofstream nomequalsiasi;
  if(indice==0)  nomequalsiasi.open("correlazionim.txt");
  else nomequalsiasi.open("correlazionie.txt");
  //nomequalsiasi<< "dimensione del sistema e'" <<Lx<<"x"<<Ly<<"x"<<Lz<<endl
    //    <<"numero di sweep e' "<<nmontecarlo<<endl;
	//nomequalsiasi<< "x _ x' _ deltay delta z _ corr" <<endl;*/
  for (int x = 1; x <= 15; x++)
  {
    for (int x1 =x; x1 <= 16-x; x1++)
    {
      for (int deltay = 0; deltay < 15; deltay++)
      {
        for (int deltaz = deltay; deltaz < 15; deltaz++)
        { 
          corrp=correlazioni[x][x1][deltay][deltaz]/N;
          nomequalsiasi<< fr*x<<" "<<fr*x1<<" "<<fr*deltay<<" "<<fr*deltaz<<" "<< 
          corrp<<" "<<sqrt(   (meansqcorr[x][x1][deltay][deltaz]/N-pow(corrp,2))/(worldsize-1)    )<<endl; 
        }// errore alla fine e' ancora diviso per radice di worldsize-1 perchè errore sulla media è varianza/N, qui N=worldsize-1 ma potrebbero essere diversi, N includerebbe tipo numero di giri se non si fosse già diviso prima per quello
      }
    }
  }
  nomequalsiasi.close();
}










/*
 for (int x = 1; x <= 15; x++)
  {
    for (int x1 =x; x1 <= 16-x; x1++)
    {
      for (int deltay = 0; deltay < 15; deltay++)
      {
        for (int deltaz = deltay; deltaz < 15; deltaz++)
        { 
          corrp=correlazioni[x][x1][deltay][deltaz]/N;
          nomequalsiasi<<"{"<< fr*x<<", "<<fr*x1<<", "<<fr*deltay<<", "<<fr*deltaz<<", "<< 
          corrp<<", "<<sqrt(   (meansqcorr[x][x1][deltay][deltaz]/N-pow(corrp,2))/(worldsize-1)    )<<"},"<<endl; 
        }// errore alla fine e' ancora diviso per radice di worldsize-1 perchè errore sulla media è varianza/N
      }
    }
  }*/