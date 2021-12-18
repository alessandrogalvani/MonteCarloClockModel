
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


typedef std::mt19937 rng_type;
std::uniform_int_distribution<rng_type::result_type> randombit(1,2);
std::uniform_int_distribution<rng_type::result_type> randomspin(1,p);
std::uniform_int_distribution<rng_type::result_type> randomx(1,Lx-2);
std::uniform_int_distribution<rng_type::result_type> randomy(0,Ly-1);

rng_type rng2;
/*int randint(int max){return 1+(int) (double(max)*rand()/(RAND_MAX+1.0)); }
int randombit(int r){return randint(2);}
int randomspin(int r){return randint(p);}
int randomx(int r){return randint(Lx-2);}
int randomy(int r){return randint(Ly)-1;}
*/



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
int idan(int n){return n/(Ly2);}
int jdan(int n){return n%(Ly2)/Ly;}
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
int vicini[N][6];
void settavicini()
{
  for(int n=0;n<N;n++)
    {
      int i=idan(n);
      int j=jdan(n);
      int k=kdan(n);
      vicini[n][0]= (i==0? -1: ndaijk(i-1,j,k) ) ;
      vicini[n][1]= (i==Lx-1? -1: ndaijk(i+1,j,k) );
      vicini[n][2]= ndaijk(i, (j+1)%Ly , k );
      vicini[n][3]= ndaijk(i, (j+Ly-1)% Ly, k);
      vicini[n][4]= ndaijk(i, j, (k+1)%Ly );
      vicini[n][5]= ndaijk(i, j, (k+Ly-1)%Ly );

    }


}

double energialocale(int n)
{
  double E=0;
  int sito=m[n];
  int i=idan(n);

  if(sito!=0)
  {
    E+=Dfrattobeta;
    for(int w=0;w<6;w++)
      if(vicini[sito][w]!=-1)
        { 
          E+= v[abs(sito- m[vicini[sito][w] ] ) ];
        }
  }
  return E;
}
double magnx(int i)    
{
  double mediax=0; double mediay=0;

       
  for(int u=0;u<Ly2;u++)
  {
  int mtemp=m[i*Ly2+u];
    if(mtemp!=0) 
      {
        mediax+=v[mtemp-1];
        mediay+=vseni[mtemp-1];
      }
  }
    
  mediax/=(Ly2); mediay/=(Ly2);
  return mediax*v[megaspin-1]+mediay*vseni[megaspin-1];

}


double energiax(int i)    
{
  double media=0;
       
  for(int u=0;u<Ly2;u++)
      media+=energialocale(i*Ly2+u);

  return media/(Ly*Lz);

}

void spinflip1(int n)      // svuota o riempie sito
{
  int angolonuovo;
  double expmenodeltaE;
  int mtemp=m[n];
  if(mtemp!=0)
    {
      angolonuovo= 0;

      expmenodeltaE=expmD;
      for(int w=0;w<6;w++)
      if(vicini[n][w]!=-1)
        { 
          if(m[vicini[n][w] ]!=0) expmenodeltaE *= pesim[abs(mtemp- m[vicini[n][w] ] ) ];
        }

    }

  if(mtemp==0)
    {
      angolonuovo= randomspin(rng2);

      expmenodeltaE=expD;             

      for(int w=0;w<6;w++)
      if(vicini[n][w]!=-1)
        { 
          if(m[vicini[n][w] ]!=0) expmenodeltaE *= pesi[abs(m[vicini[n][w] ] - angolonuovo ) ];
        }
    }

  if(randomreal(& rng)<expmenodeltaE) {m[n]=angolonuovo; }
  
}
void spinflip2(int n)      // ergodico
{
  int angolonuovo;
  double expmenodeltaE=1;
  int mtemp=m[n];
  if(randombit(rng2)==1) angolonuovo=0;    // SE SI FA ISING SI DEVE PROPORRE 0 CON PROBABILITA' 1/3, NON 1/2
  else angolonuovo=randomspin(rng2);

            
  if(angolonuovo!=0)
    {       // copincollato da spinflip1, potrebbe essere sbagliato

      expmenodeltaE*=expD; 

      for(int w=0;w<6;w++)
      if(vicini[n][w]!=-1)
        { 
          if(m[vicini[n][w] ]!=0) expmenodeltaE *= pesi[abs(m[vicini[n][w] ] - angolonuovo ) ];
        }

    }

  if(mtemp!=0)
    {   
      expmenodeltaE*=expmD; 

      for(int w=0;w<6;w++)
      if(vicini[n][w]!=-1)
        { 
          if(m[vicini[n][w] ]!=0) expmenodeltaE *= pesim[abs(mtemp- m[vicini[n][w] ] ) ];
        }

    }




      
  if(randomreal(& rng)<expmenodeltaE) {m[n]=angolonuovo;}
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
  for (int n = 0; n < N; n++)
    cluster[n] = 0;


  int i0 = randomx(rng2);
  int j0 = randomy(rng2);
  int k0 = randomy(rng2);
  int n0= ndaijk(i0,j0,k0);
  int r= randomspin(rng2)-1;
 // per evitare di fare roba ricorsiva si fa un while su siti da fare. ogni terna che indica un sito viene convertita in un numero solo, che è scritto in una base in cui una cifra è in base Lx e le altre due in base Ly. poi le tre funzioni ida,jdan,kdan riconvertono questo numero in i,j,k


  vector<int> sitidafare;
  vector<int> siticluster;

  if(m[n0]!=0)
    {
     cluster[n0]=1;
     sitidafare.push_back(n0);
     siticluster.push_back(n0);
    }
  int aa=0;
  while(!sitidafare.empty())
  {   //prima si determina chi va in cluster, poi si flippa tutti


    int sito=sitidafare[0];
    int i=idan(sito);
    if( (i==0||i==Lx-1)&&aa==0 )
      { aa=1;
        megaspin=flipangolo(megaspin,r);
       
          for (int u = 0; u < Ly2; u++)
          {
            if(!cluster[u])//if(i!=0||jj!=j||kk!=k)
            {
              cluster[u]=1;
              sitidafare.push_back(u); 
              siticluster.push_back(u);
            }
            if(!cluster[(Lx-1)*Ly2+u])//if(i!=Lx-1||jj!=j||kk!=k)
            {
              cluster[(Lx-1)*Ly2+u]=1;
              sitidafare.push_back( (Lx-1)*Ly2+u ); 
              siticluster.push_back( (Lx-1)*Ly2+u );
            }
          }
         
        
      }
      
      int vicinotemp;
      for(int w=0;w<6;w++)
      {vicinotemp=vicini[sito][w];
        if(vicinotemp!=-1)
        { 
          if(!cluster[vicinotemp ] && m[ vicinotemp ] !=0 ) 
            if (randomreal(& rng) < pbond[abs(2*(m[sito]-1)-r)][abs(2*(m[vicinotemp ] -1)-r) ] )
          { 
          cluster[vicinotemp ]=1;
          sitidafare.push_back( vicinotemp );
          siticluster.push_back( vicinotemp );
          }

        }
      }
      
    sitidafare.erase(sitidafare.begin());
  }// fine del while
//cout<<siticluster.size()<<endl;
  for(vector<int>::iterator q=siticluster.begin();q!=siticluster.end();q++)
  { 
    int n=*q;
   //cout<<is<<" "<<js<<" "<<ks<<endl;

    m[n]=flipangolo(m[n],r);
  }//cout<<"taglia "<<siticluster.size()<<endl;
  
}



void SweepEWolff()
{
  for(int n=Ly2;n<N-Ly2;n++)
      {
        spinflip2(n);
      }
      
    for(int rr=0;rr<6;rr++)
      {
        for(int n=Ly2;n<N-Ly2;n++)
          {
            spinflip1(n);
          }
       
          wolff();
       
      }
}
//ARRIVATI FINO QUA


//int rifletti(int x){return 16-x;}

double c2(int n, int n1)
{
  if( m[n]==0 || m[n1]==0 ) return 0;
  else  return v[ abs( m[n]-m[n1] ) ] ;
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
                    correlazioni[x][x1][deltay][deltaz]+=( c2(ndaijk( fr*x,j0,k0),ndaijk(fr*x1,j0+fr*deltay,k0+fr*deltaz) ) + 
                      c2( ndaijk(fr*(16-x),j0,k0 ) , ndaijk(fr*(16-x1),j0+fr*deltay,k0+fr*deltaz) )+
                      c2( ndaijk(fr*x,j0,k0) , ndaijk( fr*x1,j0+fr*deltaz,k0+fr*deltay) ) + 
                      c2( ndaijk( fr*(16-x),j0,k0) , ndaijk( fr*(16-x1),j0+fr*deltaz,k0+fr*deltay) ) )/4;
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
                    Ecorrelazioni[x][x1][deltay][deltaz]+=( energialocale(ndaijk(fr*x,j0,k0) )*energialocale( ndaijk(fr*x1,j0+fr*deltay,k0+fr*deltaz) ) + 
                      energialocale( ndaijk(fr*(16-x),j0,k0) )*energialocale( ndaijk(fr*(16-x1),j0+fr*deltay,k0+fr*deltaz) )+
                      energialocale( ndaijk(fr*x,j0,k0) )*energialocale( ndaijk(fr*x1,j0+fr*deltaz,k0+fr*deltay) ) + 
                      energialocale( ndaijk(fr*(16-x),j0,k0) )*energialocale( ndaijk(fr*(16-x1),j0+fr*deltaz,k0+fr*deltay) ) )/4;
                  }
                }
              }
            }
        }
    }
}




double Ecorrvicine(double Ecorrelazioni[16][16][16][16])    // ci sono 7672 correlazioni indi, per ora non si trasla parallelamente su y e z ma poi va fatto
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
                    Ecorrelazioni[x][x1][deltay][deltaz]+=(
                    energialocale( ndaijk(frSUX*x,j0,k0) )*energialocale( ndaijk(frSUX*x1,j0+fr*deltay,k0+fr*deltaz) ) + 
                    energialocale( ndaijk(Lx-1-(x*frSUX),j0,k0) )*energialocale( ndaijk(Lx-1-(x1*frSUX),j0+fr*deltay,k0+fr*deltaz) )+
                    energialocale( ndaijk(frSUX*x,j0,k0) )*energialocale( ndaijk(frSUX*x1,j0+fr*deltaz,k0+fr*deltay) ) + 
                    energialocale( ndaijk(Lx-1-(x*frSUX),j0,k0) )*energialocale( ndaijk(Lx-1-(x1*frSUX),j0+fr*deltaz,k0+fr*deltay) ) )/4;
                  }
                }
              }
            }
        }
    }
}


void tabella(double correlazioni[16][16][16][16], double meansqcorr[16][16][16][16], int N, int worldsize, bool indice, string nomem, string nomee)
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
  if(indice==0)  nomequalsiasi.open(nomem);
  else nomequalsiasi.open(nomee);
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
