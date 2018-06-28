/*
  In deze versie bepaald iedere processor de priemgetallen tot en met sqrt(N)
  Iedere processor maakt daarna een cyclische verdeling van de getallen sqrt(N) t/m N
  
  Vervolgens wegstrepen wat hij heeft.
*/

#include "mcbsp.hpp"     //bsp library
#include "bspedupack.c"  //vecalloc, sizes
#include <iostream>
#include <math.h>    //oa floor
#include <sstream>   //voor goede output naar het scherm
#include <stdlib.h>  //atoi: string to int

using namespace std;

/*
  Commandline argumenten:
  1. Het aantal threads
  2. De bovengrens voor het grootste priemgetal
  3. Menuitem:
     1. Zeven en tellen
     2. Tweeling priemen
     3. Het vermoeden van Goldbach
*/

unsigned int P;
unsigned long long N;


unsigned long long nloc(int p, int s, unsigned long long M){
  /* calculate the number of local components for processor s
     for a vector of length M cyclically distributed over p processors*/
  return (M + p - s - 1) / p;
}


/*void crossout(int *numbers, int *start, unsigned long long *nl){
  //input is de vector met getallen die mijn processor heeft
  //en het aantal getallen
  //output is dezelfde lijst waarin niet priemen door 0 vervangen zijn
  cout << "hoi" << endl;
}
*/
void sieve(){
  bsp_begin(P);
  /* Vanaf nu parallel */
  
  /* SUPERSTEP 1*/  
  unsigned long long p = bsp_nprocs();   // p = number of processors
  unsigned long long s = bsp_pid();      // s = processor number
  double time0 = 0, time1 = 0; // voor timings
  unsigned long long i = 0, j = 0, k = 0;
  unsigned long long count = 0, count2 = 0; /*2 zit niet in de lijst, bij het schrijven van de streepstartwaarde tellen.*/
  ostringstream oss1, oss2, oss3, oss4;
  string part1, part2, part3, part4, out;
  /* alleen de int's tm sqrt(N) proberen */
  unsigned long long M = ceil((sqrt(N) - 1) / 2);
  int *A, *start;

  A = vecalloci(M + 1);
  start = vecalloci(M + 1); //om de startwaarde voor het strepen te bewaren

  /* De vector met integers initialiseren
     0 = false
     1 = true
  */

  A[0] = 0;
  for(i = 1; i < M + 1; i++){
    A[i] = 1;
  }

  /*bsp_sync(); //aanzetten voor iets lager tijdsgebruik, zie ook regel met time1=bsp_time();*/
  if(s == 0){
    time0 = bsp_time();
  }
  //strepen, tellen en bijhouden welke startwaarde er voor het strepen geld
  start[count] = 4;
  count++;
  for(i = 1; i < M + 1; i++){
    if(A[i] == 1){
      for(j = 2*i*i + 2 * i; j < M + 1; j += 2*i + 1){
        A[j] = 0;
      }
      start[count]=(2*i+1)*(2*i+1);
      count++;
    }
  }



//if(s==0){cout << ceil(sqrt(N)) << ", " << A[M] << endl;
//cout << count << endl;}



  /*bsp_sync();*/
  if(s == 0){
    time1 = bsp_time();
  }

  /* vector maken met daarin alle priemen op volgorde */
  int *primes;
  primes = vecalloci(count);
  primes[0] = 2;
  i = 1;
  for(j = 0; j < M + 1; j++){
    if(A[j] == 1){
      primes[i] = 2*j+1;
//if(s==0){cout<<primes[i] << endl;}
      i++;
    }
  }

  /* de priemen zijn bekend, dus verwijder geheugen van vector A */
  vecfreei(A);

  /* Maak de vector integers aan die op mijn processor te vinden zijn,
     hiervoor is ceil(sqrt(N)) als laatste onderzocht
     ceil(sqrt(N)) + 1 is dus eerstvolgende te onderzoek getal
  */
  int *numbers;
  unsigned long long L = ceil(sqrt(N));
  unsigned long long K; //hulpvariabele
  M = N - L; //aantal totaal dat nog gedaan moet worden
  unsigned long long nl = nloc(p, s, M);    //aantal voor mijzelf
  numbers = vecalloci(nl);




//if(s==0){cout << "L " << L << ", K " << K << ", M " << M << ", nl " << nl << endl;}

  /* Het echte aanmaken van de vector elementen */
  K = L + 1 + s; //K is het eerste getal op processor s

//if(s==0){cout << "L " << L << ", K " << K << ", M " << M << ", nl " << nl << endl;}


  for(i = 0; i < nl; i++){

    /*if i * p + s + K == 1 mod 2 invullen?*/


    numbers[i] = i * p + K;  //uitrekenen welke getallen ik heb, K is mijn eerste nummer
 //if(i==nl-1){cout<< numbers[i] << endl;}
    //if(numbers[i] == 75){
//cout << "probleem s,i " << s << ", " << i << ", ";}
  }

  /* De vector met integers na ceil(sqrt(N)) is gemaakt. Nu per processor verder strepen
     volgens de cyclische verdeling.

     Let op: we zijn hier geswitched naar alle getallen ipv de onevens!
  */






//============================================

  /* Correctie bepalen; is natuurlijk voor iedere processor gelijk */
  int tmp = 0, correctie = -1;
  while(correctie == -1){
    if((K + tmp) % p == s){
      correctie = tmp;
    }
    else{
      tmp++;
    }
  }



//if(s==0){cout << "L " << L << ", K " << K << ", M " << M << ", nl " << nl << endl;}


  //K = L + s;


//if(s==0){cout << "L " << L << ", K " << K << ", M " << M << ", nl " << nl << endl;}


//if(s==0){cout << "cor " << correctie << ", K " << K << endl;}




  for(i = 0; i < count; i++){
    /* vanaf start[i] met stappen primes[i] uitvegen als ik het getal heb */
//if(s==0){cout << i << endl;}
    for(j = start[i]; j <= N; j = j + primes[i]){
//if(s==0){cout << j << endl;}

//if(j==75){cout << "i: " << i  << " j: " << j << endl;}
//if(s==0){cout << "i,j " << i << ", " << j << endl;}


if(j > L){//if(s==3){cout << "ja" << endl;}
      /* het uit te vegen getal is j en is nog niet aan bod geweest */
      if((j + correctie) % p == s){// || (j + correctie) % p == s - p){
        /* processor s heeft het getal j */
        /* en het staat op plek */
        k = (j - K) / p;
  //if(s==0){cout << "k: " << k << endl; }
        /* daar uitvegen, dat wil zeggen, op waarde op 0 zetten als dat
           nog niet gedaan was, en tellen*/
        if(numbers[k] != 0){
          numbers[k] = 0;
         //if(j==75){cout << "s,i,j,k " << s << ", " << i << ", " << j << ", " << k << ", ";}
          count2++;
          //if(s==0){  cout << "count++ " << count2 << endl;}
        }


}
      }
    }
  }
  //cout << "count2 " << count2 << endl;
  /*crossout(numbers, start, &nl);*/


    /*oss1 << s;
    oss2 << p;
    oss3 << nl-count2;
    oss4 << count;
    part1 = oss1.str();
    part2 = oss2.str();
    part3 = oss3.str();
    part4 = oss4.str();
    out = "Thread " + part1 + " out of " + part2 + " has found " + part3 + " prime numbers between sqrt(N) and N, also before, eq sqrtN " + part4 + " seconds!\n";
    cout << out;
    // Set strings to be empty and clear all states
    oss1.str("");
    oss2.str("");
    oss3.str("");
    oss4.str("");
    oss1.clear();
    oss2.clear();
    oss3.clear();
    oss4.clear();

bsp_sync();
int count3 = 0;
if (s==1){

  for(i=0; i<nl; i++){
    if(numbers[i] != 0){
      cout << i << " "<< numbers[i] << ", " << endl;
      count3++;
    }
  }
}
cout << endl;
bsp_sync();
if (s==3){

  for(i=0; i<nl; i++){
    if(numbers[i] != 0){
      cout << i << " " << numbers[i] << ", " << endl;
      count3++;
    }
  }
}*/
for(int t= 0; t<p; t++){
  bsp_sync();
if (s==t){
cout << "tm wortel(N) heb ik (" << t << " uit " << p << " processoren) gevonden " << count << " en in het cyclisch verdeelde volgende stuk tm N nog eens " << nl-count2 << " priemgetallen." << endl;}
}
  bsp_end();
}



int main(int arg_count, char ** arg_vector){
  /* sequential menu part */
  int keuze; //menukeuze
  int i = 1; // afhandeling inputvector, i=0 is de naam van het programma

  /* Aantal te gebruiken processoren instellen */
  if(arg_count > i){
    P = atoi(arg_vector[i]);
    i++;
  }
  else{
    cout << "How many threads do you want to start? There are " <<  bsp_nprocs() << " cores available." << endl;
    cin >> P;
  }
  if(P == 0 || P > bsp_nprocs()){
    cout << "Cannot start " << P << " threads." << endl;
    return 1;
  }

  /* Bovengrens N vaststellen */
  if(arg_count > i){
    N = atoi(arg_vector[i]);
    i++;
  }
  else{
    cout << "What is the positive upperbound for finding primes?" << endl;
    cin >> N;
  }
  if(N <= 1){
    cout << "There are no primes smaller than " << N << "." << endl;
    return 1;
  }

  /* Menukeuze keuze vaststellen */
  if(arg_count > i){
    keuze = atoi(arg_vector[i]);
    i++;
  }
  else{
    cout << "What do you want to do? Choose an option from:" << endl
         << "[1] Sieve and count" << endl 
         << "[2] Search for twinprimes" << endl
         << "[3] Check the Goldbach conjecture" << endl
         << "[0] Stoppen" << endl;
    cin >> keuze;
  }

  /* menukeuze uitvoeren, eventueel als dit niet werkt, dan if s==0 onderscheiden! */
  switch(keuze){
    case 0:
      cout << "The program will end now." << endl;
      return 1;
      break;
    case 1:
      bsp_init(sieve, arg_count, arg_vector);
      sieve();
      break;
    default:
      cout << "That was an invalid option." << endl
           << "The program will end now." << endl;
      return 1;
      break;
  }
  cout << endl;
  return 0;
}






