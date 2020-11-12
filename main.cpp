/**
 * Point packing code.  Best from each run is in best.pcid including the number
 * of parameter settings and the actual settings.
 */

/**
 * ITERATIVE VERSION
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cmath>

using namespace std;

#include "stat.h"

//random number seed
#define RNS 91207819

//code dimensions
#define dim 9

//minimum distance and its square
#define mind 0.525
double md2=mind*mind;

//run controls - number of runs, mating events per run and rerporting interval
#define runs 30
#define mevs 1000000
#define RI 10000
//running report
#define verbose 1

//population specifiers, size, random material rate, upper bound on code size
#define popsize 4000
#define RMR 20
#define MAX 1000
//conway crossover buffer size
#define CCB 2020

//Conway variables
double pot[CCB][dim];  //buffer for crossover
double res[CCB][dim];  //result for crossover
int Cz;  //resulting crossover size

//population variables
double pop[popsize][MAX][dim];  //holds population members
int csize[popsize];             //current code sizes
double fit[popsize];            //holds fitness value
int dx[popsize];                //sorting index

void initalg(); //initialize the algorithm
int okay(double *pt,double code[MAX][dim],int z);  //is a new point okay? 
void rpoint(double *pt); //generate a random point
void cpoint(double *a,double *b);  //copy a point
void spoint(double *a,double *b);  //copy a point
void create(double gene[MAX][dim],int &z);  //make an initial population member
void initpop(); //intialize a population
void ccross(int a,int b);      //perform Conway crossover of popmembers a and b
void matingevent();             //perform a mating event
void report(ostream &aus);      //make a statistical report
void reportbest(ostream &aus);  //report the best structure

int main(){//main program

    int run,mev;         //run and mating event look counters
    fstream stat,crit;   //statistics and structure reporting
    char fn[60];         //file name construction buffer

    initalg();      //initialize the algorithm
    crit.open("best.pcld",ios::out);  //open the good structure channel
    for(run=0;run<runs;run++){//loop over runs
        sprintf(fn,"run%03d.dat",run);  //create statistics filename
        stat.open(fn,ios::out);  //open the statistics channel
        initpop();      //initialize a population
        report(stat);  //make an intial report
        for(mev=0;mev<mevs;mev++){//loop over mating events
            matingevent();  //do a mating event
            if((mev+1)%RI==0){//make a report if its time
                if(verbose==1)cout << run << " " << (mev+1)/RI << " ";
                report(stat);  //to the report
            }
        }
        stat.close();  //close the statistics channel
        reportbest(crit);  //report the best structure
    }

    crit.close(); //close the reporting channel

    return(0);  //keep the system happy

}

void initalg(){//initialize the algorithm

    srand48(RNS); //seed the random number generator

}

double D2(double *a,double *b){//find the squared distance between points

    double delta,accu;    //distance scratch variables
    int i;                //loop index

    accu=0.0;  //zero the accumulator
    for(i=0;i<dim;i++){//loop over the coordinates
        delta=a[i]-b[i];  //compute the difference
        accu+=(delta*delta); //add up the contribution for that coordinate
    }

    return(accu);  //report the value

}

int okay(double *pt,double code[MAX][dim],int z){//is a new point okay? 

    int i,j;              //loop index variables

    for(i=0;i<z;i++){//loop over poitns currently in code
        if(D2(code[i],pt)<md2)return(0); //too close
    }

    return(1); //sruvived all checks

}

void rpoint(double *pt){//generate a random point

    double iv[dim]; //interval to be partitioned
    double sw;      //swap variable
    int fl;         //flag
    int i;          //loop index

    for(int i=0;i<dim;i++)iv[i]=drand48();  //partition the interval
    do {
        fl=0;  //no flag
        for(i=0;i<dim-1;i++)//sort!
            if(iv[i]>iv[i+1]){sw=iv[i];iv[i]=iv[i+1];iv[i+1]=sw;fl=1;}
    }while(fl==1);//until in order
    pt[0]=iv[0];
    for(i=1;i<dim-1;i++)pt[i]=iv[i]-iv[i-1];
    pt[dim-1]=1-iv[dim-2];

}

void cpoint(double *a,double *b){//copy a point

    for(int i=0;i<dim;i++)a[i]=b[i];  //copy the values

}

void spoint(double *a,double *b){//copy a point

    double sw; //swap variables

    for(int i=0;i<dim;i++){sw=a[i];a[i]=b[i];b[i]=sw;}  //swap the values

}

void create(double gene[MAX][dim],int &z){//make an initial population member

    double pt[dim];  //tentative new point
    int i,j;         //loop index variables

    z=0;
    for(i=0;i<MAX;i++){//loop over samples
        rpoint(pt); //make a random point
        if(okay(pt,gene,z)){//found an acceptable point
            for(j=0;j<dim;j++)gene[z][j]=pt[j];  //transfer the point
            z++;  //register its existence
        }
    }
}

void initpop(){//intialize a population

    int i;  //loop index variable

    for(i=0;i<popsize;i++){//loop over population
        create(pop[i],csize[i]);    //go make the population member
        fit[i]=((double)csize[i]);  //convert the fitness to a real number
        dx[i]=i;                    //refresh the sorting index
        cout << fit[i] << " ";
    }
    cout << endl;
}

void ccross(int a,int b){//perform Conway crossover of popmembers a and b

    int i,j,k;       //index variables
    double pt[dim];  //new point

    //Assemble the buffer
    k=0;  //number of things in the pot buffer
    for(i=0;i<csize[a];i++){//copy the first guy
        cpoint(pot[k++],pop[a][i]);
    }
    for(i=0;i<csize[b];i++){//copy the first guy
        cpoint(pot[k++],pop[b][i]);
    }
    for(i=0;i<RMR;i++){//add the random material
        rpoint(pt);
        cpoint(pot[k++],pt);
    }

    //Shuffle the buffer
    for(i=0;i<k;i++){//loop over points in buffer
        j=lrand48()%k; //select random coordinate
        if(i!=j)spoint(pot[j],pot[i]);  //swap them if they are different
    }

    //Ready for lexicode algorithm
    Cz=0;  //zero the size
    for(i=0;i<k;i++)if(okay(pot[i],res,Cz)){//if the point fits
            cpoint(res[Cz],pot[i]);  //copy it in
            Cz++;  //and record its existence
        }

}

void matingevent(){//perform a mating event

    int i;   //index variable

    tselect(fit,dx,3,popsize);  //select the three guys
    ccross(dx[1],dx[2]);        //perform conway crossover
    if(Cz>=csize[dx[0]]){//result at least as good?
        csize[dx[0]]=Cz;  //update size
        for(i=0;i<Cz;i++)cpoint(pop[dx[0]][i],res[i]); //copy the code
        fit[dx[0]]=((double)Cz); //update fitness
    }

}

void report(ostream &aus){//make a statistical report

    dset D;

    D.add(fit,popsize);  //build data set
    if(verbose)cout << D.Rmu() << " " << D.RCI95() << " "
                    << D.Rsg() << " " << D.Rmax() << endl;
    aus << D.Rmu() << " " << D.RCI95() << " "                //print out
        << D.Rsg() << " " << D.Rmax() << endl;

}

void reportbest(ostream &aus){//report the best structure

    int i,j,k,b;
    double q;

    b=0;  //initialize best pointer
    for(i=1;i<popsize;i++)if(fit[i]>fit[b])b=i;
    aus << fit[b] << " -fitness" << endl;
    k=0;
    for(i=0;i<csize[b];i++){//loop over vectors
        aus << pop[b][i][0];
        q=pop[b][i][0];
        for(j=1;j<dim;j++){aus << " " << pop[b][i][j];q+=pop[b][i][j];}
        aus << endl;

    }


}
