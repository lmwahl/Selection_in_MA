#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#define MAXTYPES 1501   // maximum number of fitness classes to keep track of
#define MAXGENS  30     // maximum number of generations of growth
#define MAXMUTES 10001  // maximum number of entries in the DFE

/* This simulates an experimental mutation accumulation protocol */

/* running from the command line, a single argument (positive integer)
   can be included to set the seed for the random number generator */

/* This version uses a gamma-distributed underlying DFE, the version using
   double exponentials is a slight modification, in ma_fitness_exp.c */

void main(int argc, char **argv)

{

  FILE *fpout,*fps,*fpfit;
  char filename[100];
  long Nvisible;   // minimum number of indls in a visible colony
  int i, igen, ngens, iline, nlines, itransfer, ntransfers;
  int initialanc, threshold, muthreshold, nNtypes, mtwt;
  long Nend, nlucky, ilucky, Ntmp, new, nmu;
  double oneplus2s;
  double Wsave, r, dfe[3][MAXMUTES], fmax, gammax, gamshift, delW;
  int gamk, muput, jsave = 0;
  double mutotal;
  int nsinglets;   // total number of single-mutation entries in DFE
  /*  Ntypes stores the population size of each fitness class
   The zeroth class has relative fitness in the bin 0 to delW,
   the first class has relative fitness in the bin delW to 2*delW, etc.
   There are two columns in Ntypes.
   In the zeroth column we count members of the colony that differ
   from the founding indl of that colony (immediate "ancestor")
   by exactly one mutation (or do not differ at all).
   In the first column we keep any members of the colony that
   differ from the founding indl of that colony by 2 or more mutations. 
   nextNtypes is a working copy of Ntypes that we'll need later.  */
  long **Ntypes = (long **)malloc(2 * sizeof(long *));
  for (i=0;i<=1;i++) Ntypes[i] = (long *)malloc(MAXTYPES * sizeof(long));
  long **nextNtypes = (long **)malloc(2 * sizeof(long *));
  for (i=0;i<=1;i++) nextNtypes[i] = (long *)malloc(MAXTYPES * sizeof(long));
   
  int ianc;  // fitness class of the ancestor of the currently growing colony
  long N[MAXGENS];  // total population size of the colony
  int j, iold, imu, indmu;
  float poidev(float,long *), ran1(long *), gamdev(int, long *);
  long seed=-1,seedsave;
  /* We will also sample the underlying DFE nsdist times and save the
     the frequency of each fitness class in sdist.  When the number of
     mutations is larger than muthreshold, we will use this distribution
     to distribute them among classes in the next generation, rather than
     choosing each mutation's fitness value stochastically  */
  double *sdist = (double *)malloc(MAXTYPES * sizeof(double));
  int nsdist = 100000, isdist;

  // set the random number generator seed according to the command line
  // argument, or else it stays at -1
  if (argc>1) seed = -((long)(atof(argv[1])));
  seedsave = seed;
  Nvisible = 1e5;  // colonies chosen below this threshold trigger a warning
  threshold = 500;  // lineage above this size grows deterministically
  muthreshold = 1000;  // more than 1000 mutations from one lineage?  then distribute s values deterministically
  ngens = 20;   // generations of growth = number of doublings btwn transfers
  fprintf(stdout,"%d generations between transfers\n",ngens);
  // do everything for the wt parameters and then repeat for the mutant
  for (mtwt = 0;mtwt<=0;mtwt++) {
   seed = seedsave;  // use the same random seed for mt and wt, shouldn't matter
   nlines = 1000;  // number of replicate lines to simulate
   if (mtwt) {  //mutant
    fprintf(stdout,"Using mutant params\n");
    ntransfers = 30;   // transfers per line, 30 mutn and 300 wt
    mutotal= 0.00043;   // need genome-wide mutation rate, 4.3e-4 vs. 0.0041
   }
   else {   //wildtype
    fprintf(stdout,"Using wildtype params\n");
    ntransfers = 5;   // transfers per line, 30 mutn and 300 wt
    mutotal= 0.04;   // need genome-wide mutation rate, 4.3e-4 vs. 0.0041}
   }   
   gamk = 10;   // shape parameter for gamma distn of underlying DFE
   gammax = 0.5;  // maximum s value in Gamma distribution
                  // gamdev(gamk)/gamk gives gamma with mean 1, shape gamk
                  // gammax*gamdev(gamk)/gamk gives gamma withe mean gammax
		  // we then subtract this from 1+gammax to give reflected
		  // Gamma distn with maximum value 1+gammax and
		  // mean 1 (for relative fitness)
   gamshift = 0.2; // We can then shift the mean to below 1 by this amount
   // resulting DFE: reflected Gamma with shape gamk, mean 1-gamshift
   // and maximum value 1+gammax-gamshift
  fmax = 0.5;   // highest relative fitness in the fitness bins will be 1+fmax
  //  This is different from max value of DFE because beneficial mutations
  //  could possibly accrue over the course of the experiment.
  nNtypes = MAXTYPES-1;
  delW = (double)((1.0+fmax)/nNtypes);   // width of fitness bins
  initialanc = (int)(1.0/delW);    // initial ancestor has relative fitness 1
                                   // so this gives its bin index
  nsinglets = 0;  // total number of entries in dfe (so far)

  // check that we can open the output file properly before doing the work
  sprintf(filename, "madata/bd_%dtrans_%dlines_%dk_%1.2fM_%1.5fmu_%dgens_%dseed.out",ntransfers,nlines,gamk,gammax-gamshift,mutotal,ngens,-seed);
  fpout = fopen(filename,"w");
  if (fpout==NULL) { fprintf(stdout,"Unable to open output file\n"); exit(1);}
  sprintf(filename, "madata/bd_%dtrans_%dlines_%dk_%1.2fM_%1.5fmu_%dgens_%dseed_fitness.out",ntransfers,nlines,gamk,gammax-gamshift,mutotal,ngens,-seed);
  fpfit = fopen(filename,"w");
  if (fpfit==NULL) { fprintf(stdout,"Unable to open output file\n"); exit(1);}

  /* set up sdistn so that when there are many mutants,
     we can distribute them deterministically  */
  for (i=0;i<nNtypes;i++) sdist[i] = 0;
  for (isdist = 1; isdist<=nsdist; isdist++) {
    r = 1+ gammax -gamshift- gammax*gamdev(gamk,&seed)/gamk;
    indmu = (int)(r/delW);
    if (indmu>=0 && indmu<nNtypes) sdist[ (int)(r/delW)]++;
  }
  for (i=0;i<nNtypes;i++) sdist[i]/=nsdist;
  // save the underlying DFE to an output file
  fps = fopen("madata/sdist.out","w");
  for (i=0;i<nNtypes;i++) fprintf(fps,"%f\n",sdist[i]);
  fclose(fps);
  

  // loop for each line
  for (iline=1;iline<=nlines;iline++) { 
    // initially start with one indl with fitness 1
    for (i=0;i<nNtypes;i++) {
     Ntypes[0][i] = 0; Ntypes[1][i]=0;}
     ianc = initialanc;
     Ntypes[0][ianc] = 1;
     Wsave = ianc*delW;   // save the ancestor fitness here
// loop to do ngens of growth and then sampling, repeated ntransfers times
    for (itransfer=1;itransfer<=ntransfers;itransfer++) {
      for (igen=1;igen<=ngens;igen++) {  
        // calculate the current size of the colony
        N[igen] = 0;
	for (i=0;i<nNtypes;i++) N[igen] += Ntypes[0][i] + Ntypes[1][i];
        if (N[igen]==0) { // this can't happen with the pure birth process
	  // but we'll leave this check in the code anyway.
	  // If extinction happens, back up and regrow
	  // from the ancestor and choose again...
          fprintf(stdout,"Population went extinct at generation %d\n",igen);
          N[igen] = 1;
	  ianc = (int)(Wsave/delW);
          Ntypes[0][ianc] = 1;
	  igen = igen-1;
	}

/* Poisson reproduction, pure birth process
Note:  2*(1+s) = 2 + 2s = 1 + (1+2s) = parents + offspring
Each fitness class gets copied into the next genn (new = Ntypes).
We allow mutations to happen to these "parent" individuals as well.
The parents have a Poisson-generated number of offspring with mean Ntypes*(1+2s)
If a fitness class has over 500, it grows deterministically (Ntypes*2*(1+s)).
Mutation can independently occur with each birth.
If the number of mutants > 1000, they are added deterministically. */
        for (iold=0;iold<nNtypes;iold++)
	  nextNtypes[0][iold] = nextNtypes[1][iold] = 0;
	for (iold = 0;iold<nNtypes;iold++) {
	  for (j=0;j<=1;j++) {
          if (Ntypes[j][iold]>0) {
	    if (Ntypes[j][iold]>threshold) {
	      new = 2*Ntypes[j][iold]*delW*iold;
	      if (new<0) {  // this happens if we exceed MAX_LONG
	        fprintf(stderr,"neg new %d %d %d\n\n",iline,itransfer,igen);
	        exit(1);
	      }
	    }
	    else {
	      new = Ntypes[j][iold];
	      oneplus2s = 1.0 + 2.0*(delW*iold -1.0);
	      if (oneplus2s>0)
	        new += (int)poidev((float)(Ntypes[j][iold]*oneplus2s),&seed);
	    }
	    // Now that we know how many new individuals to add, let them mutate
	    if (new>0)
	      nmu  = (long)poidev((float)(new*mutotal),&seed);
	    else nmu = 0;
	    //  The next 2 lines tell us to put any new indls into column 1
	    // of nextNtypes, unless they are new mutations of the ancestor,
	    // in which case they will be stored in column 0
	    muput = 1;
	    if (iold==ianc) muput = 0;
	    if (nmu>muthreshold) {  //determ mutns
	      for (i=0;i<nNtypes;i++) {  // i gives the fitness class of DFE
		// newfitness = old fitness * r, which in other words is:
		//       newW = (iold*delW)*(delW*i);  
		// Then newindex = newfitness/delW, which in other words is:
	        //          indmu = (int)(newW/delW);
		// We combine the two steps above as follows:
		indmu = (int)(iold*delW*i);
	        if (indmu>=0 && indmu<nNtypes)
		  nextNtypes[muput][indmu] += (int)(sdist[i]*nmu);
	      }
	    }
	    else {
	    for (imu = 1; imu<=nmu;imu++) {  //stochastic mutn
	      r = 1.0 + gammax - gamshift - gammax*gamdev(gamk,&seed)/gamk;
		//newW = iold*delW*r;  indmu = newW/delW;
		indmu = (int)(iold*r);
	        if (indmu>=0 && indmu<nNtypes)
		  nextNtypes[muput][indmu] += 1;
	    }}
	     //for new-nmu births, mutation did not occur
	        nextNtypes[j][iold] += new-nmu;
	  } // if non-zero Ntypes
	  } // loop on i and j, 1 step and multistep mutations
	} // loop on iold

	// new becomes the old
	for (i=0;i<nNtypes;i++) {
	  Ntypes[0][i] = nextNtypes[0][i];
	  Ntypes[1][i] = nextNtypes[1][i];
	}
	} // loop on igen, growth is finished

      /* now we randomly sample a single individual.
This is the individual that lands closest to the marked spot on the
next plate.  We will only simulate the growth process for this colony. */

      // Nend gives total size of colony at end of growth
      Nend =0;
      for (i=0;i<nNtypes;i++) Nend += Ntypes[0][i] + Ntypes[1][i];
      if (ngens>20 && Nend<Nvisible)
        fprintf(stderr,"Warning: we chose a colony with only %d indls\n",Nend);
      nlucky = (long)(ran1(&seed)*Nend);
      Ntmp = 0; jsave = 1;
      ilucky = ianc;
      for (i=0;i<nNtypes;i++) {
        for (j=0;j<=1;j++) {
	  Ntmp += Ntypes[j][i];  
	  if (Ntmp>nlucky) {
	    ilucky = i;
	    jsave = j;
	    i = nNtypes+1; j = 2;  // jump out of both loops
	  }
        }
      }
	/* add this mutation to the list if previous ancestor 
	   and this lucky differ by  a single mutation  */
      if ((ilucky != ianc) && (jsave == 0)) { 
        nsinglets++;
	if (nsinglets>MAXMUTES) fprintf(stderr,"Error: MAXMUTES exceeded\n"); 
        dfe[0][nsinglets] = ilucky*delW/Wsave;
	dfe[1][nsinglets] = (double)itransfer;
	dfe[2][nsinglets] = ilucky*delW;
      }
      // now lucky becomes the new ancestor for the next growth
      ianc = (int)ilucky;
      for (i=0;i<nNtypes;i++) { Ntypes[0][i] =Ntypes[1][i]= 0; }   // this type is the new ancestor
      Ntypes[0][ianc] = 1;
      Wsave = ianc*delW;
      fprintf(fpfit," %f ",Wsave);
    }  // loop on itransfers, do the next growth/sampling now
    fprintf(fpfit,"\n");
   fprintf(stdout,"%d ",iline);

  } // at this point, this line is done, return to do a new line

  // at this point, all the lines are done.  Save results.
  fprintf(stdout,"\n");
  for (i=1;i<=nsinglets;i++)
    fprintf(fpout,"%f %f %f\n",dfe[0][i],dfe[1][i],dfe[2][i]);
  fclose(fpout);
  }  // return to the top and repeat for the mutant or wildtype
}
