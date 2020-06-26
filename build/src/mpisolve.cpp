/*
 
  Parallelized Finite Difference Time Domain Schrodinger Eq Solver
 
  mpisolve.c
 
  Copyright (c) Michael Strickland, Rafael L. Delgado
 
  GNU General Public License (GPLv3)
  See detailed text in license directory
 
*/

#include <iostream>
#include <iomanip>
#include <fstream>
#include <climits>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <cerrno>
#include <complex>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/types.h>
#include <unistd.h>

using namespace std;

#include "mpi.h"
#include "fftw3-mpi.h"
#include "mpisolve.h"
#include "grid.h"
#include "initialconditions.h"
#include "potential.h"
#include "outputroutines.h"
#include "paramreader.h"

// these global vars are initialized from parameters file
// defaults set here are overridden by that file
int    NUMX=20,NUM=20,STEPS=40000,UPDATE=100,SNAPUPDATE=1000;
int    POTENTIAL=0,INITCONDTYPE=0,INITCONDAXIS=0,INITSYMMETRY=0,SAVEWAVEFNCS=0,DUMPSNAPS=0,SNAPDUMP=0;
double A=0.05,EPS=0.001,SIG=0.06,MASS=1.0,SIGMA=0.223,TOLERANCE=-1;
int    KINTERM=0; //KINTERM GLOBAL VARIABLE (FLAG FOR THE USAGE OF RELATIVISTIC KINTETIC TERM)
char   *EXTPOT=NULL; //EXTERNAL POTENTIAL FILE
char   *DATAFOLD=NULL; //EXTERNAL DATA FOLDER
int    SAVEPOT=0; //SAVE POTENTIAL
int    SAVEDECAY=0; //SAVE DECAY TABLE
double POTCRITR=3., POTFLATR=-1.; //CRITIAL R
// mpi vars
int nodeID,tmpnodeID, numNodes;

// files
fstream debug_out;

// debug flag; options are DEBUG_{OFF,ON,FULL}
int debug = DEBUG_OFF;

// used for MPI non-blocking sends
double *leftSendBuffer,*rightSendBuffer;
MPI_Status leftSendStatus,rightSendStatus;
MPI_Status leftSendMessageStatus,rightSendMessageStatus;
MPI_Request leftSend,rightSend;
MPI_Request leftMessageSend,rightMessageSend;

// used for MPI non-blocking receives
double *leftReceiveBuffer,*rightReceiveBuffer;
MPI_Status leftReceiveStatus,rightReceiveStatus;
MPI_Status leftReceiveMessageStatus,rightReceiveMessageStatus;
MPI_Request leftReceive,rightReceive;
MPI_Request leftMessageReceive,rightMessageReceive;

// used for worker only intercommunication
MPI_Comm workers_comm;
MPI_Group workers_group;

// used for fftw3 library
struct fftw_par_t {
  fftw_plan plan_fftw, plan_ifftw;
  fftw_complex *in, *out;
  ptrdiff_t local_n0, local_0_start;
  ptrdiff_t alloc_local;
} fftw_par;

inline void dcomp_to_fftw(const dcomp &in, fftw_complex *out){(*out)[0] = in.real(); (*out)[1] = in.imag(); }
inline void fftw_to_dcomp(const fftw_complex &in, dcomp *out){*out = dcomp(in[0], in[1]); }
inline double pow2(double a){return a*a; }
inline dcomp pow2(dcomp a){return a*a; } 

// variables which will be loaded and reduced across nodes
dcomp energy=0;			// the local node energy
dcomp energyCollect=0;		// the total energy
dcomp normalization=0;		// the local node normalization squared
dcomp normalizationCollect=0;  	// the total normalization squared
dcomp vInfinity=0;		// the local node expectation value of v_infty
dcomp vInfinityCollect=0;      	// the total expectation value of v_infty
dcomp rRMS2=0;                 	// the local node <r^2>  
dcomp rRMS2Collect=0;		// the total <r^2>  
dcomp xAvg=0;			// the local <x>
dcomp xAvgCollect=0;		// the total <x>
dcomp yAvg=0;			// the local <y>
dcomp yAvgCollect=0;		// the total <y>
dcomp zAvg=0;			// the local <z>
dcomp zAvgCollect=0;		// the total <z>

// ground state energy and final time saved in global var after convergence
dcomp	EGrnd, timef;

// counter used for recording snapshots
int snapcnt = -1;

void computeObservables(dcomp*** wfnc);

int main( int argc, char *argv[] ) 
{ 
  int done=0,checksum=0;
  char message[1024]; 
  char fname[1024];
  int exclude[1];
  struct tms starttime,endtime;

  MPI_Init(&argc,&argv); 
  MPI_Comm_size(MPI_COMM_WORLD,&numNodes); 
  MPI_Comm_rank(MPI_COMM_WORLD,&tmpnodeID); 
  MPI_Status status;
  MPI_Group all_group;

  //For dealing with old Michael's setup, with a control node on nodeID=0
  nodeID=tmpnodeID+1;
	
  // setup group consisting of computational nodes only for internal communications
  MPI_Comm_group(MPI_COMM_WORLD, &all_group);
  //exclude[0]=0;
  //MPI_Group_excl(all_group, 1, exclude, &workers_group);
  //MPI_Comm_create(MPI_COMM_WORLD, workers_group, &workers_comm);
  //
  MPI_Comm dup_comm_world;
  MPI_Comm_dup(MPI_COMM_WORLD, &dup_comm_world);
  MPI_Comm_group(dup_comm_world, &workers_group);
  MPI_Comm_create(dup_comm_world, workers_group, &workers_comm);

  int workers_rank;
  MPI_Comm_rank(workers_comm, &workers_rank );

  MPI_Barrier(workers_comm);

  cout << "Process nodeID=" << nodeID << "\t workers_rank="<<workers_rank<<endl;
	
  if (debug) {
    sprintf(fname,"debug/debug_%d.txt",nodeID);
    debug_out.open(fname, ios::out);
    debug_out << "==> Node " << nodeID << " is ready" << endl; 
  }
	
  // node 1 is the master
  if (nodeID == 1) {
    times(&starttime); // load start time into starttime structure
    print_line();
    cout << "Parameters from file" << endl;
    print_line();
    readParametersFromFile((char *)"input/params.txt",1);
    if (argc>1) {
      print_line();
      cout << "Parameters from commandline" << endl;
      print_line();
      readParametersFromCommandLine(argc,argv,1);
    }
  }
  else {
    readParametersFromFile((char *)"input/params.txt",0);
    readParametersFromCommandLine(argc,argv,0);
  }

  if(!DATAFOLD) DATAFOLD=strdup("data");
	
  if (NUM%numNodes!=0) {
    if (nodeID==1)
      print_line();
    cout << "ERROR: Unable to partition lattice ... exiting" << endl;
    print_line();
    if (debug) {
      debug_out << "==> Goodbye from node " << nodeID << endl; 
      debug_out.close();
    }
    MPI_Finalize(); 
    exit(0);
  } else {
    NUMX = NUM/numNodes;	
  }
	
  if (nodeID == 1) {
    //	
    // master node
    //	
    for (int node=1;node<numNodes;node++) {
      sprintf(message,"Hello to node %d",node); 
      MPI_Send(message, strlen(message), MPI_CHAR, node, HELLO, MPI_COMM_WORLD); 
    }
		
    // master loops and waits for children to report that they are ready to start
    cout << "Receiving hello from nodes on master node..." << endl;
    checksum=0; 
    while( checksum < numNodes-1 ){
      MPI_Recv(&done, 1, MPI_INT, MPI_ANY_SOURCE, HELLO, MPI_COMM_WORLD, &status); 
      checksum += done;
      if (debug) debug_out << "Received: hello from computational node" << endl;
    }
		
    // cluster is ready
    if (debug) debug_out << "==> Cluster ready" << endl;
    cout << "==> Master node: cluster ready..." << endl;
		
    // Currently the master process
    // simply starts and turns into
    // a usual computational process.
    //
  }else{
    //	
    // computational node
    //	
    MPI_Recv(message, 20, MPI_CHAR, 0, HELLO, MPI_COMM_WORLD, &status); 
    if (debug == DEBUG_FULL) debug_out << "==> Received : "<< message << endl;
    done=1;
    MPI_Send( &done, 1, MPI_INT, 0, HELLO, MPI_COMM_WORLD );
  }
		
  double time;
      	
  // set initial conditions and get ready for computation
  // send message to master that we are starting
  MPI_Barrier(workers_comm);   
  solveInitialize();

  // the master said hello so now let's get to work
  MPI_Barrier(workers_comm);   
  solve();
      	
  // done with main computation, now do any analysis required
  MPI_Barrier(workers_comm);   
  solveFinalize();
      	
  cout.flush();
  MPI_Barrier(workers_comm);   

  // send message to master that we are done
  if (nodeID != 1){
    done = 1;
    MPI_Send( &done, 1, MPI_INT, 0, DONE, MPI_COMM_WORLD );

    if (debug) { 
      debug_out << "==> Goodbye from node " << nodeID << endl; 
      debug_out.close();
    }
  }else{
    // master loops and waits for children to report that they are done
    checksum=0;
    while( checksum < numNodes-1 ){
      MPI_Recv(&done, 1, MPI_INT, MPI_ANY_SOURCE, DONE, MPI_COMM_WORLD, &status); 
      checksum += done;
      if (debug) debug_out << "Received: checkout from computational node" << endl;
      sleep(1.0); // sleep 0.1 seconds between checks in order to reduce CPU usage of master
    } 
		
    times(&endtime); // load end time into endtime structure
    cout << "==> User time: " << (endtime.tms_utime - starttime.tms_utime)/((double)sysconf(_SC_CLK_TCK)) << " seconds"<< endl;
    cout << "==> System time: " << (endtime.tms_stime - starttime.tms_stime)/((double)sysconf(_SC_CLK_TCK)) << " seconds" << endl;
    print_line();
    cout << "Done." << endl;
    print_line();
  }

  MPI_Group_free(&workers_group);
  MPI_Finalize(); 

  return 0; 
} 

// solve initialize
void solveInitialize() {

  char label[64];
	
  // allocate memory
  allocateMemory();
	
  // load the potential
  if (nodeID==1) {
    print_line();
    cout << "==> Loading Potential Arrays" << endl;
    flush(cout);
  }
  loadPotentialArrays();
	
  if (nodeID==1) print_line();
	
  // set initial conditions
  setInitialConditions(nodeID);

  // initalize fftw3 library
  if (KINTERM!=0){
	  
    fftw_mpi_init();
	  
    fftw_par.alloc_local =
      fftw_mpi_local_size_3d(NUM, NUM, NUM, workers_comm,
			     &(fftw_par.local_n0), &(fftw_par.local_0_start));
	  
    fftw_par.in = fftw_alloc_complex(fftw_par.alloc_local);
    fftw_par.out = fftw_alloc_complex(fftw_par.alloc_local);

    fftw_par.plan_fftw =
      fftw_mpi_plan_dft_3d(NUM, NUM, NUM, fftw_par.in, fftw_par.out, workers_comm,
			   FFTW_FORWARD,  FFTW_MEASURE); //PATIENT | FFTW_DESTROY_INPUT);

    fftw_par.plan_ifftw =
      fftw_mpi_plan_dft_3d(NUM, NUM, NUM, fftw_par.in, fftw_par.out, workers_comm,
			   FFTW_BACKWARD, FFTW_MEASURE);//FFTW_PATIENT | FFTW_DESTROY_INPUT);

  }

  if(SAVEPOT){
    sprintf(label,"%d",nodeID); 
    outputPotential(label);
  } 
}

// reduce observables across nodes to first worker node
void computeObservables(dcomp*** wfnc) {
	
  // sum energy across nodes
  double energy_re=0.,energy_im=0.;
  double energy_re_collect=0.,energy_im_collect=0.;
  energy = wfncEnergy(wfnc);
  energy_re = real(energy);
  energy_im = imag(energy);
  MPI_Reduce(&energy_re,&energy_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  MPI_Reduce(&energy_im,&energy_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  energyCollect = dcomp(energy_re_collect,energy_im_collect);
	
  // sum normalization squared across nodes
  double normalization_re=0.,normalization_im=0.;
  double normalization_re_collect=0.,normalization_im_collect=0.;
  normalization = wfncNorm2(wfnc);
  normalization_re = real(normalization);
  normalization_im = imag(normalization);
  MPI_Reduce(&normalization_re,&normalization_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  MPI_Reduce(&normalization_im,&normalization_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  normalizationCollect = dcomp(normalization_re_collect,normalization_im_collect);
	
  // sum expectation across nodes
  double vInfinity_re=0.,vInfinity_im=0.;
  double vInfinity_re_collect=0.,vInfinity_im_collect=0.;
  vInfinity = vInfinityExpectationValue(wfnc);
  vInfinity_re = real(vInfinity);
  vInfinity_im = imag(vInfinity);	
  MPI_Reduce(&vInfinity_re,&vInfinity_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  MPI_Reduce(&vInfinity_im,&vInfinity_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  vInfinityCollect = dcomp(vInfinity_re_collect,vInfinity_im_collect);
	
  // sum r-squared across nodes
  double rRMS2_re=0.,rRMS2_im=0.;
  double rRMS2_re_collect=0.,rRMS2_im_collect=0.;
  rRMS2 = r2ExpectationValue(wfnc);
  rRMS2_re = real(rRMS2);
  rRMS2_im = imag(rRMS2);
  MPI_Reduce(&rRMS2_re,&rRMS2_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  MPI_Reduce(&rRMS2_im,&rRMS2_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  rRMS2Collect = dcomp(rRMS2_re_collect,rRMS2_im_collect);

  // sum x across nodes
  double x_re=0.,x_im=0.;
  double x_re_collect=0.,x_im_collect=0.;
  xAvg = xExpectationValue(wfnc);
  x_re = real(xAvg);
  x_im = imag(xAvg);
  MPI_Reduce(&x_re,&x_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  MPI_Reduce(&x_im,&x_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  xAvgCollect = A*dcomp(x_re_collect,x_im_collect);

  // sum y across nodes
  double y_re=0.,y_im=0.;
  double y_re_collect=0.,y_im_collect=0.;
  yAvg = yExpectationValue(wfnc);
  y_re = real(yAvg);
  y_im = imag(yAvg);
  MPI_Reduce(&y_re,&y_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  MPI_Reduce(&y_im,&y_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  yAvgCollect = A*dcomp(y_re_collect,y_im_collect);

  // sum z across nodes
  double z_re=0.,z_im=0.;
  double z_re_collect=0.,z_im_collect=0.;
  zAvg = zExpectationValue(wfnc);
  z_re = real(zAvg);
  z_im = imag(zAvg);
  MPI_Reduce(&z_re,&z_re_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  MPI_Reduce(&z_im,&z_im_collect,1,MPI_DOUBLE,MPI_SUM,0,workers_comm);
  zAvgCollect = A*dcomp(z_re_collect,z_im_collect);
	
}

// main computational solve routine
void solve() {
	
  dcomp energytot,lastenergy = 1.0e10;
  int step=0,done=1;
  bool keep=true;
  char label[64]; 
	
  leftSendBuffer = (double *)malloc( 2 * (NUM+2) * (NUM+2) * sizeof(double) );
  rightSendBuffer = (double *)malloc( 2 * (NUM+2) * (NUM+2) * sizeof(double) );
  leftReceiveBuffer = (double *)malloc( 2 * (NUM+2) * (NUM+2) * sizeof(double) );
  rightReceiveBuffer = (double *)malloc( 2 * (NUM+2) * (NUM+2) * sizeof(double) );

  // Check barrier...
  MPI_Barrier(workers_comm);   

  // evolve lattice in steps of UPDATE and output data along the way
  do {
    keep = (step<=STEPS-UPDATE);
		
    // sync boundaries
    if (KINTERM==0) syncBoundaries(w);

    // reduce observables across nodes
    computeObservables(w);

    // output 2d snapshots of the wavefunction for inspection
    // and check convergence of ground state energy
    if (step%SNAPUPDATE==0 || step%SNAPDUMP==0) {
      // broadcast observables
      MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
      MPI_Bcast(&energyCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
      // force symmetry
      symmetrizeWavefunction();
      // normalize wavefunction
      normalizeWavefunction(w);
      // save ground state energy and tau_f in global variables
      energytot = energyCollect/normalizationCollect;
      EGrnd = energytot;
      timef = step*EPS;

      if (nodeID==1) outputMeasurements(step*EPS, step);
 
      if (DUMPSNAPS>0 && step>2*UPDATE && step%SNAPDUMP==0){
        sprintf(label,"%d_0_%d",step,nodeID); 
        outputSnapshot(w,label);
 
        if (nodeID==1) outputSummaryData("Ground State", step*EPS, step);
        if (DUMPSNAPS>2) findExcitedStates(step*EPS, step);
      }

      // check convergence and break if tolerance is achieved
      // otherwise, record snapshot for use in excited state 
      // computation and keep going
      //
      if (abs(energytot-lastenergy)<TOLERANCE) {
	break;
      } else if (keep){
	lastenergy = energytot;
	// record and output snapshot
	snapcnt = (snapcnt+1)%2; // assume only two snapshots for now, so cycle
	recordSnapshot(w,snapcnt);
      }
    }
    if (keep){
      evolve(UPDATE);
      step += UPDATE;
    }
  } while (keep);
	
  computeObservables(w);
  MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);
  MPI_Bcast(&energyCollect, 1, MPI_DOUBLE_COMPLEX, 0, workers_comm);

  if (nodeID==1) outputSummaryData("Ground State", step*EPS, step);
	
  if (debug) {
    debug_out << "==> Unnormalized Energy : " << energy << endl;
    debug_out << "==> Normalization2 : " << normalization << endl;
  }
	
  free(rightSendBuffer);
  free(leftSendBuffer);
  free(rightReceiveBuffer);
  free(leftReceiveBuffer);

  return;
	
}

// solve finalize
void solveFinalize() {
	
  // this routine currently computes the first excited state energy and wavefunction
  findExcitedStates(real(timef), -1.);

  // this cleans the fftw3 staff
  if (KINTERM!=0){
    fftw_destroy_plan(fftw_par.plan_fftw);
    fftw_destroy_plan(fftw_par.plan_ifftw);
    fftw_free(fftw_par.in);
    fftw_free(fftw_par.out);
  }
}

// evolves solution nsteps
void evolve(int nsteps) {
	
  int leftTest,rightTest;

  for (int i=1;i<=nsteps;i++) {
		
    if (KINTERM==0){
      // receive boundary sync
      leftTest=1; 
      rightTest=1;
      if (i>1) {
        if (nodeID < numNodes) { 
          receiveRightBoundary(); 
          rightTest = 0;
        }
        if (nodeID-1 >= 1) {
          receiveLeftBoundary();
          leftTest = 0;
        }
        while (!leftTest || !rightTest) {
          if (!rightTest) {
            MPI_Test(&rightReceive,&rightTest,&rightReceiveStatus);
            if (rightTest) loadRightBoundaryFromBuffer(w);
          }
          if (!leftTest) {
            MPI_Test(&leftReceive,&leftTest,&leftReceiveStatus);
            if (leftTest) loadLeftBoundaryFromBuffer(w);
          }
        }
      }
    }

    // compute kinetic term from w, for update routines.
    wfncKinetic(w);

    // first update boundary so that the send can be happening while we update the interior
    updateBoundaries(EPS);

    if (KINTERM==0){    
          	
      // wait and make sure send buffers are ready
      if (i>1) {
        if (nodeID < numNodes) MPI_Wait(&rightSend,&rightSendStatus);
        if (nodeID > 1 ) MPI_Wait(&leftSend,&leftSendStatus);
      }
          	
      // send boundary sync 
      leftTest=1; 
      rightTest=1;
      if (i==1) {
        if (nodeID < numNodes) sendRightBoundary(W);
        if (nodeID > 1) sendLeftBoundary(W);
      }
      else if (i!=nsteps && i>1) {
        if (nodeID < numNodes) rightTest=0;
        if (nodeID > 1) leftTest=0;
        while (!leftTest || !rightTest) {
          if (!rightTest) {
            MPI_Test(&rightSend,&rightTest,&rightSendStatus);
            if (rightTest) sendRightBoundary(W);
          }
          if (!leftTest) {
            MPI_Test(&leftSend,&leftTest,&leftSendStatus);
            if (leftTest) sendLeftBoundary(W);
          }
        }
      }
    }
		
    // evolve interior of the grid forward in time by EPS storing updated fields in capital vars
    updateInterior(EPS);
		
    // copy fields from capital vars (updated) down to lowercase vars (current)
    copyDown();
		
  }
	
}

/*---------------------------------------------------------------------------*/
/* Compute Wavefunction Overlap                                              */
/*---------------------------------------------------------------------------*/

dcomp computeOverlap(dcomp*** wfnc1, dcomp*** wfnc2){
  dcomp overlap=0,overlapCollect=0;
  for (int sx=1;sx<=NUMX;sx++)
    for (int sy=1;sy<=NUM;sy++)
      for (int sz=1; sz<=NUM;sz++)
	overlap += conj(wfnc1[sx][sy][sz])*wfnc2[sx][sy][sz];

  MPI_Reduce(&overlap,&overlapCollect,1,MPI_DOUBLE_COMPLEX,MPI_SUM,0,workers_comm);
  MPI_Bcast(&overlapCollect,1,MPI_DOUBLE_COMPLEX,0,workers_comm);
  return overlapCollect;
}

/*---------------------------------------------------------------------------*/
/* Find excited states                                                       */
/*---------------------------------------------------------------------------*/

void findExcitedStates(const double time, int step) {

  char label[64]; 
  dcomp ener,dEtau;

  /*---------------------------------------------------------------------------*/
  /* Find first excited state                                                  */
  /*---------------------------------------------------------------------------*/
	
  int snap = snapcnt;
	
  // compute overlap
  dcomp overlap = computeOverlap(w,wstore[snap]);

  // subtract overlap
  for (int sx=0;sx<NUMX+2;sx++) 
    for (int sy=0;sy<NUM+2;sy++)
      for (int sz=0; sz<NUM+2;sz++) 
	w1[sx][sy][sz] = wstore[snap][sx][sy][sz] - overlap*w[sx][sy][sz];
	
  // compute observables
  computeObservables(w1);
	
  // normalize
  MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE, 0, workers_comm);
  normalizeWavefunction(w1);

  if (DUMPSNAPS>0 && step>0) {	
    // output snapshot of excited states
    sprintf(label,"%d_1_%d",step,nodeID); 
    outputSnapshot(w1,label);
  }
	
  if (nodeID==1) {

    outputSummaryData("First Excited State",time, step);

    // consistency check and warning if fail
    ener = energyCollect/normalizationCollect;
    dEtau = (ener-EGrnd)*timef;     
    if (real(dEtau)<4.) {
      print_line();
      cout << "==> WARNING: states nearly degenerate, tau_f too small!" << endl;
      print_line();
    }
		
  }
  EGrnd = ener; // redfine EGrnd to hold energy of first excited state for comparison below

  /*---------------------------------------------------------------------------*/
  /* Find second excited state                                                 */
  /*---------------------------------------------------------------------------*/
	
  int save = snap;
  snap = (snapcnt+1)%2;
	
  // compute overlap
  overlap = computeOverlap(wstore[snap],w);
  dcomp overlap2 = computeOverlap(wstore[snap],w1);

  // subtract overlap
  for (int sx=0;sx<NUMX+2;sx++) 
    for (int sy=0;sy<NUM+2;sy++)
      for (int sz=0; sz<NUM+2;sz++) 
	w2[sx][sy][sz] = wstore[snap][sx][sy][sz] - overlap*w[sx][sy][sz] - overlap2*w1[sx][sy][sz];
	
  // compute observables
  computeObservables(w2);
	
  // normalize
  MPI_Bcast(&normalizationCollect, 1, MPI_DOUBLE, 0, workers_comm);
  normalizeWavefunction(w2);

  if (DUMPSNAPS>0 && step>0) {	
    // output snapshot of excited states
    sprintf(label,"%d_2_%d",step,nodeID); 
    outputSnapshot(w2,label);
  }
	
  if (nodeID==1) {

    outputSummaryData("Second Excited State", time, step);

    // consistency check and warning if fail
    ener = energyCollect/normalizationCollect;
    dEtau = (ener-EGrnd)*timef;     
    if (real(dEtau)<4.) {
      print_line();
      cout << "==> WARNING: states nearly degenerate, tau_f too small!" << endl;
      print_line();
    }
		
  }

  /*---------------------------------------------------------------------------*/
  /* Save extracted wavefunctions                                              */
  /*---------------------------------------------------------------------------*/
	
  if (SAVEWAVEFNCS && step<0) {
    // save 3d wavefunction for extracted states
    sprintf(label,"0_%d",nodeID); 
    outputWavefunction(w,label);
    // For now only output ground state, NFS causes delays for writes and this eats time
		
    sprintf(label,"1_%d",nodeID); 
    outputWavefunction(w1,label);
    sprintf(label,"2_%d",nodeID); 
    outputWavefunction(w2,label);
  }
	
  return;
}

/*---------------------------------------------------------------------------*/
/* Boundary sync routines                                                    */
/*---------------------------------------------------------------------------*/

void syncBoundaries(dcomp ***wfnc) {
	
  // initiate sends and receives
  if (nodeID < numNodes) { 
    sendRightBoundary(wfnc);
    receiveRightBoundary(); 
  }
  if (nodeID > 1) {
    sendLeftBoundary(wfnc);
    receiveLeftBoundary(); 
  }
	
  // now wait for communications to complete and sync wfnc when they do
  if (nodeID < numNodes) { 
    MPI_Wait(&rightReceive,&rightReceiveStatus);
    loadRightBoundaryFromBuffer(wfnc);
    MPI_Wait(&rightSend,&rightSendStatus);
  }
  if (nodeID > 1) {
    MPI_Wait(&leftReceive,&leftReceiveStatus);
    loadLeftBoundaryFromBuffer(wfnc);
    MPI_Wait(&leftSend,&leftSendStatus);
  }
}

void sendRightBoundary(dcomp*** wfnc) {
  char message[64]; 
  if (debug == DEBUG_FULL) {
    sprintf(message,"%d -> %d",nodeID,nodeID+1); 
    debug_out << "==> Sending : " << message << endl;
    MPI_Isend(message, strlen(message), MPI_CHAR, nodeID, SYNC_RIGHT_MESSAGE, MPI_COMM_WORLD, &rightMessageSend); 
  }
  for (int sy=0;sy<NUM+2;sy++)
    for (int sz=0;sz<NUM+2;sz++) {
      rightSendBuffer[sy*(NUM+2)+sz%(NUM+2)] = real(wfnc[NUMX][sy][sz]);
      rightSendBuffer[sy*(NUM+2)+sz%(NUM+2) + (NUM+2)*(NUM+2)] = imag(wfnc[NUMX][sy][sz]);
    }
  MPI_Isend(rightSendBuffer, 2*(NUM+2)*(NUM+2), MPI_DOUBLE, nodeID, SYNC_RIGHT, MPI_COMM_WORLD, &rightSend); 
}

void sendLeftBoundary(dcomp*** wfnc) {
  char message[64]; 
  if (debug == DEBUG_FULL) {
    sprintf(message,"%d -> %d",nodeID,nodeID-1); 
    debug_out << "==> Sending : " << message << endl;
    MPI_Isend(message, strlen(message), MPI_CHAR, nodeID-2, SYNC_LEFT_MESSAGE, MPI_COMM_WORLD, &leftMessageSend); 
  }
  for (int sy=0;sy<NUM+2;sy++)
    for (int sz=0;sz<NUM+2;sz++) { 
      leftSendBuffer[sy*(NUM+2)+sz%(NUM+2)] = real(wfnc[1][sy][sz]);
      leftSendBuffer[sy*(NUM+2)+sz%(NUM+2) + (NUM+2)*(NUM+2)] = imag(wfnc[1][sy][sz]);
    }
  MPI_Isend(leftSendBuffer, 2*(NUM+2)*(NUM+2), MPI_DOUBLE, nodeID-2, SYNC_LEFT, MPI_COMM_WORLD, &leftSend);
}

void receiveRightBoundary() {
  char message[64]; 
  if (debug == DEBUG_FULL) {
    MPI_Irecv(message, 255, MPI_CHAR, nodeID, SYNC_LEFT_MESSAGE, MPI_COMM_WORLD, &rightMessageReceive); 
    debug_out << "==> Received : " << message << endl;
  }
  MPI_Irecv(rightReceiveBuffer, 2*(NUM+2)*(NUM+2), MPI_DOUBLE, nodeID, SYNC_LEFT, MPI_COMM_WORLD, &rightReceive); 
}

inline void loadRightBoundaryFromBuffer(dcomp ***wfnc) {
  // update w array right boundary
  for (int sy=0;sy<NUM+2;sy++)
    for (int sz=0;sz<NUM+2;sz++)
      wfnc[NUMX+1][sy][sz] = dcomp(rightReceiveBuffer[sy*(NUM+2)+sz%(NUM+2)],rightReceiveBuffer[sy*(NUM+2)+sz%(NUM+2)+(NUM+2)*(NUM+2)]);
}

void receiveLeftBoundary() {
  char message[64];
  if (debug == DEBUG_FULL) {
    MPI_Irecv(message, 255, MPI_CHAR, nodeID-2, SYNC_RIGHT_MESSAGE, MPI_COMM_WORLD, &leftMessageReceive); 
    debug_out << "==> Received : "<< message << endl;
  }
  MPI_Irecv(leftReceiveBuffer, 2*(NUM+2)*(NUM+2), MPI_DOUBLE, nodeID-2, SYNC_RIGHT, MPI_COMM_WORLD, &leftReceive); 
}

inline void loadLeftBoundaryFromBuffer(dcomp ***wfnc) {
  // update w array left boundary
  for (int sy=0;sy<NUM+2;sy++)
    for (int sz=0;sz<NUM+2;sz++)
      wfnc[0][sy][sz] = dcomp(leftReceiveBuffer[sy*(NUM+2)+sz%(NUM+2)],leftReceiveBuffer[sy*(NUM+2)+sz%(NUM+2)+(NUM+2)*(NUM+2)]);
}

void wfncKinetic(dcomp ***wfnc)
{
  bool N_EVEN;

  switch(KINTERM){
    
  case 0: // NON RELATIVISTIC, FINITE DIREFENCES
    for (int sx=1;sx<=NUMX;sx++)
      for (int sy=1;sy<=NUM;sy++)
  	for (int sz=1;sz<=NUM;sz++)
  	  t_kin[sx][sy][sz] = - ( wfnc[sx+1][sy][sz] + wfnc[sx-1][sy][sz] +
				  wfnc[sx][sy+1][sz] + wfnc[sx][sy-1][sz] +
				  wfnc[sx][sy][sz+1] + wfnc[sx][sy][sz-1] -
				  ((dcomp) 6.)*wfnc[sx][sy][sz])/(((dcomp) 2.)*A*A*MASS);
    
    break;

  case 1: // FFTW: NON RELATIVISTIC
  case 2: // FFTW: NON RELATIVISTIC, WITH SIN(k)
  case 3: // FFTW: RELATIVISTIC

    N_EVEN = NUM%2==0;

    if (fftw_par.local_0_start!=NUMX*(nodeID-1)){
      cout << "ERROR: BAD ALIGNMENT!!!" << endl;
      exit(1);
    }


    MPI_Barrier(workers_comm);   
     
    for (int i=0; i<fftw_par.local_n0; ++i)
      for (int j=0; j<NUM; ++j)
	for (int k=0; k<NUM; ++k)
	  dcomp_to_fftw(wfnc[i+1][j+1][k+1], &fftw_par.in[i*NUM*NUM+j*NUM+k]);
    
    MPI_Barrier(workers_comm);   
    fftw_execute(fftw_par.plan_fftw);
    MPI_Barrier(workers_comm);   

    for (int i=0; i<fftw_par.local_n0; ++i)
      for (int j=0; j<NUM; ++j) 
	for (int k=0; k<NUM; ++k) {
	  dcomp v_ft;
	  double tmp;
	  int ki, kj, kk;
	  
	  ki = fftw_par.local_0_start+i;
	  kj = j;
	  kk = k;
	 
	  if (2*ki > NUM) ki-= NUM;
	  if (2*kj > NUM) kj-= NUM;
	  if (2*kk > NUM) kk-= NUM;
	  
	  // if (N_EVEN){
	  //   if (2*ki==NUM) ki=0;
	  //   if (2*kj==NUM) kj=0;
	  //   if (2*kk==NUM) kk=0;
	  // }
	  
	  fftw_to_dcomp(fftw_par.out[i*NUM*NUM+j*NUM+k], &v_ft);

	  switch(KINTERM){
	    
	  case 1:
	    tmp = pow2(2.*M_PI*((double)ki)/((double)NUM))
	      + pow2(2.*M_PI*((double)kj)/((double)NUM))
	      + pow2(2.*M_PI*((double)kk)/((double)NUM));

	    tmp /= (A*A)*((double) (NUM*NUM*NUM));
	    tmp /= (double)2.*MASS;
	    break;

	  case 2:
	    tmp = 4.*pow2(sin(M_PI*((double)ki)/((double)NUM)))
	      + 4.*pow2(sin(M_PI*((double)kj)/((double)NUM)))
	      + 4.*pow2(sin(M_PI*((double)kk)/((double)NUM)));

	    tmp /= (A*A)*((double) (NUM*NUM*NUM));
	    tmp /= (double)2.*MASS;
	    break;

	  case 3:
	    tmp = 4.*pow2(sin(M_PI*((double)ki)/((double)NUM)))
	      + 4.*pow2(sin(M_PI*((double)kj)/((double)NUM)))
	      + 4.*pow2(sin(M_PI*((double)kk)/((double)NUM)));

	    tmp = sqrt(tmp+MASS*MASS*A*A)/(A*((double)(NUM*NUM*NUM)));
	  
	    break;

	  default:
	    cerr << " >> FATAL ERROR INSIDE void wfncKinetic(dcomp ***wfnc) " << endl;
	    exit(1);
	  }

	  v_ft = v_ft*dcomp(tmp,0.);
	  
	  dcomp_to_fftw(v_ft, &(fftw_par.in[i*NUM*NUM+j*NUM+k]) );
	}

    MPI_Barrier(workers_comm);
    fftw_execute(fftw_par.plan_ifftw);
    MPI_Barrier(workers_comm);   

    for (int i=0; i<fftw_par.local_n0; ++i)
      for (int j=0; j<NUM; ++j)
	for (int k=0; k<NUM; ++k){
	  fftw_to_dcomp(fftw_par.out[i*NUM*NUM+j*NUM+k], &(t_kin[i+1][j+1][k+1]) );
        }
  
    break;

  default:
    cerr << " >> FATAL ERROR!!! UNVALID KINTERM = " << KINTERM << " PARAMETER." << endl;
    exit(1);
  }

  MPI_Barrier(workers_comm);   
}

