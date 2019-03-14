// -*- C++ -*-
//
// MamboPhasespace.cc is a part of Herwig - A multi-purpose Monte Carlo event generator
// Copyright (C) 2002-2017 The Herwig Collaboration
//
// Herwig is licenced under version 3 of the GPL, see COPYING for details.
// Please respect the MCnet academic guidelines, see GUIDELINES for details.
//
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MamboPhasespace class.
//

#include "MamboPhasespace.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/Interface/Switch.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "Herwig/Utilities/GSLBisection.h"
#include "Herwig/Utilities/Kinematics.h"
#include "ThePEG/Cuts/Cuts.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MamboPhasespace::MamboPhasespace() :   _maxweight(10.), _a0(10,0.), _a1(10,0.) {}

MamboPhasespace::~MamboPhasespace() {}

IBPtr MamboPhasespace::clone() const {
  return new_ptr(*this);
}

IBPtr MamboPhasespace::fullclone() const {
  return new_ptr(*this);
}

double MamboPhasespace::generateKinematics(vector<Lorentz5Momentum>& P,
						    Energy Ecm,
						    const double* r) const {
  
  double weight = calculateMomentum(P,Ecm);

  //cout << "MAMBO weight = " << weight << endl;
  
  return weight;

}


double MamboPhasespace::calculateMomentum(vector<Lorentz5Momentum> & partmomenta,
					  Energy comEn) const {
  const int N = partmomenta.size()-2;
  vector<Lorentz5Momentum> mom(N);

  
  Energy rmtot(ZERO);
  Energy2 rm2tot(ZERO);
  for(int i = 0;i < N;++i) {
    rmtot += mom[i].mass();
    rm2tot += mom[i].mass2();
  }
  Energy wb = sqrt( (N*sqr(comEn) - sqr(rmtot))/N/N/(N - 1.) ) - rmtot/N;
  Energy wmax = (2.0/3.0)*wb;
  const Energy tol(1e-16*MeV);
  long double r(0.), sf1(0.);
  Energy2 sm2f2(ZERO);
  Energy sf(ZERO), sff1(ZERO), w(ZERO), 
    wold(wmax), err(ZERO);
  unsigned int iter(0), maxiter(50);
  do {
    sf = ZERO; sf1 = 0.; sff1 = ZERO; sm2f2 = ZERO;        
    for(int i = 0;i < N;++i) {
      r = abs(mom[i].mass()/wold);
      Energy f(ZERO);
      long double f1(0.);
      if (r == 0.0) {
	f=2.*wold;
	f1=2.;
      }
      else {
	long double fk0(0.), fkp(0.);
	BesselFns(r, fk0, fkp);
	f = wold*(2.0 + r*fk0);
	f1 = 2.- r*r*fkp;
      }
      sf += f; 
      sf1 += f1; 
      sff1 += f*f1; 
      sm2f2 -= f*f;
    }
    Energy u1 = 2.*(sf*sf1 - sff1);
    Energy2 u0 = sf*sf + sm2f2 + rm2tot - sqr(comEn);
    w = wold - u0/u1;
    err = abs(w - wold);
    wold  = w;
    ++iter;
  }
  while(err > tol && iter < maxiter);
  long double xu,xv;
  vector<long double> alpha(N),um(N),vm(N);
  for(int i = 0;i < N;++i) { 
    alpha[i] = 2.*(mom[i].mass()/w);
    xu = (1.-alpha[i]+sqrt(1.+alpha[i]*alpha[i]))/2.;
    xv = (3.-alpha[i]+sqrt(9.+ 4.*alpha[i]+alpha[i]*alpha[i]))/2.;
    um[i] = exp(-xu/2.)*pow((xu*(xu+alpha[i])),0.25l);
    vm[i] = xv*exp(-xv/2.)*pow((xv*(xv+alpha[i])),0.25l);
  }
  
  //start k-momenta generation
  long double u(0.),v(0.),x(0.);
  vector<Lorentz5Momentum> qk(N);
  Lorentz5Momentum qktot;
  do {
    qktot=LorentzMomentum();
    for(int i=0;i<N;++i) {
      long double usq(0.),bound(0.);
      do {
	u  =  UseRandom::rnd()*um[i];
	v  =  UseRandom::rnd()*vm[i];
	x  =  v/u;
	usq  =  u*u;
	bound  =  exp(-x)*sqrt(x*(x+alpha[i]));
      }
      while(usq>bound);
      double ck,phi;
      Kinematics::generateAngles(ck,phi);
      double sk  =  sqrt(1.0-ck*ck);
      Energy qkv  =  w*sqrt(x*(x+alpha[i]));
      qk[i] = Lorentz5Momentum(qkv*sk*sin(phi),qkv*sk*cos(phi),qkv*ck,
			       mom[i].mass()+w*x);
      qktot += qk[i];
    }
    qktot.rescaleMass();
    x = sqrt(comEn*comEn/qktot.mass2());
  }
  while(x>1.0);

  //Perform lorentz boost from k to q
  vector<Lorentz5Momentum> q(N);
  Energy q0=ZERO, q1=ZERO, q2=ZERO, q3=ZERO;
  long double t=0.;
  vector<Energy2> qsq(N);
  for(int i = 0;i<N;++i){
    q3 = (qk[i]*qktot)/qktot.mass(); 
    t = (q3+qk[i].e())/(qktot.e()+qktot.mass());
    q2 = qk[i].z()-qktot.z()*t; 
    q1 = qk[i].y()-qktot.y()*t;
    q0 = qk[i].x()-qktot.x()*t; 
    q[i] = Lorentz5Momentum(x*q0,x*q1,x*q2,x*q3);
    qsq[i] = sqr(q[i].e())-x*x*mom[i].mass2();
  }
  
  long double xiold(1.),xi(0.); 
  vector<Energy> en(N);
  iter = 0;
  do {
    Energy f = -comEn;
    Energy f1 = ZERO;
    for(int i = 0; i < N; ++i)	    {
      en[i] = sqrt((xiold*xiold*qsq[i]) + mom[i].mass2());
      f += en[i];
      f1 += qsq[i]/en[i];
    }
    xi = xiold - f/(xiold*f1);
    err = abs(xi-xiold)*MeV;
    xiold = xi;
    ++iter;
  }
  while(err > tol && iter < maxiter);
  //Now have desired momenta
  double sumx(0.), sumy(0.), sumz(0.), sume(0.);
  for(int i = 0;i < N;++i) { 
    mom[i] = Lorentz5Momentum(xi*q[i].x(),xi*q[i].y(),xi*q[i].z(),en[i]);
    sumx += mom[i].x()/GeV;
    sumy += mom[i].y()/GeV;
    sumz += mom[i].z()/GeV;
    sume += mom[i].e()/GeV;
  }
  //cout << "sum(x,y,z,e) = " << sumx << ", " << sumy << ", "  << sumz << ", " << sume << endl;
  //Calculate weight of distribution
  double s1(1.);
  Energy s2(ZERO),s3(ZERO);
  double wxi(0.);
  for(int i=0;i<N;++i) {
    s1 *= q[i].e()/mom[i].e();
    s2 += mom[i].mass2()/q[i].e();
    s3 += mom[i].mass2()/mom[i].e();
  }
  wxi = pow(xi,(3*N-3))*s1*(comEn-x*x*s2)/(comEn-s3);
  
  /*sumx = 0.; 
  sumy = 0.; sumz = 0.; sume = 0.;
  Lorentz5Momentum Q = lastPartons().first->momentum() + lastPartons().second->momentum();
  Boost toLab = (lastPartons().first->momentum() + 
		 lastPartons().second->momentum()).boostVector();
  for(int i = 0;i < N;++i) { 
    mom[i].boost(toLab);
    sumx += mom[i].x()/GeV;
    sumy += mom[i].y()/GeV;
    sumz += mom[i].z()/GeV;
    sume += mom[i].e()/GeV;
    //mom[i].rescaleEnergy();
  }
  cout << "toLab(x,y,z)  = " << toLab.x()<< ", " << toLab.y() << ", "  << toLab.z() <<  endl;
  cout << "Q(x,y,z,e)  = " << Q.x()/GeV << ", " << Q.y()/GeV << ", "  << Q.z()/GeV << ", " << Q.e()/GeV <<  endl;

  cout << "sum(x,y,z,e) (after boost) = " << sumx << ", " << sumy << ", "  << sumz << ", " << sume << endl;
  */
  for ( int i = 0; i < N; i++ ) {
    partmomenta[i+2] = mom[i];
  }
  return wxi;
}

double MamboPhasespace::generateTwoToNKinematics(const double* r,
							  vector<Lorentz5Momentum>& momenta) {

  double weight = 1.;

  weight *= generateKinematics(momenta,sqrt(lastXCombPtr()->lastSHat()),r);

  return weight;

}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MamboPhasespace::persistentOutput(PersistentOStream& os) const {
  os << _maxweight << _a0 << _a1 ;
}

void MamboPhasespace::persistentInput(PersistentIStream& is, int) {
  is >> _maxweight >> _a0 >> _a1 ;
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MamboPhasespace,MatchboxPhasespace>
  describeHerwigMamboPhasespace("Herwig::MamboPhasespace", "Herwig.so");

void MamboPhasespace::Init() {

  static ClassDocumentation<MamboPhasespace> documentation
    ("MamboPhasespace implements the Mambo alogorithm for phase-space generation.");

  static Parameter<MamboPhasespace,double> interfaceMaximumWeight
    ("MaxWeight",
     "Maximum phase-space weight",
     &MamboPhasespace::_maxweight, 10.0, 1.0, 50.,
     false, false, true);

  
}

