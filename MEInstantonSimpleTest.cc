// -*- C++ -*-
//
// This is the implementation of the non-inlined, non-templated member
// functions of the MEInstantonSimpleTest class.
//

#include "MEInstantonSimpleTest.h"
#include "ThePEG/Interface/ClassDocumentation.h"
#include "ThePEG/EventRecord/Particle.h"
#include "ThePEG/Repository/UseRandom.h"
#include "ThePEG/Repository/EventGenerator.h"
#include "ThePEG/Utilities/DescribeClass.h"
#include "ThePEG/Interface/Parameter.h"
#include "ThePEG/Interface/Switch.h"

#include "ThePEG/Persistency/PersistentOStream.h"
#include "ThePEG/Persistency/PersistentIStream.h"

using namespace Herwig;

MEInstantonSimpleTest::MEInstantonSimpleTest() : theNQuarkPair(4), theColourConnections(0) {}

MEInstantonSimpleTest::~MEInstantonSimpleTest() {}

IBPtr MEInstantonSimpleTest::clone() const {
  return new_ptr(*this);
}

IBPtr MEInstantonSimpleTest::fullclone() const {
  return new_ptr(*this);
}

double MEInstantonSimpleTest::me2() const {
  //cout << "meMomenta().size() = " << meMomenta().size() << endl;
  // the square of the matrix element 
  return exp(7*(meMomenta().size()-2-2*nQuarkPair())*log(10));
}

void MEInstantonSimpleTest::doinit() {
  //the number of maximum gluons is given by the nAdditional() number of extra partons 
  ngluonmax(this->nAdditional());
}

void MEInstantonSimpleTest::doinitrun() {
  //the number of maximum gluons is given by the nAdditional() number of extra partons 
  ngluonmax(this->nAdditional());
}

multimap<tcPDPair,tcPDVector> MEInstantonSimpleTest::processes() const {
  //the processmap to return 
  multimap<tcPDPair,tcPDVector> processmap;
  
  //define the particles to be used in this process
  tcPDPtr g = getParticleData(ParticleID::g);
  vector<tcPDPtr> q, qb;
  for(int i = 1; i <= 5; ++i) {
    tcPDPtr quark = getParticleData(i);
    q.push_back(quark);
    qb.push_back(quark->CC());
  }

  //make the incoming partons
  tcPDPair incoming = make_pair(g, g);

  //the vectors that hold the outgoing partons
  tcPDVector outgoing;

  //push all the quarks and anti-quarks (one of each for the instanton processes). 
  for(int i = 0; i < nQuarkPair(); ++i) {
    outgoing.push_back(q[i]);
    outgoing.push_back(qb[i]);
  }

  // make the process using the incoming and outgoing particles so far (i.e. no additional gluons)
  processmap.insert(make_pair(incoming,outgoing));

  //insert (1 to ngluon_max) gluons into the process
  for(unsigned int jj=0; jj < ngluon_max; jj++) { 
    outgoing.push_back(g);
    processmap.insert(make_pair(incoming,outgoing));
  }

  return processmap;
}

list<BlobMEBase::ColourConnection> MEInstantonSimpleTest::colourConnections() const {
  list<BlobMEBase::ColourConnection> res;
  
  //count the number of gluons in the given event 
  int ngluon = (meMomenta().size()-2-nQuarkPair()*2);


  if(theColourConnections==0) {
    /*  a simple choice:qqbar pairs apart from the last one in the case of odd number of gluons.
        Then the last pair is connected to the first gluon. The rest of the gluons are 
        paired together. 
    */
    //make the colour connections for the incoming gluons
    //initial-state gluons are connected to each other (i.e. colour singlet)
    BlobMEBase::ColourConnection first; BlobMEBase::ColourConnection second;
    first.addColour(0); first.addAntiColour(1);
    second.addColour(1); second.addAntiColour(0);
    res.push_back(first);
    res.push_back(second);
  

    //make the colour connections for some of the quark lines
    vector<BlobMEBase::ColourConnection> quark_lines;
    quark_lines.resize(4);
    for(unsigned int cc = 0; cc < nQuarkPair()-1; cc++) {
      quark_lines[cc].addColour(2*cc+2);
      quark_lines[cc].addAntiColour(2*cc+3);
      res.push_back(quark_lines[cc]);
    }
    //if even number of gluons then connect the remaining quarks as well
    //and connect the gluons in pairs
    if(ngluon%2 == 0) {
      BlobMEBase::ColourConnection quark_lastpair;
      quark_lastpair.addColour(2*(nQuarkPair()-1)+2);
      quark_lastpair.addAntiColour(2*(nQuarkPair()-1)+3);
      res.push_back(quark_lastpair);
      //loop over gluons and pair adjacent ones 
      for(unsigned int gg = 1; gg <= ngluon/2; gg++) {
        BlobMEBase::ColourConnection firstgl;
        BlobMEBase::ColourConnection secondgl;
        firstgl.addColour(2*nQuarkPair()+1+2*gg-1);
        firstgl.addAntiColour(2*nQuarkPair()+2+2*gg-1);
        secondgl.addColour(2*nQuarkPair()+2+2*gg-1);
        secondgl.addAntiColour(2*nQuarkPair()+1+2*gg-1);
        res.push_back(firstgl);
        res.push_back(secondgl);      
      }  
    } else {
      //if odd number of gluons, connect the last quark pair to one of them
      // and the rest in pairs as before. 
      BlobMEBase::ColourConnection quark_lastpair1;
      BlobMEBase::ColourConnection quark_lastpair2;
    
      quark_lastpair1.addColour(2*(nQuarkPair()-1)+2);
      quark_lastpair1.addAntiColour(2*(nQuarkPair()-1)+2+2);
      res.push_back(quark_lastpair1);
      quark_lastpair2.addAntiColour(2*(nQuarkPair()-1)+3);
      quark_lastpair2.addColour(2*(nQuarkPair()-1)+2+2);
      res.push_back(quark_lastpair2);
      for(unsigned int gg = 1; gg <= (ngluon-1)/2; gg++) {
        BlobMEBase::ColourConnection firstg;
        BlobMEBase::ColourConnection secondg;
        firstg.addColour(nQuarkPair()*2+2+2*gg-1);
        firstg.addAntiColour(nQuarkPair()*2+3+2*gg-1);
        secondg.addColour(nQuarkPair()*2+3+2*gg-1);
        secondg.addAntiColour(nQuarkPair()*2+2+2*gg-1);
        res.push_back(firstg);
        res.push_back(secondg);      
      }
    } 
  } else if(theColourConnections==1) { //purely random selection (colour singlet gg only)
    //construct the array of colour and anticolour numbers
    vector<int> colours; vector<int> anticolours; vector<int> colourmap;
    // cout << "number of particles = " << meMomenta().size() << " number of gluons = " << ngluon << endl;
    //loop over the qqbar
    for(unsigned int cc = 0; cc < nQuarkPair(); cc++) {
      colours.push_back(2*cc+2);
      anticolours.push_back(2*cc+3);
    }
    //loop over the gluons
    for(int gg = 0; gg < ngluon; gg++) {
      colours.push_back(2+nQuarkPair()*2+gg);
      anticolours.push_back(2+nQuarkPair()*2+gg);
    }
    //pick a random element of the anticolour array that has not already been used. 
    int col = UseRandom::rnd(0, int(colours.size()));
    for(int pp = 0; pp < colours.size(); pp++) {
      col = UseRandom::rnd(0, int(anticolours.size()));
      if(anticolours.size()>1) { 
        while(colours[pp]==anticolours[col]) {
          col = UseRandom::rnd(0, int(anticolours.size()));
          /*cout << "anticolours.size() = " << anticolours.size() << endl;
          cout << "col chosen = " << col << endl;
          cout << "trying to connect " << colours[pp] << " to -" << anticolours[col] << endl;*/
        }
        colourmap.push_back(anticolours[col]);
        anticolours.erase(anticolours.begin()+col);
        continue; 
      } else if(anticolours.size()==1) { 
        if(colours[pp]==anticolours[col]) {
          // cout << "only one colour/anti-colour left and they belong to the same gluon!" << endl;
          //pick a random element of the colourmap 
          int switchcol = UseRandom::rnd(0, int(colourmap.size())-1);
          //save the previous anticolour that was placed there
          int acolourold = colourmap[switchcol];
          /*cout << "switching a randomly-chosen colour, element " << switchcol << " corresponding to colour " <<  colours[switchcol] << " to connect to anticolour " << anticolours[col] << endl;
            cout << "this used to correspond to anticolour " << acolourold << endl;*/
          //push back the new connection 
          colourmap.push_back(acolourold);
          //and change the old one to the last anticolour
          colourmap[switchcol] = anticolours[col];
          continue;
        } else {
          colourmap.push_back(anticolours[col]);
          continue; 
        }
      }
    }
    for(int ii = 0; ii < colourmap.size(); ii++) {
      BlobMEBase::ColourConnection conline;
      conline.addColour(colours[ii]);
      conline.addAntiColour(colourmap[ii]);
      //  cout << "connecting " << colours[ii] << " to -" << colourmap[ii] << endl;
      res.push_back(conline);
    }
    
    //make the colour connections for the incoming gluons
    BlobMEBase::ColourConnection first; BlobMEBase::ColourConnection second;
    first.addColour(0); first.addAntiColour(1);
    second.addColour(1); second.addAntiColour(0);
    res.push_back(first);
    res.push_back(second);
    // cout << "done connecting randomly" << endl;
  } else if(theColourConnections==2) { //purely random selection (singlet + octet gg)
    //construct the array of colour and anticolour numbers
    vector<int> colours; vector<int> anticolours; vector<int> colourmap;
    //   cout << "number of particles = " << meMomenta().size() << " number of gluons = " << ngluon << endl;
    //loop over the qqbar
    colours.push_back(0);
    anticolours.push_back(1);
    for(unsigned int cc = 0; cc < nQuarkPair(); cc++) {
      colours.push_back(2*cc+2);
      anticolours.push_back(2*cc+3);
    }
    //loop over the gluons
    for(int gg = 0; gg < ngluon; gg++) {
      colours.push_back(2+nQuarkPair()*2+gg);
      anticolours.push_back(2+nQuarkPair()*2+gg);
    }
    //pick a random element of the anticolour array that has not already been used. 
    int col = UseRandom::rnd(0, int(colours.size()));
    for(int pp = 0; pp < colours.size(); pp++) {
      col = UseRandom::rnd(0, int(anticolours.size()));
      if(anticolours.size()>1) { 
        while(colours[pp]==anticolours[col]) {
          col = UseRandom::rnd(0, int(anticolours.size()));
          /*cout << "anticolours.size() = " << anticolours.size() << endl;
          cout << "col chosen = " << col << endl;
          cout << "trying to connect " << colours[pp] << " to -" << anticolours[col] << endl;*/
        }
        colourmap.push_back(anticolours[col]);
        anticolours.erase(anticolours.begin()+col);
        continue; 
      } else if(anticolours.size()==1) { 
        if(colours[pp]==anticolours[col]) {
          // cout << "only one colour/anti-colour left and they belong to the same gluon!" << endl;
          //pick a random element of the colourmap 
          int switchcol = UseRandom::rnd(0, int(colourmap.size())-1);
          //save the previous anticolour that was placed there
          int acolourold = colourmap[switchcol];
          /*cout << "switching a randomly-chosen colour, element " << switchcol << " corresponding to colour " <<  colours[switchcol] << " to connect to anticolour " << anticolours[col] << endl;
            cout << "this used to correspond to anticolour " << acolourold << endl;*/
          //push back the new connection 
          colourmap.push_back(acolourold);
          //and change the old one to the last anticolour
          colourmap[switchcol] = anticolours[col];
          continue;
        } else {
          colourmap.push_back(anticolours[col]);
          continue; 
        }
      }
    }
    for(int ii = 0; ii < colourmap.size(); ii++) {
      BlobMEBase::ColourConnection conline;
      if(colourmap[ii] == 1 && colours[ii]!= 0) {
        conline.addColour(colourmap[ii]);
        conline.addColour(colours[ii]);
	//        cout << "connecting " << colours[ii] << " to " << colourmap[ii] << endl;

      }
      if(colours[ii] == 0 && colourmap[ii]!=1) {
        conline.addAntiColour(colours[ii]);
        conline.addAntiColour(colourmap[ii]);
        //cout << "connecting -" << colours[ii] << " to -" << colourmap[ii] << endl;
                
      }
      if (colourmap[ii] != 1 && colours[ii] != 0){ 
        conline.addColour(colours[ii]);
        conline.addAntiColour(colourmap[ii]);
        //cout << "connecting " << colours[ii] << " to -" << colourmap[ii] << endl;

      }
      if (colourmap[ii] == 1 && colours[ii] == 0){ 
        conline.addAntiColour(colours[ii]);
        conline.addColour(colourmap[ii]);
	//        cout << "connecting -" << colours[ii] << " to " << colourmap[ii] << endl;
        
      }
      res.push_back(conline);
    }
    
    //make the colour connections for the incoming gluons
    BlobMEBase::ColourConnection first; 
    first.addColour(0); first.addAntiColour(1);
    //    cout << "connecting 0 to -1" << endl;

    res.push_back(first);
    // cout << "done connecting randomly" << endl;
  }

  return res;
}

size_t MEInstantonSimpleTest::nOutgoing() const {
  return (2*nQuarkPair());
}

// If needed, insert default implementations of virtual function defined
// in the InterfacedBase class here (using ThePEG-interfaced-impl in Emacs).


void MEInstantonSimpleTest::persistentOutput(PersistentOStream & os) const {
  // *** ATTENTION *** os << ; // Add all member variable which should be written persistently here.
  os << theNQuarkPair << ngluon_max << theColourConnections;

}

void MEInstantonSimpleTest::persistentInput(PersistentIStream & is, int) {
  // *** ATTENTION *** is >> ; // Add all member variable which should be read persistently here.
  is >> theNQuarkPair >> ngluon_max >> theColourConnections;
    
}


// *** Attention *** The following static variable is needed for the type
// description system in ThePEG. Please check that the template arguments
// are correct (the class and its base class), and that the constructor
// arguments are correct (the class name and the name of the dynamically
// loadable library where the class implementation can be found).
DescribeClass<MEInstantonSimpleTest,Herwig::BlobME>
  describeHerwigMEInstantonSimpleTest("Herwig::MEInstantonSimpleTest", "Instantons.so");

void MEInstantonSimpleTest::Init() {
  
  static ClassDocumentation<MEInstantonSimpleTest> documentation
    ("There is no documentation for the MEInstantonSimpleTest class");

  static Parameter<MEInstantonSimpleTest,size_t> interfaceNQuarkPair
    ("NQuarkPair",
     "The number of quark pairs to consider.",
     &MEInstantonSimpleTest::theNQuarkPair, 4, 1, 6,
     false, false, Interface::limited);

  static Switch<MEInstantonSimpleTest,unsigned int> interfaceColourConnections
    ("ColourConnections",
     "How to connect the colour lines",
     &MEInstantonSimpleTest::theColourConnections, 0, false, false);
  static SwitchOption interfaceColourConnectionsSimple
    (interfaceColourConnections,
     "Simple",
     "A very simple arbitrary choice.",
     0);
  static SwitchOption interfaceColourConnectionsRandom
    (interfaceColourConnections,
     "Random",
     "Completely randomized colour connections. Singlet gg only.",
     1);
    static SwitchOption interfaceColourConnectionsRandom2
    (interfaceColourConnections,
     "Random2",
     "Completely randomized colour connections. Singlet and octet gg.",
     2);

}

