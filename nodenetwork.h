#ifndef NODE_H
#define NODE_H

#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include "pixelnetwork.h"

using namespace std;


// Number of beads making a rod (BELOW)
// IMPORTANT !!!!!!!!!
#define NB 36

class rod;
class cluster;
class strand;
class amorphon;

class nodenetwork: public trajectory
{

 public:
  
 nodenetwork(string fname): trajectory(fname) { }
 
  void do_network(int,int); 
  void do_classify(int);
  void outputgnu(int,int);
  void set_touchperc(double);
  void set_anglethreshstrand(double);
  void set_anglethreshclust(double);
  void set_rodsthreshstrand(int);
  void add_particletype(int);

private:
  
  const string dir="../nodenetworkout/";


  ofstream clusters_out;
  ofstream strands_out;
  ofstream amorphons_out;

  vector<cluster> clusters;
  vector<strand> strands;
  vector<amorphon> amorphons;

  vector<rod> rods;

  vector<int> ptypes; // particle types from LAMMPS  

  void reset();
  void sort_strands();
 
  
  
  int rodsthreshstrand;// for use in classifying strands
  double anglethreshstrand; // for use in classifying strands
  double anglethreshclust; // for use in classifying clusters
  double touchperc;// number to times the diameter inorder to see if two beads are touching
  

};


class cluster
{

  friend class strand;
  friend class amorphon;
 public:
  cluster(){ }
  cluster(rod,rod);

  virtual void plotcluster(int,int,int);

  int rsize();
  int csize();
  rod getrod(int );

  int getrodid(int );
  void addrod(rod );
  void addrodid(int);
  void setall_ids(int);

  vector<int> getrod_ids();
  vector<rod> getrods();

  bool operator==(cluster& a);
  cluster operator=(cluster a);
  cluster operator+(cluster a);

  void blacklist();
  bool status();

 private:

  bool blacklisted=false;  

  vector<rod> rods;
  

  vector<int> rod_ids;
};


class bead
{

 public:

  void setpos(double,double,double);
  double getx();
  double gety();
  double getz();

 private:

  double x,y,z;

};

class rod
{

 public:
  
  void setbead(int,double,double,double);
  
  double dot(rod); // scalar product
  double ang(rod,int);// angle between rods

  bool connected(rod,double);
  bool thesame(rod);

  void setCOM();
  void setBox(
	      double,double,double,
	      bool,bool,bool
);
  

  void addneigh(int);
  int getneigh();
  
  void belong();
  bool getbelong();
  int mycluster_id();
  void setcluster_id(int);


  void belongstrand();
  bool getbelongstrand();
  int mystrand_id();
  void setstrand_id(int);

  void belongamorph();
  bool getbelongamorph();
  int myamorph_id();
  void setamorph_id(int);



  double PBCx(double);
  double PBCy(double);
  double PBCz(double);


  double getlen();

  double getcomx();
  double getcomy();
  double getcomz();

  double vecx();
  double vecy();
  double vecz();

  double topx(); //top of rod is bead # NB
  double topy();
  double topz();
  
  double botx(); // bottom of rod is bead # 1
  double boty();
  double botz();

  double getbeadx(int);
  double getbeady(int);
  double getbeadz(int);

  static int nbeads(){return NB;};
 private:

  vector<int> neighs; // list of rods that are connected to me

  //int NB; // # beads in a rod is a constant

  static double mag(rod);

  bool belongs=false;
  bool belongsstrand=false;
  bool belongsamorphon=false;

  int cluster_id;
  int strand_id;
  int amorphon_id;

  double xcom,ycom,zcom; // center of mass
  double Dx,Dy,Dz;
  bool xper,yper,zper;
  bead beads[NB];

};


class strand: public cluster
{
 public:
  
 strand(): cluster() { }
 strand(rod a,rod b): cluster(a,b) { }
  
  void plotcluster(int,int,int);

 private:
};

class amorphon: public cluster
{

 public:

 amorphon(): cluster() { }
 amorphon(rod a,rod b): cluster(a,b) { }

  void plotcluster(int,int,int);

 private:


};

#endif
