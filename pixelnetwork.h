#ifndef PIXELNETWORK_H
#define PIXELNETWORK_H

#include <fstream>
#include <vector>
#include <string>


using namespace std;

class trajectory
{// Basic class for a simple xyz trajectory from LAMMPS

 friend class nodenetwork; // nodenetwork can access private parts

 public:
  
  trajectory(string);
  ~trajectory();
  void setbox(double,double,double,double,double,double);
  void setperiodic(bool,bool,bool);
  void set_pixl(double);
  void set_dia(double);
  void set_threshdens(double);
  void set_xyzgridname(string);
  void do_grid(int,int); // the gridding mechanism    

 private:
  
  ifstream xyzfile;
  ofstream xyzgridfile;

  double dia; // bead diameter from sims

  double xlo,xhi,ylo,yhi,zlo,zhi;
  double vol,Dx,Dy,Dz;
  double threshdens;
  double pixl=0,pixvol; // pixel dimensions
  int Natoms;
  bool xperiodic,yperiodic,zperiodic;
  vector<int> types;
  
  // Grid stuff
  
  int *grid;
  double thresh_dens;
  int Nx, Ny, Nz;
  int gridsize;
  int index1D(int,int,int);
  int indexX(int);
  int indexY(int);
  int indexZ(int);
  int binindexX(double );
  int binindexY(double );
  int binindexZ(double );
  void output_gridxyz(int);
  void output_grid(int);
  void zero_grid();
  void set_grid();

};



#endif
