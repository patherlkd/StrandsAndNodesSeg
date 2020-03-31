#include "pixelnetwork.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

using namespace std;

trajectory::trajectory(string filename){

  xyzfile.open(filename);

}

trajectory::~trajectory(){

  if(xyzfile.is_open())
    xyzfile.close();
  
  if(xyzgridfile.is_open())
    xyzgridfile.close();

  grid = new int[1];
  delete[] grid;

}


void trajectory::do_grid(int interval,int Ntraj){

  cout << " ================= \n ";
  cout << " do_grid("<<interval<<")     \n";
  cout << " ================= \n ";

  if(pixl==0){
    cout<<"Error: pixel length not set\n";
    return;
  }
  
  int type;
  double x,y,z;
  int traj=0;
  string line;


  int atom=1;

  set_grid();

  while(getline(xyzfile,line))
    {

      stringstream ss(line);
      ss >> type >> x >> y >> z;
      
      if(traj==0){
	cout<<"Natoms ="<<type<<"\n";
	Natoms=type;
	traj++;
	continue;
      }
      else if(Natoms==type){
	atom=1;
	zero_grid();
	traj++;
      }
      
      
      
      if(traj%interval==0){

      	
	if(type!=Natoms && type!=0){
	 
	   grid[index1D(binindexX(x),binindexY(y),binindexZ(z))]++;
	   
	   if(atom>=Natoms)
	     {
	       output_gridxyz(traj);
	       output_grid(traj);
	     } 
	  atom++;
	}

	
      }
      

     
    }
  
}

void trajectory::set_xyzgridname(string filename){

  xyzgridfile.open(filename);

}

void trajectory::output_gridxyz(int traj){

  cout<<"Updating xyz grid file"<<'\n';

  int N=0;
   for(int i=0;i<gridsize;i++){
     if(((double)grid[i]/pixvol)>threshdens)
       N++;
   }


  xyzgridfile <<Natoms<<'\n';
  xyzgridfile <<"Atoms. Timestep: "<<traj<<'\n';
  double pixh=pixl*0.5;
  
  double x,y,z;
  
  for(int i=0;i<gridsize;i++){
    
    x=pixl*(double)indexX(i)+pixh;
    y=pixl*(double)indexY(i)+pixh;
    z=pixl*(double)indexZ(i)+pixh;

    if(((double)grid[i]/pixvol)>threshdens)
      xyzgridfile<<1<<' '<<x+xlo<<' '<<y+ylo<<' '<<z+zlo<<'\n';
    }

  for(int n=0;n<(Natoms-N);n++){ // pad xyz file
   xyzgridfile<<2<<' '<<0.0<<' '<<0.0<<' '<<0.0<<'\n';
  }

  

}

void trajectory::zero_grid(){

   for(int i=0;i<gridsize;i++){
     grid[i]=0;
   }


}

void trajectory::output_grid(int traj){

  ostringstream os;

  os<<"heatgrid_pixl"<<pixl<<"_boxside"<<Dx<<"_traj"<<traj<<".txt";

  string file = os.str();
  
  ofstream heatgridfile(file);
  double pixh=pixl*0.5;
  
  cout << " ================= \n ";
  cout << " outputting grid     \n";
  cout << " ================= \n ";

  double x,y,z;
  
  for(int i=0;i<gridsize;i++){
    
    x=pixl*(double)indexX(i)+pixh;
    y=pixl*(double)indexY(i)+pixh;
    z=pixl*(double)indexZ(i)+pixh;

    
    //      heatgridfile<<x+xlo<<'\t'<<y+ylo<<'\t'<<z+zlo<<'\t'<<((double)grid[i])<<'\n';
    heatgridfile<<x+xlo<<'\t'<<y+ylo<<'\t'<<z+zlo<<'\t'<<(double)grid[i]/pixvol<<'\n';
    
    }

  heatgridfile.close();

}

void trajectory::set_grid(){

  gridsize=Nx*Ny*Nz;
  grid = new int[gridsize];

}

void trajectory::set_threshdens(double dens)
{

  if(dens>1.0 || dens<0.0)
    {
      cout<<"Error: set_threshdens can only accept densities between 0.0 and 1.0\n";
      exit(1);
    }

  threshdens=dens;

}

void trajectory::setbox(double xlo,double xhi,double ylo,double yhi,double zlo, double zhi){

  this->xlo=xlo;
  this->ylo=ylo;
  this->zlo=zlo;
  this->xhi=xhi;
  this->yhi=yhi;
  this->zhi=zhi;

  Dx=xhi-xlo;
  Dy=yhi-ylo;
  Dz=zhi-zlo;

  vol=Dx*Dy*Dz;
  cout<<"Vol: "<<vol<<'\n';
}

void trajectory::setperiodic(bool x,bool y,bool z){

  xperiodic=x;
  yperiodic=y;
  zperiodic=z;

}

void trajectory::set_pixl(double pixl){

  this->pixl=pixl;
  cout<<"Pixel Length: "<<pixl<<'\n';
  
  pixvol=pixl*pixl*pixl;

  Nx=(int)(Dx/pixl);
  Ny=(int)(Dy/pixl);
  Nz=(int)(Dz/pixl);

  cout<<"(Nx,Ny,Nz): ("<<Nx<<","<<Ny<<","<<Nz<<")"<<'\n';
  
}

void trajectory::set_dia(double dia){

  this->dia=dia;

}

int trajectory::index1D(int ix,int iy,int iz){

  return (ix + Nx*(iy + iz*Ny)); 
}

int trajectory::indexX(int i){

  return i%Nx;

}

int trajectory::indexY(int i){

  return (i/Nx)%Ny;

}

int trajectory::indexZ(int i){

  return i/(Nx*Ny);

}

int trajectory::binindexX(double x){

  int mode;

  if(x-xlo > xhi-x)
    mode=1; // closer to xhi
  else
    mode=0;

  x-=xlo;

 
  double prev;
  double next;
  int bin=0;
  
  if(mode=1){

    prev=Dx;
    next=Dx-pixl;

  while(x>=0){

    if(x<prev && x>=next)
      return bin;

    prev-=pixl;
    next-=pixl;
    bin++;
  }
  }
  else{
    prev=0.0;
    next=pixl;

    while(x<=Dx){
      
    if(x>prev && x<=next)
      return bin;

    prev+=pixl;
    next+=pixl;
    bin++;
    }

  }


  cout<<"No bin value found in X exitting\n";

}

int trajectory::binindexY(double x){


  int mode;

  if(x-ylo > yhi-x)
    mode=1; // closer to xhi
  else
    mode=0;

  x-=ylo;

 
  double prev;
  double next;
  int bin=0;
  
  if(mode=1){

    prev=Dy;
    next=Dy-pixl;

  while(x>=0){

    if(x<prev && x>=next)
      return bin;

    prev-=pixl;
    next-=pixl;
    bin++;
  }
  }
  else{
    prev=0.0;
    next=pixl;

    while(x<=Dy){
      
    if(x>prev && x<=next)
      return bin;

    prev+=pixl;
    next+=pixl;
    bin++;
    }

  }

}


int trajectory::binindexZ(double x){


  int mode;

  if(x-zlo > zhi-x)
    mode=1; // closer to xhi
  else
    mode=0;

  x-=zlo;

 
  double prev;
  double next;
  int bin=0;
  
  if(mode=1){

    prev=Dz;
    next=Dz-pixl;

  while(x>=0){

    if(x<prev && x>=next)
      return bin;

    prev-=pixl;
    next-=pixl;
    bin++;
  }
  }
  else{
    prev=0.0;
    next=pixl;

    while(x<=Dz){
      
    if(x>prev && x<=next)
      return bin;

    prev+=pixl;
    next+=pixl;
    bin++;
    }

  }

  cout<<"No bin value found in Z for z:"<<x<<"for zlo: "<<zlo<<"and zhi: "<<zhi<<"exitting\n";

}
