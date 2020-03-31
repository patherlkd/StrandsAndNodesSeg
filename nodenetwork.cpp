#define _USE_MATH_DEFINES

#include "nodenetwork.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>

using namespace std;



void nodenetwork::add_particletype(int type)
{
  ptypes.push_back(type);
}




void nodenetwork::do_network(int interval, int Ntraj){

  clusters_out.open(dir+"clusters.dat");
  strands_out.open(dir+"strands.dat");
  amorphons_out.open(dir+"amorphons.dat");


  string line;
  int traj=0;  
  int numbead=0;
  int numrods=0;

  while(getline(xyzfile,line))
    {
      double x,y,z;
      int a;
      stringstream ss(line);
      bool linetest=false;      
      

      ss >> a >> x >> y >>z;
      
      //      cout<< x<< '\t' << y<< '\t' << z << '\n';
	if(traj==0){
	  Natoms=a;
	  cout<<"a: "<<a<<'\n';
	  reset();
	  traj++;
	  continue;
	}
	
	if(traj>Ntraj)
	  return;


	
	
	if(Natoms==a){

	  numrods=0;
	  numbead=0;
	  if(traj%interval==0 && traj>0){// actually do network stuff here
	    //   outputgnu(traj,1); // send a file to check rod COM


	    // MAIN GUTS OF PROGRAM
	    do_classify(traj);
	    
	    
	    reset();
	  }
	  traj++;
	}
	
	
	int pcnt=0;
	while(pcnt<ptypes.size())
	  {
	    if(a==ptypes[pcnt])
	      {linetest=true; break;}
	    pcnt++;
	  }

	
      if(traj%interval==0 && a!=Natoms && linetest){

	
	if(numbead>rod::nbeads()-1)
	  {
	    numbead=0;
	    //cout<<"numrods: "<<numrods<<'\n';
	    rods[numrods].setCOM();
	    rods[numrods].setBox(Dx,Dy,Dz,xperiodic,yperiodic,zperiodic);
	    numrods++;
	  }
	
	//cout<<"numbeads: "<<numbead<<'\n';

	rods[numrods].setbead(numbead,x,y,z);
	numbead++;
	
      }


    }
    

  clusters_out.close();
  strands_out.close();
  amorphons_out.close();

}

void nodenetwork::reset()
{

  rods.clear();
  clusters.clear();
  amorphons.clear();
  strands.clear();

  rods.resize(Natoms/rod::nbeads());
  clusters.resize(0);
  strands.resize(0);
  amorphons.resize(0);
}

void nodenetwork::outputgnu(int traj,int mode){

  ostringstream os;
  os<<dir<<"rods_com_traj"<<traj<<".txt";
 
  ofstream out(os.str());
  if(mode==0){
    for(int i=0;i<rods.size();i++){
      out<<rods[i].getcomx()<<'\t'<<rods[i].getcomy()<<'\t'<<rods[i].getcomz()<<'\n';
    }
    
  }
  else if(mode==1){
    
    
    for(int i=0;i<rods.size();i++)
      for(int j=0;j<rod::nbeads();j++){
	
	out<<rods[i].getbeadx(j)<<'\t'<<rods[i].getbeady(j)<<'\t'<<rods[i].getbeadz(j)<<'\n';
	
      }
    
  }
  else{
    out<<"WRONG MODE in outputgnu(traj,MODE).\n";
  }
  
  out.close();
  
}

void nodenetwork::set_rodsthreshstrand(int numrods){

  rodsthreshstrand=numrods;
}

void nodenetwork::set_touchperc(double touchperc){

  this->touchperc=touchperc*dia;
}

void nodenetwork::set_anglethreshclust(double anglethresh)
{
  this->anglethreshclust=anglethresh;
}


void nodenetwork::set_anglethreshstrand(double anglethresh)
{
  this->anglethreshstrand=anglethresh;
}


void nodenetwork::sort_strands()
{

  for(int i=0;i<strands.size();i++)
    {
      if(strands[i].csize()<rodsthreshstrand)
	strands[i].blacklist();
     
    }

}

cluster::cluster(rod a,rod b){
  rods.push_back(a);
  rods.push_back(b);
}

void strand::plotcluster(int stra,int mode, int extra)
{
 ostringstream os;
  if(mode==0)
    {// plot beads in cluster
      os<<"../nodenetworkout/strand_"<<stra<<"traj_"<<extra<<"beads_xyz.txt";

  ofstream out(os.str());

  for(int i=0;i<rods.size();i++)
    for(int j=0;j<rods[i].nbeads();j++){
      out << rods[i].getbeadx(j) << '\t' << rods[i].getbeady(j) << '\t' <<rods[i].getbeadz(j)<<'\n';
    }

  out.close();
    }
  else if(mode==1)
    {// plot COM of rods in cluster
      os<<"../nodenetworkout/strand_"<<stra<<"traj_"<<extra<<"COM_xyz.txt";;

      ofstream out(os.str());

      for(int i=0;i<rods.size();i++)
	{	
	  out << rods[i].getcomx() << '\t' << rods[i].getcomy() << '\t' <<rods[i].getcomz()<<'\n';
	}

      out.close();    
    }
  else if(mode==2)
    {// plot system as vectors

      os<<"../nodenetworkout/strand_"<<stra<<"traj_"<<extra<<"VEC_xyz.txt";;

      ofstream out(os.str());

      for(int i=0;i<rods.size();i++)
	{	
	  out << rods[i].botx() << '\t';
	  out << rods[i].boty() << '\t';
	  out << rods[i].botz() << '\t';
	  
	  out << rods[i].PBCx(rods[i].topx()-rods[i].botx()) << '\t';
	  out << rods[i].PBCy(rods[i].topy()-rods[i].boty()) << '\t';
	  out << rods[i].PBCz(rods[i].topz()-rods[i].botz()) << '\n';
	  
	}
      
      out.close();


    }

}

void amorphon::plotcluster(int stra,int mode, int extra)
{
 ostringstream os;
  if(mode==0)
    {// plot beads in cluster
      os<<"../nodenetworkout/amorphon_"<<stra<<"traj_"<<extra<<"beads_xyz.txt";

  ofstream out(os.str());

  for(int i=0;i<rods.size();i++)
    for(int j=0;j<rods[i].nbeads();j++){
      out << rods[i].getbeadx(j) << '\t' << rods[i].getbeady(j) << '\t' <<rods[i].getbeadz(j)<<'\n';
    }

  out.close();
    }
  else if(mode==1)
    {// plot COM of rods in cluster
      os<<"../nodenetworkout/amorphon_"<<stra<<"traj_"<<extra<<"COM_xyz.txt";;

      ofstream out(os.str());

      for(int i=0;i<rods.size();i++)
	{	
	  out << rods[i].getcomx() << '\t' << rods[i].getcomy() << '\t' <<rods[i].getcomz()<<'\n';
	}

      out.close();    
    }
  else if(mode==2)
    {// plot system as vectors

      os<<"../nodenetworkout/amorphon_"<<stra<<"traj_"<<extra<<"VEC_xyz.txt";;

      ofstream out(os.str());

      for(int i=0;i<rods.size();i++)
	{	
	  out << rods[i].botx() << '\t';
	  out << rods[i].boty() << '\t';
	  out << rods[i].botz() << '\t';
	  
	  out << rods[i].PBCx(rods[i].topx()-rods[i].botx()) << '\t';
	  out << rods[i].PBCy(rods[i].topy()-rods[i].boty()) << '\t';
	  out << rods[i].PBCz(rods[i].topz()-rods[i].botz()) << '\n';
	  
	}
      
      out.close();


    }

}


void cluster::plotcluster(int clust,int mode,int extra)
{

  ostringstream os;
  if(mode==0)
    {// plot beads in cluster
      os<<"../nodenetworkout/cluster_"<<clust<<"traj_"<<extra<<"beads_xyz.txt";

  ofstream out(os.str());

  for(int i=0;i<rods.size();i++)
    for(int j=0;j<rods[i].nbeads();j++){
      out << rods[i].getbeadx(j) << '\t' << rods[i].getbeady(j) << '\t' <<rods[i].getbeadz(j)<<'\n';
    }

  out.close();
    }
  else if(mode==1)
    {// plot COM of rods in cluster
      os<<"../nodenetworkout/cluster_"<<clust<<"traj_"<<extra<<"COM_xyz.txt";;

      ofstream out(os.str());

      for(int i=0;i<rods.size();i++)
	{	
	  out << rods[i].getcomx() << '\t' << rods[i].getcomy() << '\t' <<rods[i].getcomz()<<'\n';
	}

      out.close();    
    }
  else if(mode==2)
    {// plot system as vectors

      os<<"../nodenetworkout/cluster_"<<clust<<"traj_"<<extra<<"VEC_xyz.txt";;

      ofstream out(os.str());

      for(int i=0;i<rods.size();i++)
	{	
	  out << rods[i].botx() << '\t';
	  out << rods[i].boty() << '\t';
	  out << rods[i].botz() << '\t';
	  
	  out << rods[i].PBCx(rods[i].topx()-rods[i].botx()) << '\t';
	  out << rods[i].PBCy(rods[i].topy()-rods[i].boty()) << '\t';
	  out << rods[i].PBCz(rods[i].topz()-rods[i].botz()) << '\n';
	  
	}
      
      out.close();


    }
}

void cluster::addrodid(int id){

  rod_ids.push_back(id);

}

void cluster::setall_ids(int id){

  for(int i=0;i<rods.size()-1;i++)
    {
      rods[i].setcluster_id(id);
    }

}

void cluster::addrod(rod a){

  rods.push_back(a);

}

vector<rod> cluster::getrods(){

  return rods;

}

void cluster::blacklist(){

  blacklisted=true;
}


bool cluster::status(){
  return blacklisted;
}


cluster cluster::operator=(cluster a){

  this->rods=a.getrods();
  this->rod_ids=a.getrod_ids();
  return *this;

}

cluster cluster::operator+(cluster a){

  for(int i=0;i<a.csize();i++){
    rods.push_back(a.getrod(i));
    rod_ids.push_back(a.getrodid(i));
  }

  return *this;

}

vector<int> cluster::getrod_ids(){
  return rod_ids;
}

int cluster::getrodid(int i){

  return rod_ids[i];

}

bool cluster::operator==(cluster& a){
  

  int chk=0;
  for(int i=0;i<csize();i++)
    for(int j=0;j<a.csize();j++)
      {
	if(rods[i].thesame(rods[j])){
	  chk++;

	  break;}
      }

  if(chk==csize()-1)
    return true;
  else
    return false;
}


rod cluster::getrod(int i){
  return rods[i];
}


int cluster::rsize(){

  return rod_ids.size();
}
int cluster::csize(){

  return rods.size();
}

void bead::setpos(double x,double y,double z){
  
  this->x=x;
  this->y=y;
  this->z=z;
  
}

double bead::getx(){return x;}
double bead::gety(){return y;}
double bead::getz(){return z;}


void rod::belongstrand(){

  belongsstrand=true;

}

void rod::belongamorph(){
  
  belongsamorphon=true;

}

bool rod::getbelongstrand(){
  return belongsstrand;
}

bool rod::getbelongamorph(){
  return belongsamorphon;
}



bool rod::getbelong(){

return belongs;

}

void rod::addneigh(int rod_id){

  neighs.push_back(rod_id);

}

void rod::setcluster_id(int id){
  cluster_id=id;
}

void rod::setstrand_id(int id){
  strand_id=id;
}

void rod::setamorph_id(int id){
  amorphon_id=id;
}

int rod::mycluster_id(){

  return cluster_id;

}

int rod::myamorph_id(){
  return amorphon_id;
}

int rod::mystrand_id(){
  return strand_id;
}

void rod::belong(){

  belongs=true;

}

int rod::getneigh(){

  return neighs.size();

}

void rod::setbead(int beadid,double x,double y,double z){
  beads[beadid].setpos(x,y,z);
}

void rod::setBox(double Dx,double Dy,double Dz,bool xper,bool yper,bool zper){
  
  this->Dx=Dx;
  this->Dy=Dy;
  this->Dz=Dz;
  
  this->xper=xper;
  this->yper=yper;
  this->zper=zper;
  
}

double rod::PBCx(double x){
  
  if(xper)
    {
      while(x<=-Dx*0.5)
	{x+=Dx;}

      while(x>Dx*0.5){
	x-=Dx;
      }

      return x;
    }
  else
    return x;
}

double rod::PBCy(double y){
  
  if(yper)
    {
      while(y<=-Dy*0.5)
	{y+=Dy;}
      
      while(y>Dy*0.5){
	y-=Dy;
      }

      return y;
    }
  else
    return y;
}

double rod::PBCz(double z){
  
  if(zper)
    {
      while(z<=-Dz*0.5)
	{z+=Dz;}
      while(z>Dz*0.5){
	z-=Dz;
      }

      return z;
    }
  else
    return z;
}

bool rod::thesame(rod a){
  // check if two rods are the same by some arbitrary measure

  if(vecx()==a.vecx() && vecy()==a.vecy() && vecz()==a.vecz())
    return true;
  else
    return false;

}


bool rod::connected(rod b,double r){
  
  double R,x,y,z;
  
  for(int i=1;i<nbeads();i++)
    for(int j=b.nbeads()-1;j>=0;j--)
      {
	
	x=PBCx(getbeadx(i)-b.getbeadx(j));
	y=PBCy(getbeady(i)-b.getbeady(j));
	z=PBCz(getbeadz(i)-b.getbeadz(j));
	
	R=sqrt((x*x)+(y*y)+(z*z));
	
	if(R<=r)
	  {// Can leave function as you have shown the rods touch

	    /*double xc=PBCx(getcomx()-b.getcomx());
	    double yc=PBCy(getcomy()-b.getcomy());
	    double zc=PBCz(getcomz()-b.getcomz());
	    double comdiff=sqrt(xc*xc + yc*yc + zc*zc);
	    
	    cout<<"bead i: "<<i<<" and bead j: "<<j<<'\n';
	    cout<<"COM difference: "<<comdiff<<'\n';
	    */
	    return true;	    
	  }
      }
  return false; // haven't found any touching beads
  
}



double rod::dot(rod b){
  // computes dot product and divides by product of lengths
  double result;
  
  double x=(this->vecx()*b.vecx());
  double y=(this->vecy()*b.vecy());
  double z=(this->vecz()*b.vecz());
  
  result= x + y + z;
  result/=(this->getlen()*b.getlen());

  return result;  

}

double rod::getlen(){
  
  double x,y,z;
  
  x=vecx();
  y=vecy();
  z=vecz();
  
  return sqrt(((x*x)+(y*y)+(z*z)));
  
}

double rod::ang(rod b,int mode){
  
  double theta=acos(dot(b));;
  if(mode==0)
    return theta;
  else
    return theta*(180.0/M_PI);
}

void rod::setCOM(){
  
  xcom=ycom=zcom=0.0;
  
  for(int i=0;i<nbeads();i++)
    {
      xcom+=beads[i].getx();
      ycom+=beads[i].gety();
      zcom+=beads[i].getz();
    }
  
  xcom/=(double)nbeads();
  ycom/=(double)nbeads();
  zcom/=(double)nbeads();

}


double rod::getbeadx(int i){return beads[i].getx();}
double rod::getbeady(int i){return beads[i].gety();}
double rod::getbeadz(int i){return beads[i].getz();}


double rod::getcomx(){return xcom;}
double rod::getcomy(){return ycom;}
double rod::getcomz(){return zcom;}

double rod::vecx(){return PBCx(topx()-botx());}
double rod::vecy(){return PBCy(topy()-boty());}
double rod::vecz(){return PBCz(topz()-botz());}

double rod::topx(){return beads[nbeads()-1].getx();}
double rod::topy(){return beads[nbeads()-1].gety();}
double rod::topz(){return beads[nbeads()-1].getz();}

double rod::botx(){return beads[0].getx();}
double rod::boty(){return beads[0].gety();}
double rod::botz(){return beads[0].getz();}


void nodenetwork::do_classify(int traj){
  cout<<"Classifying rods\n";
  cout<<"# of rods: "<<rods.size()<<'\n';

  for(int i=0;i<rods.size()-1;i++)
    for(int j=i+1;j<rods.size();j++)
      {

	

	if(rods[i].connected(rods[j],touchperc)){
	
	  double angij=rods[i].ang(rods[j],1);
  
	  rods[i].addneigh(j);

	  
	  if( (angij >= (90.0 - anglethreshclust)) && ((90.0 + anglethreshclust ) >= angij))
	    {// Two rods do not align according to angle threshold
 
	      if(!rods[i].getbelong() && !rods[j].getbelong()){
		// two rods do not yet belong to a cluster
		
		
		cluster c(rods[i],rods[j]);
		clusters.push_back(c);
		
		clusters.back().addrodid(i);
		clusters.back().addrodid(j);
		
		rods[i].belong();
		rods[j].belong();
		
		rods[i].setcluster_id(clusters.size()-1);
		rods[j].setcluster_id(clusters.size()-1);
		
	      }
	      else if(rods[i].getbelong() && !rods[j].getbelong()){
		// i belongs but j doesnt, add j to i
		
		
		
		int id=rods[i].mycluster_id();
		
		clusters[id].addrod(rods[j]);
		clusters[id].addrodid(j);
		rods[j].belong();
		rods[j].setcluster_id(id);
		
		
	      }
	      else if(!rods[i].getbelong() && rods[j].getbelong()){
		// j belongs but i doesnt, add i to j
		
		
		int id=rods[j].mycluster_id();
		clusters[id].addrod(rods[i]);
		clusters[id].addrodid(i);
		
		rods[i].belong();
		rods[i].setcluster_id(id);
		
		
		
	      }
	      else if(rods[j].getbelong() && rods[i].getbelong()){// both belong so join the clusters
		
		int idi=rods[i].mycluster_id();
		int idj=rods[j].mycluster_id();
		
		if(clusters[idi]==clusters[idj])
		  {
		    
		    clusters[idi].blacklist();
		    continue;
		  }
		
		if(idi==idj)
		  {
		    continue;
		  }   
		
		
		
		
		cluster added;
		
		for(int h=0;h<clusters[idi].csize();h++)
		  {
		    added.addrodid(clusters[idi].getrodid(h));
		    added.addrod(clusters[idi].getrod(h));
		  }
		
		for(int h=0;h<clusters[idj].csize();h++)
		  {
		    added.addrodid(clusters[idj].getrodid(h));
		    added.addrod(clusters[idj].getrod(h));
		  }
		
		clusters[idi].blacklist();
		clusters[idj].blacklist();
		
		clusters.push_back(added);
		
		for(int k=0;k<clusters.back().rsize();k++)
		  {
		    int rodid=clusters.back().getrodid(k);
		    rods[rodid].setcluster_id(clusters.size()-1);
		  }
		
		
	      }
	      
	    }
	  else if(angij <= anglethreshstrand || angij >= (180-anglethreshstrand) )
	    {// rods do align


	      if(!rods[i].getbelongstrand() && !rods[j].getbelongstrand()){
		// two rods do not yet belong to a cluster
		
		
		strand c(rods[i],rods[j]);
		strands.push_back(c);
		
		strands.back().addrodid(i);
		strands.back().addrodid(j);
		
		rods[i].belongstrand();
		rods[j].belongstrand();
		
		rods[i].setstrand_id(strands.size()-1);
		rods[j].setstrand_id(strands.size()-1);
		
	      }
	      else if(rods[i].getbelongstrand() && !rods[j].getbelongstrand()){
		// i belongs but j doesnt, add j to i
		
		
		
		int id=rods[i].mystrand_id();
		
		strands[id].addrod(rods[j]);
		strands[id].addrodid(j);
		rods[j].belongstrand();
		rods[j].setstrand_id(id);
		
		
	      }
	      else if(!rods[i].getbelongstrand() && rods[j].getbelongstrand()){
		// j belongs but i doesnt, add i to j
		
		

		int id=rods[j].mystrand_id();
		
		strands[id].addrod(rods[i]);
		strands[id].addrodid(i);
		rods[i].belongstrand();
		rods[i].setstrand_id(id);
				
		
	      }
	      else if(rods[j].getbelongstrand() && rods[i].getbelongstrand()){// both belong so join the clusters
		
		int idi=rods[i].mystrand_id();
		int idj=rods[j].mystrand_id();
		
		if(strands[idi]==strands[idj])
		  {// gets rid of identical strands
		    
		    strands[idi].blacklist();
		    continue;
		  }
		
		if(idi==idj)
		  {
		    continue;
		  }   
		
		
		
		
		strand added;
		
		for(int h=0;h<strands[idi].csize();h++)
		  {
		    added.addrodid(strands[idi].getrodid(h));
		    added.addrod(strands[idi].getrod(h));
		  }
		
		for(int h=0;h<strands[idj].csize();h++)
		  {
		    added.addrodid(strands[idj].getrodid(h));
		    added.addrod(strands[idj].getrod(h));
		  }
		
		strands[idi].blacklist();
		strands[idj].blacklist();
		
		strands.push_back(added);
		
		for(int k=0;k<strands.back().rsize();k++)
		  {
		    int rodid=strands.back().getrodid(k);
		    rods[rodid].setstrand_id(strands.size()-1);
		  }
		
		
	      }



	    }
	  else
	    {// Create an amorphon - an amorhpous cluster

	      /*
	      if(!rods[i].getbelongamorph() && !rods[j].getbelongamorph()){
		// two rods do not yet belong to an amorphon
		
		
		amorphon c(rods[i],rods[j]);
		amorphons.push_back(c);
		
		amorphons.back().addrodid(i);
		amorphons.back().addrodid(j);
		
		rods[i].belongamorph();
		rods[j].belongamorph();
		
		rods[i].setamorph_id(amorphons.size()-1);
		rods[j].setamorph_id(amorphons.size()-1);
		
	      }
	      else if(rods[i].getbelongamorph() && !rods[j].getbelongamorph()){
		// i belongs but j doesnt, add j to i
		
		
		
		int id=rods[i].myamorph_id();
		
		amorphons[id].addrod(rods[j]);
		amorphons[id].addrodid(j);
		rods[j].belongamorph();
		rods[j].setamorph_id(id);
		
		
	      }
	      else if(!rods[i].getbelongamorph() && rods[j].getbelongamorph()){
		// j belongs but i doesnt, add j to i
		
		
		
		int id=rods[j].myamorph_id();
		
		amorphons[id].addrod(rods[i]);
		amorphons[id].addrodid(i);
		rods[i].belongamorph();
		rods[i].setamorph_id(id);
		
		
	      }
	      else if(rods[j].getbelongamorph() && rods[i].getbelongamorph()){// both belong so join the clusters
		
		int idi=rods[i].myamorph_id();
		int idj=rods[j].myamorph_id();
		
		if(amorphons[idi]==amorphons[idj])
		  {
		    
		    amorphons[idi].blacklist();
		    continue;
		  }
		
		if(idi==idj)
		  {
		    continue;
		  }   
		
		
		
		
		amorphon added;
		
		for(int h=0;h<amorphons[idi].csize();h++)
		  {
		    added.addrodid(amorphons[idi].getrodid(h));
		    added.addrod(amorphons[idi].getrod(h));
		  }
		
		for(int h=0;h<amorphons[idj].csize();h++)
		  {
		    added.addrodid(amorphons[idj].getrodid(h));
		    added.addrod(amorphons[idj].getrod(h));
		  }
		
		amorphons[idi].blacklist();
		amorphons[idj].blacklist();
		
		amorphons.push_back(added);
		
		for(int k=0;k<amorphons.back().rsize();k++)
		  {
		    int rodid=amorphons.back().getrodid(k);
		    rods[rodid].setamorph_id(amorphons.size()-1);
		  }
		
		
	      }
	      */
	    }
	  

	  
	  
	  
	    }
	
	
      }
  
  


  

  int clust=0;
  int nrods=0;
  for(int i=0;i<clusters.size();i++)
    {
      if(!clusters[i].status())
	{
	  //clusters[i].plotcluster(i,1,traj);
	  clusters[i].plotcluster(i,0,traj);
	  //clusters[i].plotcluster(i,2,traj);
	  nrods+=clusters[i].csize();
	  clust++;
	}	
    }

  cout<<"clusters: "<<clust<<'\n';
  clusters_out<<traj<<'\t'<<clust<<'\n';
  cout<<"rods: "<<nrods<<'\n';

  int stra=0;
  nrods=0;
  for(int i=0;i<strands.size();i++)
    {
     
      if(strands[i].csize()<rodsthreshstrand)
	strands[i].blacklist();
     
      if(!strands[i].status())
	{
	  //strands[i].plotcluster(i,1,traj);
	  strands[i].plotcluster(i,0,traj);
	  //strands[i].plotcluster(i,2,traj);
	  nrods+=strands[i].csize();
	  stra++;
	}	
    }

  cout<<"strands: "<<stra<<'\n';
  strands_out<<traj<<'\t'<<stra<<'\n';
  cout<<"rods: "<<nrods<<'\n';
  

  int amor=0;
  nrods=0;
  for(int i=0;i<amorphons.size();i++)
    {
      if(!amorphons[i].status())
	{
	  //amorphons[i].plotcluster(i,1,traj);
	  amorphons[i].plotcluster(i,0,traj);
	  //amorphons[i].plotcluster(i,2,traj);
	  nrods+=amorphons[i].csize();
	  amor++;
	}	
    }
  
  cout<<"amorphons: "<<amor<<'\n';
  amorphons_out<<traj<<'\t'<<amor<<'\n';
  cout<<"rods: "<<nrods<<'\n';  
}
