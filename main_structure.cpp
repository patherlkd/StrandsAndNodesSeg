#include "main.h"
#include "nodenetwork.h"

/* command line arguments

1st = angle for aligned rods
2nd = angle for non-aligned rods
3rd = compute clusters every "this number" of intervals
4th = maximum trajctory to analyse

*/

int main(int argc, char* argv[] ){

  string xyzname=argv[5];

  nodenetwork test(xyzname);

  test.add_particletype(1);
  test.add_particletype(2);


  double dia = 1.12;  //bead diameter
  double rodthresh=5; // threshold number of rods in a strand (e.g. need at least 5)
  double length=35; // box length
  double touchperc=1.0; // 100% of diameter
  double ang=atof(argv[1]); // in degrees. vectors @ < 'ang' and 180 - 'ang' < @ < 180  are aligned
  double ang2=atof(argv[2]); // to sort non-aligned rods

  if(ang< 0 || ang2 < 0 || ang > 180 || ang2 > 180)
    {
      cout << "angles must be  0 < angle < 180\n";
      exit(1);
    }

  int inter=atoi(argv[3]); // do every inter trajectory
  int maxtraj=atoi(argv[4]); // max number of trajectories

  if(inter>=maxtraj)
    {
      cout<<"inter must be less than max trajectories (maxtraj)\n";
      exit(1);
    }

  test.set_dia(dia); // set bead diameter
  test.setbox(-length,length,-length,length,-length,length);
  test.setperiodic(true,true,true); // periodic boundaries have been used? ( x , y , z ) set to all false for unwrapped coordinates

  test.set_rodsthreshstrand(rodthresh);
  test.set_touchperc(touchperc);
  test.set_anglethreshstrand(ang);
  test.set_anglethreshclust(ang2);

  // Main program! 
  test.do_network(inter,maxtraj);

  return 0;

}
