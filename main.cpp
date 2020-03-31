#include "main.h"
#include "pixelnetwork.h"

/* DEPRECATED.. use only if you understand pixelnetwork .h and .cpp */

int main( ){

  trajectory test("/storage/users/ldavis/Simulation/mol_dyn/Collagen/xyz/kick_out.xyz");
  
  test.set_xyzgridname("/storage/users/ldavis/Simulation/mol_dyn/Collagen/xyz/grid.xyz");

  double length=56;
  double pixlen=8.0;
  int inter=5;
  int maxtraj=87;
  double thresh=0.0;

  test.set_dia(1.12);
  test.setbox(-length,length,-length,length,-length,length);
  test.setperiodic(true,true,true);
  test.set_pixl(pixlen);
  test.set_threshdens(thresh);

  test.do_grid(inter,maxtraj);

  return 0;

}
