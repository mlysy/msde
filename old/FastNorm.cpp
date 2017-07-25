#include <R_ext/Random.h>
#include <math.h>

static bool doCalc = true;
static double yy1, yy2;
static double unival, norval;

void user_unif_init(Int32 seed_in) {
  doCalc = true;
  return;
}

double * user_unif_rand() {
  unival = 0.0;
  return &unival;
}

extern "C" double *user_norm_rand() {
  double x1, x2, w;
  if(doCalc) {
    do {
      x1 = 2.0 * unif_rand() - 1.0;
      x2 = 2.0 * unif_rand() - 1.0;
      w = x1 * x1 + x2 * x2;
    } while ( w >= 1.0 );
    w = sqrt( (-2.0 * log( w ) ) / w );
    yy1 = x1 * w;
    yy2 = x2 * w;
    norval = yy1;
    doCalc = false;
  }
  else {
    norval = yy2;
    doCalc = true;
  }
  return &norval;
}
