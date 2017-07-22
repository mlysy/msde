#ifndef myCobj_h
#define myCobj_h

class myCobj {
  double x;
  double z;
public:
  myCobj() {
    x = 19.1;
    z = -2.3;
  }
  double eval(double y) {
    return x + z * y;
  }
  double get_x(void) {return x;}
  void set_x(double _x) {
    x = _x;
    return;
  }
  ~myCobj() {
    Rprintf("objB destroyed.\n");
  }
};

#endif
