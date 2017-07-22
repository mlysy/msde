#ifndef myCobj_h
#define myCobj_h

class myCobj {
  double x;
public:
  myCobj() {
    x = 3.5;
  }
  double eval(double y) {
    return x + y;
  }
  double get_x(void) {return x;}
  void set_x(double _x) {
    x = _x;
    return;
  }
  ~myCobj() {
    Rprintf("objA destroyed.\n");
  }
};

#endif
