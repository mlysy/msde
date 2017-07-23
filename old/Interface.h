// how to interface R and C objects

// Abstract Base class Cobj
// derived class template Robj
// Robj has private members defined at compile-time.
// Robj has virtual public members (from Cobj) aliasing private members,
// called at run-time

#ifndef Interface_h
#define Interface_h

class Cobj {
 public:
  virtual double eval(double y) = 0;
  virtual double get_x(void) = 0;
  virtual void set_x(double x) = 0;
  virtual ~Cobj() = 0;
};

inline Cobj::~Cobj() {
  Rprintf("Cobj destroyed.\n");
}

template <class T>
class Robj : public Cobj {
  T *foo;
 public:
  virtual double eval(double y);
  virtual double get_x(void);
  virtual void set_x(double x);
  virtual ~Robj() {
    delete foo;
    Rprintf("Robj destroyed.\n");
  }
  Robj() {
    foo = new T;
  }
};

template <class T>
inline double Robj<T>::eval(double y) {
  return foo->eval(y);
}

template <class T>
inline double Robj<T>::get_x(void) {
  return foo->get_x();
}

template <class T>
inline void Robj<T>::set_x(double x) {
  foo->set_x(x);
  return;
}


#endif
