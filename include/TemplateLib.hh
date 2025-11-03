/*
  TemplateLib.hh

  2012/1/24 
*/

#ifndef TemplateLib_h
#define TemplateLib_h 1

struct DeleteObject {
  template <typename T>
  void operator()(T *ptr) const
  { delete ptr; }
}; 

#endif
