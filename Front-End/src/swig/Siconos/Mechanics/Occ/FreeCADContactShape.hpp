xc#ifndef FreeCADContactShape_hpp
#define FreeCADContactShape_hpp

struct FreeCADContactShape : public OccContactShape
{
  FreeCADContactShape(PyObject *o)
  {


  /* ptr to FreeCAD TopoShapePy */
  PyObject * ref;

  virtual void move(const SiconosVector& q)
  {
    this->OccContactShape::move(q);
    
  }
}

#endif
