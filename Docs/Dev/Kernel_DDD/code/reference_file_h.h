//=============== BASE CLASS =============================
#ifndef CLASS_NAME_H
#define CLASS_NAME_H

// include ...

class BaseClassName

{
public:

  // --- CONSTRUCTORS ---

  // From XML
  BaseClass::BaseClass(DSXML * dsXML);

  // From a minimum set of data
  BaseClass::BaseClass(typeObject A, typeObject* B);

  // --- DESTRUCTOR ---
  virtual BaseClass::~BaseClass();

  // --- GETTERS AND SETTERS ---

  // get sizeB
  inline const int getSizeB() const
  {
    return sizeB;
  }

  // assign newValue to sizeB
  inline void setSizeB(const int& newValue)
  {
    sizeB = newValue;
  }

  // get the value of objectMemberB
  inline const Type getObjectMemberB() const
  {
    return *objectMemberB;
  }

  // get the pointer
  inline Type* getObjectMemberBPtr() const
  {
    return objectMemberB;
  }

  // set the value of the pointed object
  inline void setObjectMemberB(const Type& newValue)
  {
    *objectMemberB = newValue;
  }

  // set the pointer
  inline void setObjectMemberBPtr(Type* newPtr)
  {
    delete objectMemberB;
    objectMemberB = 0;
    objectMemberB = newPtr;
  }
  // ...

  // --- OTHER FUNCTIONS ---
  // able to change data members
  returnType func1(...) ;

  // which is not supposed to modify data members
  returnType func2(...) const ;

  // which can not modify input arguments
  returnType func3(const typeInput& inputArg);

  // with return argument that can not be modified => can not be set as an lValue
  const returnType func4(...);

  // and combinations of these previous cases ...

  // --- PRIVATE/PROTECTED MEMBERS ---

protected:

  // -- Default constructor --
  BaseClass::BaseClass();

  // -- data members --
  typeObject ObjectMemberA ;

  typeObject *ObjectMemberB;

  int sizeB;

  // ...

};

#endif

//=============== DERIVED CLASS =============================

// warning !! Implemented in a different file !!

#ifndef DERIVED_CLASS_NAME_H
#define DERIVED_CLASS_NAME_H

// include ...

class DerivedClassName: public BaseClassName
{
public:

  // --- CONSTRUCTORS ---
  // From XML
  DerivedClass::DerivedClass(DSXML * dsXML);
  // From a minimum set of data
  DerivedClass::DerivedClass(typeObject A, typeObject* B);
  // --- DESTRUCTOR ---
  DerivedClass::~DerivedClass();

  // --- GETTERS AND SETTERS ---
  // ...

  // --- OTHER FUNCTIONS ---
  // ...

  // --- PRIVATE/PROTECTED MEMBERS ---

protected:

  // -- Default constructor --
  DerivedClass::DerivedClass();

  // -- data members --
  typeObject *ObjectMemberC ;

  int sizeC;

  // ...

};

#endif
