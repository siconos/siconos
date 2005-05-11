//=============== BASE CLASS =============================

// include ...
// --- CONSTRUCTORS ---

// From XML
BaseClass::BaseClass(DSXML * dsXML): objectMemberA(initValue), sizeB(0), objectMemberB(0), ...
{
  // xml data loading
  sizeB = ...;
  // ...
  // Memory allocation
  objectMemberB = new typeObject(sizeB) ;
}

// From a minimum set of data
BaseClass::BaseClass(typeObject A, typeObject* B, int size): objectMemberA(A), sizeB(size), objectMemberB(0), ...
{
  // Memory allocation
  objectMemberB = new typeObject(sizeB) ;
  // Assign B value to objectMemberB
  // warning! Do not set objectMemberB = B: should cause segmentation error when call to destructor if B has been declared with a new.
  *objectMemberB = *B ;
  // ...
}

// --- DESTRUCTOR ---
BaseClass::~BaseClass()
{
  delete objectMemberB;
  objectMemberB = 0;
}

// --- Default constructor ---
BaseClass::BaseClass(): objectMemberA(0), sizeB(0), objectMemberB(0), ...
{}


//=============== DERIVED CLASS =============================

// warning !! Implemented in a different file !!
// include ...
// --- CONSTRUCTORS ---

// From XML
DerivedClass::DerivedClass(DSXML * dsXML): BaseClass(dsXML), objectMemberC(0), sizeC(0) ...
{
  // xml data loading
  sizeC = ...;
  // ...
  // Memory allocation
  objectMemberC = new typeObject(sizeC) ;
}

// From a minimum set of data
DerivedClass::DerivedClass(typeObject* C, typeObject A, typeObject* B, int size): BaseClass(A, B, size), objectMemberC(C), sizeC(size) ...
{
  objectMemberC = new typeObject(sizeC) ;
  *objectMemberC = *C ;
  // ...
}

// --- DESTRUCTOR ---
DerivedClass::~DerivedClass()
{
  delete objectMemberC;
  objectMemberC = 0;
  // Warning: do not delete baseClass members, the base destructor is automatically called.
}

// ...

// --- Default constructor ---
DerivedClass::DerivedClass(): BaseClass(), objectMemberC(0), sizeC(0)...
{}

