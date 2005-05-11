
...

public:

// --- GETTERS AND SETTERS ---

// -- for non pointers members --
// get objectMemberA
inline const typeObject getObjectMemberA() const
{
  return objectMemberA;
}

// assign newValue to objectMemberA
inline void setObjectMemberA(const typeObject& newValue)
{
  objectMemberA = newValue;
}

// return a pointer on objectMemberA (optional)
inline const typeObject* getObjectMemberAPtr() const
{
  return &objectMemberA;
}

// -- for pointer type members --

// a - get the value of pointed object
inline const typeObject getObjectMemberB() const
{
  return *objectMemberB;
}

// b - get the pointer
inline typeObject* getObjectMemberBPtr() const
{
  return objectMemberB;
}

// c - set the value of the pointed object
inline void setObjectMemberB(const typeObject& newValue)
{
  *objectMemberB = newValue;
}

// d - set the pointer
inline void setObjectMemberBPtr(typeObject* newPtr)
{
  delete objectMemberB;
  objectMemberB = 0;
  objectMemberB = newPtr;
}

// --- PRIVATE/PROTECTED MEMBERS ---

protected:

typeObject ObjectMemberA ;

typeObject *ObjectMemberB;

