.. _cpp_reminder::


C++ Refresher
=============

C++ is the basic language of Siconos input files, thus you need at least a basic knowledge of C++.
This topic is covered by many books or tutorials, try for example :cite:`eckel_cpp`.
However, you can find in this first section the main C++ commands you need to know to write a first input file.

This page presents some basic commands and things you need to know to write a first C++ driver for your simulation.

   Note: In recent versions it is possible to design complete Siconos
   simulations using the Python programming language.  If you are
   using the Python interface, most of this document can be ignored.
   However, the Python interface uses the same class names.  See the
   :ref:`siconos_python_reference`.

Building/destroying and using objects
-------------------------------------

.. highlight:: c++

To describe your problem, you will mainly need to build pointers to objects.
Siconos use smart pointers so you do not have to bother with deletion.
The namespace SP enclose typedefs of all Siconos classes as smart pointers.

So to build an object of class *CHILD* derived from base class *BASE*, the syntax will be::

  SP::BASE my_object_name(new CHILD(some args ...)

For example, suppose you want to build a LagrangianDS, which belongs
to the base class DynamicalSystem, then do::

  SP::DynamicalSystem nameOfMyDS(new LagrangianDS(someParameters));

*new* call will reserve memory and build the object.

For each object, different types of constructors may be available. For an
exhaustive list see :ref:`siconos_api_reference`.

Members/Methods access
----------------------

The access to the methods of the objects is done thanks to "*" or "->"::
  
  nameOfMyDS->getType(); // return the type of the DynamicalSystem
  (*nameOfMyDS).getType(); // Same thing. 

Matrices and vectors handling
-----------------------------

The basic classes available in the platform to deal with vectors and matrices are:

* :doxysiconos:`SiconosVector` : vector of double, can be dense or sparse.
* :doxysiconos:`BlockVector` : vector of :doxysiconos:`SiconosVector`
* :doxysiconos:`SimpleMatrix` : matrix of double, can be dense, sparse, triangular, banded, symmetric, zero or identity.
* :doxysiconos:`BlockMatrix` : matrix of :doxysiconos:`SimpleMatrix` 

All these objects are just an interface to `Boost Ublas library <http://www.boost.org/libs/numeric/ublas/doc/index.htm>`_ vector and matrix. 

Notice that BlockVector or BlockMatrix are no more that a collection of pointers to SiconosVector or SimpleMatrix.
Then in most cases, to build such an object, you just need to insert some existing objects.
The usual ways of construction are described below.

::

   # build a dense vector of 4 elements
   SP::SiconosVector v1(new SimpleVector(4));

   # build a sparse vector of size 4
   SP::SiconosVector vsparse(new SimpleVector(4, Siconos::SPARSE)

   # build a block vector, which contains 3 dense blocks of size 5
   SP::BlockVector vblock(new BlockVector(3, 5));

   # build a block vector, which contains 2 sparse
   SP::SiconosVector v2(new SiconosVector(4, Siconos::SPARSE)
   SP::SiconosVector v3(new SiconosVector(7, Siconos::SPARSE)
   SP::BlockVector vblock(new BlockVector(v2, v3));

   int row = 3, col = 3;
   // row X col Dense matrix:
   SP::SiconosMatrix m(new SimpleMatrix(row,col));
   // row X col matrix, all elements initialized with a scalar value:
   double a = 4.4;
   SP::SiconosMatrix m(new SimpleMatrix(row,col,a));
   // row X row Symmetric matrix:
   SP::SiconosMatrix m(new SimpleMatrix(row,row, Siconos::SYMMETRIC));
   // Read from a file
   SP::SiconosMatrix m2(new SimpleMatrix("mat.dat",1)); // 1: ascii, 0:binary
   // Build an empty vector and insert some existing vectors.
   SP::BlockVector V0(new BlockVector());
   // Pointer insertion: 
   V0->insertPtr(v1); 
   // V0 has now one block equal to v1.
   // warning: because of pointer equality, 
   // v1 and (*V0)[0] represent the same object
   // and thus have the same memory location.
   // Copy of an existing vector:
   V0->insert(*v2); 
   // A new block has been created in V0
   // and v2 has been copied into this block.
   // Thus v2 and (*V0)[1] contain the same 
   // elements but are two different objects.

Note that a BlockVector can also contain some other BlockVector::

  SP::BlockVector V1(new BlockVector());
  V1->insertPtr(V0);
  V1->insertPtr(v1);
  
V1 has now two blocks: the first one is a block of two blocks and the second is equal to v1.

::

   // m1 ... m4 some SP::SiconosMatrix
   SP::SiconosMatrix M(new BlockMatrix(m1,m2,m3,m4));
   // M is a 2X2 blocks matrix 
   // (first row: m1, m2, second: m3, m4).

   
Keywords for constructors, in Siconos namespace: DENSE (default), TRIANGULAR, SYMMETRIC, SPARSE, BANDED, ZERO, IDENTITY.


Check the complete list of available constructors in reference documentation of each class.


Read/write vectors and matrices from/to file
""""""""""""""""""""""""""""""""""""""""""""

This is done using :doxysiconos:`ioVector` and :doxysiconos:`ioMatrix` classes.

::

   // Read/write vector/matrix from/to file
   // v is a vector, m a matrix
   ioVector myOutput ("MyData","ascii"); 
   myOutput.read(v); // read v from file MyData
   ioMatrix myMat("outMat","ascii");
   myMat.write(m); // Write m in file outMat
   
Input/Ouput Files format:

On the first line, the dimensions, with space as separator. Then the data. 

Example, for a 2-rows, 3-columns matrix:

::

   2 3
   1 2 3
   4 5 6

However, if you give as a second argument to write function "noDim", the first line with dimensions will not be written.

Methods and operations on matrices and vectors
""""""""""""""""""""""""""""""""""""""""""""""

Important note: in many of the operators described below, a boolean argument "init" can be set. If equal to true (default value) then the operator used "=" and if set to 
false, "+=".

::
   
   v->size() // return the size of the vector
   m->size(0); // number of rows in the matrix
   m->size(1), // number of columns
   m->resize(a,b); // resize m, available also for vectors

   // To compute C = A*B
   prod(A,B,C,true);
   // or
   prod(A,B,C);

   // To compute C += A*B
   prod(A,B,C,false);

   //Single elements access or assignment: operator "()" or \e get/setValue functions.
   SP::SiconosVector v(new SimpleVector(3)); // v = [0 0 0]
   SimpleVector w(4);  			 // w = [0 0 0 0]
   (*v)(0) = 4;				 // v = [4 0 0] 
   // equivalent to:
   v->setValue(0,4); 
   w(1) = 2;
   w(2) = (*v)(0);				 // w = [0 2 4 0]
   // equivalent to:
   w.setValue( 2,v->getValue(0) );

   SP::SiconosMatrix M(new SimpleMatrix(3,3)); // M = [ 0 0 0 ]
                                               //     [ 0 0 0 ]
					       //     [ 0 0 0 ]
   SimpleMatrix P(1,2);   		       // P = [ 0 0 ]

   (*M)(1,2) = 2; 
   P(0,1) = 12;				   // P = [ 0 12.0 ]
   M->setValue(2,0,3.6);		   // M = [  0  0  0  ]
					   //     [  0  0 2.0 ]
					   //     [ 3.6 0  0  ]
	
   cout << P.getValue(0,1); // display 12.0

Note: for sparse matrices, assignment with operator "()" fails. It is then necessary to use setValue function.

::
   
   SP::SiconosMatrix A(new SimpleMatrix(10,10,SPARSE));
   (*A)(0,0) = 12; // WRONG
   A->setValue(0,0,12); // OK

For BlockVector: "()" and get/setValue functions have the same action as for SimpleVectors::

  // We suppose that v1 and v2 are two pointers to SimpleVector of size 3 and 4.
  SP::SiconosVector vB(new BlockVector(v1,v2)); // vB = [ [1 2 3] [4 5 6 7] ]
  (*vB)(4) = 12; 				      // vB = [ [1 2 3] [4 12 6 7] ]
  vB->setValue(6,8.6); 		              // vB = [ [1 2 3] [4 12 6 8.6] ]	
  // Warning: the given input for position is an "absolute" one, not a block position.

Remark: get/setValue functions are equivalent to "()" operator but mainly useful in Siconos-Python, since in that case operators can not be overloaded and thus
"()" is not available. The same remark applies for "[ ]" get/setVector and in a general way for all operators overloading.

::

   // Set vector or matrix to zero or identity
   x->zero();
   A->zero();
   A->eye(); 

   // Assignment of vectors or matrices: "A = B" or "x = y"
   // Operator =
   // Ok if A and x have been built before.
   A = B;
   x = y;
   // Remark: sizes must be consistents between A/B and x/y, 
   // else it results in a Siconos Exception.

   // Else copy constructor: memory allocation and initialization with B or x
   SP::SiconosMatrix A(new SimpleMatrix(*B));
   SP::SiconosVector x(new SimpleVector(*y));

   // Addition of matrices or vectors

   // add "in place": A = A+B  or x = x+y
   A += B;
   x += y;
   
   C = A+B;
   add(A,B,C);
   A -= B;
   C = A-B;
   sub(A,B,C);

   // Multiplication by a scalar:
   A *=a;
   B = a*A;
   scal(a,A,B);
   A /=a;
   x /=a;
   B = A/a;
   scal(1.0/a,A,B);
   // matrices product
   C = A*B;
   prod(A,B,C); // Based on atlas gemm for Dense matrices and ublas::prod for others. 
                // C and A or B can be the same matrices (ie have common memory), 
	        // but that will slow down the operation.
   gemm(A,B,C); // Only for denses matrices.

   // It is also possible to compute product of sub-blocks of matrices or vectors:
   // Declare A, x, y ...
   // 
   std::vector<unsigned int> coord;
   // Set coord values ...
   bool init = false;
   subprod(A,x,y,coord,init);


Coord vector is equal to [r0A, r1A,, c0A, c1A, r0x, r1x, r0y, r1y]. The sub-matrix A is the matrix between row positions
r0A and r1A, column position between c0A and c1A. Same thing for x and y with rix, riy.
Then subprod computes suby = subA*subx if init = true, or suby += subA*subx if init = false.


::
   
   Matrix transpose:
   // in place:
   A->trans();
   // B = At
   B->trans(A);

   // inner product: a = x.y
   a = inner_prod(x,y);

   // Matrix-vector product: \f$y=A*x\f$
   y = prod(A,x);
   prod(A,x,y);

To handle a specific block, use "[ ]" or getVector and getVectorPtr functions::

  SP::SiconosVector v3(new SimpleVector(3));  // v3 = [0 0 0]
  SP::SiconosVector v4(new SimpleVector(4));  // v4 = [0 0 0 0]
  // get and copy a block:
  *v3 = *(*vB[0]); 			   // v3 = v1 = [1 2 3]
  // Equivalent to
  *v3 = *vB->getVectorPtr(0);

  // get and copy pointer to block:
  v4 = vB->getVectorPtr(1);		   // v4 = v2 = [4 12 6 8.6]
					   // AND pointer equality 
					   // between v4, vB[1] and v2
  // Equivalent to:
  v4 = (*vB)[1];			           // v4 = v2 = [4 12 6 8.6]. 

  // Assignment:
  SP::SiconosVector v5(new SimpleVector(3));  // v5 = [0 0 0]
  
  *(*vB)[0] = *v5; //  vB = [ [0 0 0] [4 5 6 7] ]
                   //  AND v1 = [0 0 0] because of pointer link between vB[0] and v1.
  // Equivalent to:
  vB->setVector(0,*v5);

  (*v5)(1) = 12;
  vB->setVectorPtr(0,v5); // vB = [ [0 0 0] [0 12 0] ]
  // Pointer equality between v5 and vB[0]. 
  // The pointer link between vB[0] and v1 has been canceled.

  // Warning: when using setVectorPtr(i,w), 
  // the vector w must be of the same size as the block[i] of v. 

About efficiency
""""""""""""""""

As you can see above, for most functionality, two solutions are available: either an overloaded operator or a function without any return value.
For example in the case of matrix addition::

  C = A + B;
  // or 
  add(A,B,C);

In a general way, if you need efficiency, always prefer functions to overloaded operators. 
The first solution is just there to give a more pleasant and readable, way of writing operations.

Try also to use pointers to objects to avoid temporary and time-consuming copies.
