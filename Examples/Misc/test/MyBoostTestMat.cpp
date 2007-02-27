#include "ioMatrix.h"


#include <iostream>
#include <vector>



int main()
{
  const std::vector<double> v0(4, 0);
  const std::vector<double> v1(4, 1);
  const std::vector<double> v3(2, 3);
  const std::vector<double> v5(25, 0);
  const std::vector<double> v6(6, 0);
  const std::vector<double> v(6, 1);

  const scalar_matrix<double> sc1(2, 2, 1);
  const scalar_matrix<double> sc0(2, 2, 0);
  mapped_matrix<double> m1(sc1, 4);
  mapped_matrix<double> m0(sc0, 4);



  MySimpleVector vect3(DENSE, v3, 2);
  MySimpleVector vect1(DENSE, v3, 2);
  MySiconosVector *vect10 = new MySimpleVector(DENSE, v3, 2);
  //DENSE
  MySiconosMatrix     *truc1   = new MySimpleMatrix(DENSE, v1, 2, 2);
  MySiconosMatrix     *truc1_2 = new MySimpleMatrix(DENSE, v1, 2, 2);
  MySiconosMatrix     *truc1_3 = new MySimpleMatrix(DENSE, v0, 2, 2);
  //TRIANGULAR
  MySiconosMatrix     *truc2   = new MySimpleMatrix(TRIANGULAR, v1, 2, 2);
  MySiconosMatrix     *truc2_2 = new MySimpleMatrix(TRIANGULAR, v1, 2, 2);
  //SYMMETRIC
  MySiconosMatrix     *truc3   = new MySimpleMatrix(SYMMETRIC, v1, 2);
  MySiconosMatrix     *truc3_2 = new MySimpleMatrix(SYMMETRIC, v1, 2);
  //SPARSE
  MySiconosMatrix     *truc4   = new MySimpleMatrix(m1);
  MySiconosMatrix     *truc4_2 = new MySimpleMatrix(m1);
  //BANDED
  MySiconosMatrix     *truc5   = new MySimpleMatrix(BANDED, v, 2, 2, 1, 1);
  MySiconosMatrix     *truc5_2 = new MySimpleMatrix(BANDED, v, 2, 2, 1, 1);

  MySimpleMatrix val1(DENSE, v1, 2, 2);
  MySimpleMatrix val1_2(DENSE, v1, 2, 2);
  MySimpleMatrix val1_3(DENSE, v0, 2, 2);

  MySimpleMatrix val2(TRIANGULAR, v1, 2, 2);
  MySimpleMatrix val2_1(TRIANGULAR, v0, 2, 2);
  MySimpleMatrix val2_2(TRIANGULAR, v0, 2, 2);

  MySimpleMatrix val3(SYMMETRIC, v1, 2);
  MySimpleMatrix val3_1(SYMMETRIC, v0, 2);
  MySimpleMatrix val3_2(SYMMETRIC, v0, 2);

  MySimpleMatrix val4(m1);
  MySimpleMatrix val4_1(m0);
  MySimpleMatrix val4_2(m0);

  MySimpleMatrix val5(BANDED, v, 2, 2, 1, 1);
  MySimpleMatrix val5_1(BANDED, v6, 2, 2, 1, 1);
  MySimpleMatrix val5_2(BANDED, v6, 2, 2, 1, 1);

  std::cout << "*truc1 ";
  (*truc1).display();

  std::cout << "*truc2 ";
  (*truc2).display();

  std::cout << "*truc3 ";
  (*truc3).display();

  std::cout << "*truc4 ";
  (*truc4).display();

  std::cout << "*truc5 ";
  (*truc5).display();

  std::cout << "val1 ";
  val1.display();

  std::cout << "val1_2 ";
  val1_2.display();

  std::cout << "val2 ";
  val2.display();

  std::cout << "val2_1 ";
  val2_1.display();

  std::cout << "val2_2 ";
  val2_2.display();

  std::cout << "val3 ";
  val3.display();

  std::cout << "val3_1 ";
  val3_1.display();

  std::cout << "val3_2 ";
  val3_2.display();

  std::cout << "val4 ";
  val4.display();

  std::cout << "val4_1 ";
  val4_1.display();

  std::cout << "val4_2 ";
  val4_2.display();

  std::cout << "val5 ";
  val5.display();

  std::cout << "val5_1 ";
  val5_1.display();

  std::cout << "val5_2 ";
  val5_2.display();

  std::cout << std::endl;

  /*
     val1 -= val1;
     std::cout<<"val1, val1-=val1";
     val1.display();
  */
  /************** OPERATIONS AVEC POINTEURS ******************/

  std::cout << "************************** operation = ***************************" << std::endl;
  *truc1 = *truc2;
  std::cout << "*truc1 = *truc2 ";
  (*truc1).display();

  *truc1 = *truc4;
  std::cout << "*truc1 = *truc2 ";
  (*truc1).display();

  *truc1 = *truc5;
  std::cout << "*truc1 = *truc5 ";
  (*truc1).display();

  *truc1 = *truc3;
  std::cout << "*truc1 = *truc3 ";
  (*truc1).display();
  std::cout << std::endl;
  std::cout << "*************************  operation /=  *************************" << std::endl;
  *truc1 /= 2;
  std::cout << "*truc1 /= 2 ";
  (*truc1).display();

  *truc2 /= 3;
  std::cout << "*truc2 /= 3 ";
  (*truc2).display();

  *truc3 /= 4;
  std::cout << "*truc3 /= 4 ";
  (*truc3).display();

  *truc4 /= 5;
  std::cout << "*truc4 /= 5 ";
  (*truc4).display();

  *truc5 /= 10;
  std::cout << "*truc5 /= 10 ";
  (*truc5).display();

  std::cout << std::endl;
  std::cout << "*************************  operation *=  *************************" << std::endl;
  *truc1 *= 2;
  std::cout << "*truc1 *= 2 ";
  (*truc1).display();

  *truc2 *= 3;
  std::cout << "*truc2 *= 3 ";
  (*truc2).display();

  *truc3 *= 4;
  std::cout << "*truc3 *= 4 ";
  (*truc3).display();

  *truc4 *= 5;
  std::cout << "*truc4 *= 5 ";
  (*truc4).display();

  *truc5 *= 10;
  std::cout << "*truc5 *= 10 ";
  (*truc5).display();

  std::cout << std::endl;
  std::cout << "*************************  operation -=  *************************" << std::endl;
  *truc1 -= *truc1;
  std::cout << "*truc1 -= *truc1 ";
  (*truc1).display();

  *truc1 -= *truc2;
  std::cout << "*truc1 -= *truc2 ";
  (*truc1).display();

  *truc1 -= *truc3;
  std::cout << "*truc1 -= *truc3 ";
  (*truc1).display();

  *truc1 -= *truc4;
  std::cout << "*truc1 -= *truc4 ";
  (*truc1).display();

  *truc1 -= *truc5;
  std::cout << "*truc1 -= *truc5 ";
  (*truc1).display();

  std::cout << std::endl;
  std::cout << "*************************  operation +=  *************************" << std::endl;
  *truc1 += *truc1;
  std::cout << "*truc1 += *truc1 ";
  (*truc1).display();

  *truc2 += *truc2;
  std::cout << "*truc2 += *truc2 ";
  (*truc2).display();

  *truc4 += *truc4;
  std::cout << "*truc4 += *truc4 ";
  (*truc4).display();

  *truc1 += *truc5;
  std::cout << "*truc1 += *truc5 ";
  (*truc1).display();

  *truc1 += *truc3;
  std::cout << "*truc1 += *truc3 ";
  (*truc1).display();
  std::cout << std::endl;
  std::cout << "*************** Restauration des Matrices initiales **************" << std::endl;
  *truc1 = *truc1_2;
  std::cout << "*truc1 ";
  (*truc1).display();

  *truc2 = *truc2_2;
  std::cout << "*truc2 ";
  (*truc2).display();

  *truc3 = *truc3_2;
  std::cout << "*truc3 ";
  (*truc3).display();

  *truc4 = *truc4_2;
  std::cout << "*truc4 ";
  (*truc4).display();

  *truc5 = *truc5_2;
  std::cout << "*truc5 ";
  (*truc5).display();

  std::cout << std::endl;
  std::cout << "******************** operation de Soustraction *******************" << std::endl;
  std::cout << "*truc2 - *truc3 " << std::endl;
  (sub(*truc2, *truc3)).display();

  std::cout << "*truc2 - *truc2 " << std::endl;
  (*truc2 - *truc2).display();

  std::cout << "*truc1 - *truc1 " << std::endl;
  (*truc1 - *truc1).display();

  std::cout << "*truc2 - *truc1 " << std::endl;
  (sub(*truc2, *truc1)).display();

  std::cout << "*truc4 - *truc5 " << std::endl;
  (sub(*truc4, *truc5)).display();

  std::cout << "*truc2 - *truc5 " << std::endl;
  (sub(*truc2, *truc5)).display();


  std::cout << "*truc5 - *truc5 " << std::endl;
  (*truc5 - *truc5).display();

  std::cout << std::endl;
  std::cout << "*********************** Somme de Matrices ************************" << std::endl;
  std::cout << "*truc2 + *truc3 " << std::endl;
  (add(*truc2, *truc3)).display();

  std::cout << "num de *truc2 + *truc3 " << (add(*truc2, *truc3)).GetNum() << std::endl;

  std::cout << "*truc2 + *truc2_2 " << std::endl;
  (*truc2 + *truc2_2).display();

  std::cout << "num de *truc2 + *truc2_2 " << (*truc2 + *truc2_2).GetNum() << std::endl;

  std::cout << "*truc1 + *truc3 " << std::endl;
  (add(*truc1, *truc3)).display();

  std::cout << "*truc2 + *truc1 " << std::endl;
  (add(*truc2, *truc1)).display();

  std::cout << "*truc5 + *truc5_2 " << std::endl;
  //(add(*truc2, *truc1)).display();
  (*truc5 + *truc5_2).display();

  std::cout << "add(*truc5, *truc5_2) " << std::endl;
  (add(*truc5, *truc5_2)).display();

  std::cout << std::endl;
  std::cout << "******************* Division avec un scalaire ********************" << std::endl;
  std::cout << "*truc1/2 " << std::endl;
  (*truc1 / 2.).display();

  std::cout << "*truc2/2. " << std::endl;
  (*truc2 / 2.).display();

  std::cout << "*truc3/3 " << std::endl;
  (*truc3 / 3.).display();

  std::cout << "*truc4/2. " << std::endl;
  (*truc4 / 2.).display();

  std::cout << std::endl;
  std::cout << "****************** Multiplication avec scalaire ******************" << std::endl;
  std::cout << "*truc1*2 " << std::endl;
  (*truc1 * 2.).display();

  std::cout << "*truc2*2. " << std::endl;
  (*truc2 * 2.).display();

  std::cout << "*truc3*3 " << std::endl;
  (*truc3 * 3.).display();

  std::cout << "3 * *truc3 " << std::endl;
  (3 * (*truc3)).display();

  std::cout << "*truc5*2. " << std::endl;
  (*truc5 * 2.).display();

  std::cout << std::endl;
  std::cout << "****************** Multiplication entre Matrices *****************" << std::endl;
  std::cout << "*truc2 * *truc2 ";
  (prod(*truc2, *truc2)).display();

  std::cout << "num de *truc2 * *truc2 " << (prod(*truc2, *truc2)).GetNum() << std::endl;

  std::cout << "*truc1 * *truc2 ";
  (prod(*truc1, *truc2)).display();

  std::cout << "num de *truc1 * *truc2 " << (prod(*truc1, *truc2)).GetNum() << std::endl;

  std::cout << "*truc2 * *truc3 " << std::endl;
  (prod(*truc2, *truc3)).display();

  std::cout << "num de *truc2 * *truc3 " << (prod(*truc2, *truc3)).GetNum() << std::endl;

  std::cout << "*truc5 * *truc5 ";
  (*truc5 * *truc5).display();

  std::cout << "*truc5 * *truc5 ";
  (prod(*truc5, *truc5)).display();

  std::cout << "*truc4 * *truc4 ";
  (*truc4 * *truc4).display();

  std::cout << std::endl;
  std::cout << "***************** Acces aux elements dune Matrice ****************" << std::endl;

  std::cout << "*truc1(0,0) " << (*truc1)(0, 0) << std::endl;

  std::cout << "*truc2(0,1) " << (*truc2)(0, 1) << std::endl;

  std::cout << "*truc3(1,0) " << (*truc3)(1, 0) << std::endl;

  std::cout << "*truc4(1,0) " << (*truc4)(1, 0) << std::endl;

  std::cout << "*truc5(1,0) " << (*truc5)(1, 0) << std::endl;

  std::cout << "val1(0,0) " << val1(0, 0) << std::endl;

  std::cout << "val2(0,1) " << val2(0, 1) << std::endl;

  std::cout << "val3(1,0) " << val3(1, 0) << std::endl;

  std::cout << "val4(1,0) " << val4(1, 0) << std::endl;

  std::cout << "val5(1,0) " << val5(1, 0) << std::endl;

  std::cout << std::endl;
  std::cout << "************************ Zero d une Matrice **********************" << std::endl;

  (*truc1).zero();
  std::cout << "*truc1 apres zero() " << std::endl;
  (*truc1).display();

  std::cout << std::endl;

  (*truc1).eye();
  std::cout << "*truc1 apres eye() " << std::endl;
  (*truc1).display();

  std::cout << std::endl;

  (*truc4).eye();
  std::cout << "*truc4 apres eye() " << std::endl;
  (*truc4).display();

  std::cout << std::endl;

  (*truc5).eye();
  std::cout << "*truc5 apres eye() " << std::endl;
  (*truc5).display();

  std::cout << std::endl;
  std::cout << "*************** Restauration des Matrices initiales **************" << std::endl;
  *truc1 = *truc1_2;
  std::cout << "*truc1 ";
  (*truc1).display();

  *truc2 = *truc2_2;
  std::cout << "*truc2 ";
  (*truc2).display();

  *truc3 = *truc3_2;
  std::cout << "*truc3 ";
  (*truc3).display();

  *truc4 = *truc4_2;
  std::cout << "*truc4 ";
  (*truc4).display();

  *truc5 = *truc5_2;
  std::cout << "*truc5 ";
  (*truc5).display();

  std::cout << std::endl;
  std::cout << "********************* transpose d une Matrice *******************" << std::endl;

  std::cout << "trans(*truc1).display() " << std::endl;
  trans(*truc1).display();
  std::cout << std::endl;
  std::cout << "********************* Row et Column d une Matrice *******************" << std::endl;

  std::cout << "Row de *truc1" << std::endl;
  (*truc1).GetRow(1, vect1);
  vect1.display();

  std::cout << "Row de *truc2" << std::endl;
  (*truc2).GetRow(1, vect1);
  vect1.display();

  std::cout << "Row de *truc5" << std::endl;
  (*truc5).GetRow(1, vect1);
  vect1.display();

  std::cout << "Column de *truc1" << std::endl;
  (*truc1).GetCol(1, vect1);
  vect1.display();

  std::cout << "Column de *truc2" << std::endl;
  (*truc2).GetCol(1, vect1);
  vect1.display();

  std::cout << "Column de *truc4" << std::endl;
  (*truc4).GetCol(1, vect1);
  vect1.display();

  std::cout << std::endl;
  std::cout << "********************* Affectation Row et Column d une Matrice *******************" << std::endl;
  std::cout << "vect3 " << std::endl;
  vect3.display();
  std::cout << std::endl;

  std::cout << "*truc1.SetRow(1,vect3) " << std::endl;
  (*truc1).SetRow(1, vect3);
  (*truc1).display();

  std::cout << "*truc4.SetRow(1,vect3) " << std::endl;
  (*truc4).SetRow(1, vect3);
  (*truc4).display();

  std::cout << std::endl;
  std::cout << "*truc2.SetCol(1,vect3) " << std::endl;
  (*truc2).SetCol(1, vect3);
  (*truc2).display();

  std::cout << std::endl;
  std::cout << "*truc5.SetCol(1,vect3) " << std::endl;
  (*truc5).SetCol(1, vect3);
  (*truc5).display();

  /*
  std::cout<<"vect1=Row(*truc1,1) "<<std::endl;
  vect1 = Row(*truc1,1);
  std::cout<<"vect1 "<<std::endl;
  vect1.display();
  */

  val1 = val2;
  std::cout << "val1, val1 = val2";
  val1.display();

  std::cout << std::endl;
  std::cout << "********************* Norme d une Matrice *******************" << std::endl;

  std::cout << "Norme de truc1: " << (*truc1).normInf() << std::endl;
  std::cout << std::endl;

  std::cout << "truc4: " << std::endl;
  (*truc4).display();
  std::cout << "Norme de truc4: " << (*truc4).normInf() << std::endl;
  std::cout << std::endl;

  (*truc2)(0, 1) = -3.L;
  std::cout << "Norme de truc2: " << (*truc2).normInf() << std::endl;
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "********************* GetBlock d une Matrice *******************" << std::endl;

  MySiconosMatrix *m = new MySimpleMatrix(DENSE, v6, 3, 2);
  DenseMat p(5, 5);
  for (unsigned i = 0; i < p.size1(); ++ i)
    for (unsigned j = 0; j < p.size2(); ++ j)
      p(i, j) = 5 * i + j;
  MySimpleMatrix val40(p);
  std::cout << "val40" << std::endl;
  val40.display();
  std::cout << std::endl;
  //std::cout<<"val40.Getblock(0,3, 0,2, *m)"<<std::endl;
  std::cout << "val40.Getblock(0, 0, *m)" << std::endl;
  //val40.GetBlock(0,3, 0,2, *m);
  val40.GetBlock(0, 0, *m);

  std::cout << "matrice m " << std::endl;
  m->display();

  std::cout << std::endl;
  std::cout << "********************* BlockMatrixCopy d une Matrice *******************" << std::endl;

  MySiconosMatrix *m10 = new MySimpleMatrix(DENSE, v6, 3, 2);
  std::cout << "matrice m10 " << std::endl;
  m10->display();

  std::cout << std::endl;
  std::cout << "val40.blockMatrixCopy(*m10, 0, 0)" << std::endl;
  val40.BlockMatrixCopy(*m10, 0, 0);

  std::cout << "val40 " << std::endl;
  val40.display();

  std::cout << "val1*2 " << std::endl;
  (val1 * 2.).display();

  if (val1 == val1)
  {
    std::cout << "ca marche" << std::endl;
  }
  else
    std::cout << "ca marche pas" << std::endl;

  /************************************* BLOCKMATRIX  ***************************************/

  std::cout << std::endl;
  std::cout << "************************************************************************" << std::endl;
  std::cout << "*********************          BlockMatrix           *******************" << std::endl;
  std::cout << "************************************************************************" << std::endl;
  std::cout << std::endl;

  std::vector<MySiconosMatrix* > vectSic1(4, truc1);
  std::vector<MySiconosMatrix* > vectSic2(4, truc2);
  std::vector<MySiconosMatrix* > vectSic3(4, truc3);

  MySiconosMatrix *block1 = new MyBlockMatrix(vectSic1, 2, 2);
  MySiconosMatrix *block1_2 = new MyBlockMatrix(vectSic1, 2, 2);
  MySiconosMatrix *block2 = new MyBlockMatrix(vectSic2, 2, 2);
  MySiconosMatrix *block3 = new MyBlockMatrix(vectSic3, 2, 2);
  MySiconosMatrix *block4 = new MyBlockMatrix(*block1);
  MyBlockMatrix block5(vectSic1, 2, 2);
  MyBlockMatrix block6(*block3);

  std::cout << std::endl;
  std::cout << "(*block1).blockMatrixCopy(*truc2, 0, 0)" << std::endl;
  (*block1).BlockMatrixCopy((*truc2), 0, 0);

  //---------------------------------------------------------------------------
  /*
  if( (*block1) == (*block1) )
    std::cout<<" SA MARCHEEEEEEEEEEEEE "<<std::endl;
  else
    std::cout<<" SA MARCHE PASSSSSSSSSSSSS"<<std::endl;
  */

  //---------------------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "********************* BlockMatrixCopy d une Matrice Block *******************" << std::endl;

  std::cout << std::endl;
  std::cout << "(*block2).blockMatrixCopy(*truc1, 0, 2)" << std::endl;
  (*block2).BlockMatrixCopy((*truc1), 0, 2);


  std::cout << std::endl;
  std::cout << "********************* GetBlock d une Matrice Block *******************" << std::endl;

  std::cout << std::endl;
  std::cout << "(*block2).Getblock(0, 0, val1)" << std::endl;
  (*block2).GetBlock(0, 0, val1);
  std::cout << "val1 " << std::endl;
  val1.display();

  std::cout << std::endl;
  std::cout << "(*block2).Getblock(0, 1, val1_2)" << std::endl;
  (*block2).GetBlock(0, 1, val1_2);
  std::cout << "val1_2 " << std::endl;
  val1_2.display();

  std::cout << std::endl;
  std::cout << "********************* acces aux dimensions d une Matrice Block *******************" << std::endl;

  std::cout << "(*block1).size1 () = " << (*block1).size1() << std::endl;
  std::cout << "(*block1).size2 () = " << (*block1).size2() << std::endl;

  std::cout << std::endl;
  std::cout << "********************* NormInf d une Matrice Block *******************" << std::endl;

  double d = (*block1_2).normInf();
  std::cout << "(*block1_2).normInf () = " << d << std::endl;

  std::cout << std::endl;
  std::cout << "********************* zero et eye d une Matrice Block *******************" << std::endl;


  // ZERO

  std::cout << "*block1 av zero " << std::endl;
  (*block1).display();
  std::cout << std::endl;
  (*block1).zero();
  std::cout << "*block1 apres zero " << std::endl;
  (*block1).display();
  std::cout << std::endl;




  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "*block1 av eye " << std::endl;
  (*block1).display();
  std::cout << std::endl;
  (*block1).eye();
  std::cout << "*block1 apres eye " << std::endl;
  (*block1).display();
  std::cout << std::endl;

  std::cout << std::endl;
  std::cout << "********************* GetRow et GetCol d une Matrice Block *******************" << std::endl;

  MySimpleVector vecteur1(DENSE);
  MySimpleVector vecteur2(DENSE);

  MySimpleVector vecteur3(DENSE, v3, 2);
  MySiconosVector *vecteur4 = new MySimpleVector(DENSE, v3, 2);



  std::cout << std::endl;
  std::cout << "********************* GetRow et GetCol d une Matrice Block *******************" << std::endl;

  MySimpleVector tmp1(DENSE, 4);
  MySimpleVector tmp2(DENSE, 4);

  std::cout << "block1 av *block1.GetCol " << std::endl;
  (*block1).display();

  (*block1).GetCol(1, tmp1);
  std::cout << "(*block1).GetCol(1, tmp1)" << std::endl;
  tmp1.display();


  (*block1).GetRow(0, tmp2);
  std::cout << "(*block1).GetRow(0, tmp2)" << std::endl;
  tmp2.display();


  (*block1).SetRow(0, tmp1);
  std::cout << "(*block1).SetRow(0, tmp1)" << std::endl;
  (*block1).display();

  (*block1).SetCol(1, tmp2);
  std::cout << "(*block1).SetCol(1, tmp2)" << std::endl;
  (*block1).display();


  std::cout << std::endl;
  std::cout << "********************* operateur ()(i, j) d une Matrice Block *******************" << std::endl;

  std::cout << "(*block1)(1, 1) = " << (*block1)(1, 1) << std::endl;



  std::cout << std::endl;
  std::cout << "********************* operateur = MyBlockMatrix, d une Matrice Block *******************" << std::endl;


  std::cout << "block6 :" << std::endl;
  block6.display();

  std::cout << "block1_2 av *block1_2 = block6" << std::endl;
  (*block1_2).display();

  std::cout << "block1_2 apres *block1 = block6" << std::endl;
  *block1_2 = block6;
  (*block1_2).display();



  std::cout << std::endl;
  std::cout << "******************** operateur = MySiconosMatrix, d une Matrice Block *******************" << std::endl;


  std::cout << "block2 :" << std::endl;
  (*block2).display();

  std::cout << "block1 av *block1 = *block2" << std::endl;
  (*block1).display();

  std::cout << "block1 apres *block1 = *block2" << std::endl;
  (*block1) = (*block2);
  (*block1).display();



  std::cout << std::endl;
  std::cout << "******************** operateur /= double, d une Matrice Block *******************" << std::endl;


  std::cout << "block1 apres *block1 /= 2" << std::endl;
  *block1 /= 2;
  (*block1).display();


  std::cout << std::endl;
  std::cout << "******************** operateur *= double, d une Matrice Block *******************" << std::endl;


  std::cout << "block1 apres *block1 *= 2" << std::endl;
  *block1 *= 2;
  (*block1).display();


  std::cout << std::endl;
  std::cout << "******************** operateur += MySiconosMatrix, d une Matrice Block *******************" << std::endl;



  std::cout << "block3 :" << std::endl;
  (*block3).display();

  std::cout << "block1_2 av *block1_2 += *block3" << std::endl;
  (*block1_2).display();

  std::cout << "block1_2 apres *block1_2 += *block3" << std::endl;
  (*block1_2) += (*block3);
  (*block1_2).display();


  std::cout << std::endl;
  std::cout << "******************** operateur -= MySiconosMatrix, d une Matrice Block *******************" << std::endl;


  std::cout << "block3 :" << std::endl;
  (*block3).display();


  std::cout << "block1_2 av *block1_2 -= *block3" << std::endl;
  (*block1_2).display();

  std::cout << "block1_2 apres *block1_2 -= *block3" << std::endl;
  (*block1_2) -= (*block3);
  (*block1_2).display();



  /************************************* VECTOR  ***************************************/

  std::cout << std::endl;
  std::cout << "************************************************************************" << std::endl;
  std::cout << "***********************          VECTOR           **********************" << std::endl;
  std::cout << "************************************************************************" << std::endl;
  std::cout << std::endl;
  /*
    const std::vector<double> v0 (4,0);
    const std::vector<double> v1 (4,1);
    const std::vector<double> v3 (2,3);
    const std::vector<double> v5 (25,0);
    const std::vector<double> v6 (6,0);

    MySimpleVector vect3 (DENSE,v3,2);
    MySimpleVector vect1 (DENSE);
    MySiconosVector *vect10 = new MySimpleVector (DENSE,v3,2);
  */


  //const std::vector<double> v0 (4,0);
  //const std::vector<double> v1 (4,1);
  /*
    const scalar_vector<double> sc1 (4,1);
    const scalar_vector<double> sc0 (4,0);
    mapped_vector<double> mVect1 (sc1);
    mapped_vector<double> mVect0 (sc0);

    MySiconosVector     *trucVect1   = new MySimpleVector(DENSE,v1,4);
    MySiconosVector     *trucVect1_2 = new MySimpleVector(DENSE,v1,4);
    MySiconosVector     *trucVect1_3 = new MySimpleVector(DENSE,v0,4);

    MySiconosVector     *trucVect2   = new MySimpleVector(mVect1);
    MySiconosVector     *trucVect2_2 = new MySimpleVector(mVect0);


    std::cout<< std::endl;
    std::cout << "****************** Multiplication entre Vecteurs *****************" << std::endl;

    std::cout<<"*trucVect1 "<<std::endl;
    (*trucVect1).display ();

    std::cout<<"*trucVect1_2 "<<std::endl;
    (*trucVect1_2).display ();

    std::cout<<"*trucVect1_3 "<<std::endl;
    (*trucVect1_3).display ();

    std::cout<<"*trucVect2 "<<std::endl;
    (*trucVect2).display ();

    std::cout<<"*trucVect2_2 "<<std::endl;
    (*trucVect2_2).display ();


    std::cout << "outer_prod(*trucVect1, *trucVect2) ";
    (outer_prod(*trucVect1, *trucVect2)).display();


    std::cout<<"inner_prod(*trucVect1, *trucVect2) "<< inner_prod(*trucVect1, *trucVect2) <<std::endl;

  */

  /************************************* ioMatrix  ***************************************/

  std::cout << std::endl;
  std::cout << "************************************************************************" << std::endl;
  std::cout << "***********************          ioMatrix           **********************" << std::endl;
  std::cout << "************************************************************************" << std::endl;
  std::cout << std::endl;


  ioMatrix io("fichierTestSimpleMat", "ascii");
  //ioMatrix io ("fichierTestSimpleMat", "binary");


  std::cout << std::endl;
  std::cout << "****************** ioMatrix::write *****************" << std::endl;


  std::cout << std::endl;
  std::cout << "		*********** write matrice simple ***********" << std::endl;

  std::cout << "*truc1.display () " << std::endl;
  (*truc1).display();

  io.write(*truc1);


  std::cout << std::endl;
  std::cout << "		*********** read matrice simple ***********" << std::endl;

  MySiconosMatrix *copy1   = new MySimpleMatrix(DENSE, 2, 2);
  io.read(*copy1);
  std::cout << "*copy1.display () " << std::endl;
  (*copy1).display();
  std::cout << "*truc1.display () " << std::endl;
  (*truc1).display();

  std::cout << "*copy1 * 2 " << std::endl;
  (*copy1) *= 2;
  std::cout << "*copy1.display () apres" << std::endl;
  (*copy1).display();
  std::cout << "*truc1.display () " << std::endl;
  (*truc1).display();



  ioMatrix io_block("fichierTestBlockMat", "ascii");
  //ioMatrix io_block ("fichierTestBlockMat", "binary");

  std::cout << std::endl;
  std::cout << "		*********** write matrice block ***********" << std::endl;

  std::cout << "*block1.display () " << std::endl;
  (*block1).display();

  //io.write(*block1);
  io_block.write(*block1);



  std::cout << std::endl;
  std::cout << "		*********** read matrice block ***********" << std::endl;

  MySiconosMatrix *truc1_4 = new MySimpleMatrix(DENSE, v0, 2, 2);
  std::vector<MySiconosMatrix* > vectSic1_bis(4, truc1_4);
  /*
  MyBlockMatrix copy2 (vectSic1_bis, 2, 2);
  io.read(copy2);
  std::cout<<"copy2.display () "<<std::endl;
  copy2.display ();
  */
  MySiconosMatrix *copy2 = new MyBlockMatrix(vectSic1_bis, 2, 2);
  std::cout << "*copy2.display () avant " << std::endl;
  (*copy2).display();

  //io.read(*copy2);
  io_block.read(*copy2);

  std::cout << "*copy2.display () apres " << std::endl;
  (*copy2).display();


  std::cout << "*truc1.display () " << std::endl;
  (*truc1).display();

  if ((*truc1) == (*truc1))
    std::cout << " SA MARCHEEEEEEEEEEEEE " << std::endl;
  else
    std::cout << " SA MARCHE PASSSSSSSSSSSSS" << std::endl;


  std::cout << "*truc1.display () " << std::endl;
  (*truc1).display();

  (*copy2) *= 2;
  std::cout << "*copy2.display () apres produit" << std::endl;
  (*copy2).display();

  std::cout << "*block1.display () " << std::endl;
  (*block1).display();

  /************** OPERATIONS AVEC VALEURS ******************/

  delete(block1);
  delete(block1_2);
  delete(block2);
  delete(block3);
  delete(block4);
  delete(vect10);

  delete(truc1);
  delete(truc1_2);
  delete(truc1_3);
  delete(truc2);
  delete(truc2_2);
  delete(truc3);
  delete(truc3_2);

  delete(vecteur4);


  delete(copy1);
  delete(copy2);
  delete(truc1_4);
  delete(m);
  delete(m10);


  /*
    delete(trucVect1);
    delete(trucVect1_2);
    delete(trucVect1_3);

    delete(trucVect2);
    delete(trucVect2_2);

    delete(&v0);
    delete(&v1);
    delete(&v3);
    delete(&v5);
    delete(&v6);
    delete(&vect3);
    delete(&vect1);


    delete(&val1);
    delete(&val1_2);
    delete(&val1_3);
    delete(&val2);
    delete(&val2_1);
    delete(&val2_2);
    delete(&val3);
    delete(&val3_1);
    delete(&val3_2);

    delete(&val40);

    delete(&p);

    delete(&vectSic1);
    delete(&vectSic2);
    delete(&vectSic3);
    delete(&block5);
    delete(&block6);

    delete(&vecteur1);
    delete(&vecteur2);
    delete(&vecteur3);


    delete(&tmp1);
    delete(&tmp2);

    delete(&sc1);
    delete(&sc0);
    delete(&mVect1);
    delete(&mVect0);
  */

};
