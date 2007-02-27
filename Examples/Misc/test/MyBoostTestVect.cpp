#include "ioVector.h"

int main()
{

  const std::vector<double> v0(4, 0);
  const std::vector<double> v1(4, 1);

  const scalar_vector<double> sc1(4, 1);
  const scalar_vector<double> sc0(4, 0);
  mapped_vector<double> m1(sc1);
  mapped_vector<double> m0(sc0);

  MySiconosVector     *truc1   = new MySimpleVector(DENSE, v1, 2);
  MySiconosVector     *truc1_2 = new MySimpleVector(DENSE, v1, 2);
  MySiconosVector     *truc1_3 = new MySimpleVector(DENSE, v0, 2);

  MySiconosVector     *truc2   = new MySimpleVector(m1);
  MySiconosVector     *truc2_2 = new MySimpleVector(m1);
  MySiconosVector     *truc3   = new MySimpleVector(m1);
  MySiconosVector     *truc3_2 = new MySimpleVector(m1);


  MySimpleVector val1(DENSE, v1, 2);
  MySimpleVector val1_2(DENSE, v1, 2);
  MySimpleVector val1_3(DENSE, v0, 2);

  MySimpleVector val2(m1);
  MySimpleVector val2_1(m0);
  MySimpleVector val2_2(m0);

  MySimpleVector val3(m1);
  MySimpleVector val3_1(m0);
  MySimpleVector val3_2(m0);


  std::cout << "*truc1 ";
  (*truc1).display();

  std::cout << "*truc2 ";
  (*truc2).display();

  std::cout << "*truc3 ";
  (*truc3).display();

  std::cout << "val1 ";
  val1.display();

  std::cout << "val1_2 ";
  val1_2.display();

  std::cout << "val3 ";
  val3.display();

  std::cout << "val3_1 ";
  val3_1.display();

  std::cout << "val3_2 ";
  val3_2.display();

  std::cout << std::endl;

  val1 -= val1;
  std::cout << "val1, val1-=val1";
  val1.display();

  /************** OPERATIONS AVEC POINTEURS ******************/

  std::cout << "************************** operation = ***************************" << std::endl;
  *truc1 = *truc2;
  std::cout << "*truc1 = *truc2 ";
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
  std::cout << std::endl;
  std::cout << "*************************  operation +=  *************************" << std::endl;
  *truc1 += *truc1;
  std::cout << "*truc1 += *truc1 ";
  (*truc1).display();

  *truc2 += *truc2;
  std::cout << "*truc2 += *truc2 ";
  (*truc2).display();

  *truc1 += *truc3;
  std::cout << "*truc1 += *truc3 ";
  (*truc1).display();
  std::cout << std::endl;
  std::cout << "*************** Restauration des Vecteurs initiaux **************" << std::endl;
  *truc1 = *truc1_2;
  std::cout << "*truc1 ";
  (*truc1).display();

  *truc2 = *truc2_2;
  std::cout << "*truc2 ";
  (*truc2).display();

  *truc3 = *truc3_2;
  std::cout << "*truc3 ";
  (*truc3).display();

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

  std::cout << std::endl;
  std::cout << "*********************** Somme de Vecteurs ************************" << std::endl;
  std::cout << "*truc2 + *truc3 " << std::endl;
  (add(*truc2, *truc3)).display();

  std::cout << "num de *truc2 + *truc3 " << (*truc2 + *truc3).GetNum() << std::endl;

  std::cout << "*truc2 + *truc2_2 " << std::endl;
  (*truc2 + *truc2_2).display();

  std::cout << "num de *truc2 + *truc2_2 " << (*truc2 + *truc2_2).GetNum() << std::endl;

  std::cout << "*truc1 + *truc3 " << std::endl;
  (add(*truc1, *truc3)).display();

  std::cout << "*truc2 + *truc1 " << std::endl;
  (add(*truc2, *truc1)).display();

  std::cout << std::endl;
  std::cout << "******************* Division avec un scalaire ********************" << std::endl;
  std::cout << "*truc1/2 " << std::endl;
  (*truc1 / 2.).display();

  std::cout << "*truc2/2. " << std::endl;
  (*truc2 / 2.).display();

  std::cout << "*truc3/3 " << std::endl;
  (*truc3 / 3.).display();

  std::cout << std::endl;
  std::cout << "****************** Multiplication avec scalaire ******************" << std::endl;
  std::cout << "*truc1*2 " << std::endl;
  (*truc1 * 2.).display();

  std::cout << "*truc2*2. " << std::endl;
  (*truc2 * 2.).display();

  std::cout << "*truc3*3 " << std::endl;
  (*truc3 * 3.).display();

  std::cout << std::endl;
  std::cout << "****************** Multiplication entre Vecteurs *****************" << std::endl;
  std::cout << "outer_prod(*truc2, *truc2) ";
  (outer_prod(*truc2, *truc2)).display();

  std::cout << "inner_prod(*truc2, *truc2) " << inner_prod(*truc2, *truc2) << std::endl;

  std::cout << "outer_prod(*truc1, *truc2) ";
  (outer_prod(*truc1, *truc2)).display();

  std::cout << "inner_prod(*truc1, *truc2) " << inner_prod(*truc1, *truc2) << std::endl;

  std::cout << "outer_prod(*truc2, *truc3) " << std::endl;
  (outer_prod(*truc2, *truc3)).display();

  std::cout << "inner_prod(*truc2, *truc3) " << inner_prod(*truc2, *truc3) << std::endl;

  std::cout << std::endl;
  std::cout << "***************** Acces aux elements d un Vecteur ****************" << std::endl;

  std::cout << "*truc1(0) " << (*truc1)(0) << std::endl;

  std::cout << "*truc2(1) " << (*truc2)(1) << std::endl;

  std::cout << "*truc3(0) " << (*truc3)(0) << std::endl;

  std::cout << "val1(0) " << val1(0) << std::endl;

  std::cout << "val2(1) " << val2(1) << std::endl;

  std::cout << "val3(0) " << val3(0) << std::endl;

  /************************************* ioVector  ***************************************/

  std::cout << std::endl;
  std::cout << "************************************************************************" << std::endl;
  std::cout << "***********************          ioVector           **********************" << std::endl;
  std::cout << "************************************************************************" << std::endl;
  std::cout << std::endl;


  ioVector io("fichierTestSimpleVect", "ascii");
  //ioVector io ("fichierTestSimpleVect", "binary");


  std::cout << std::endl;
  std::cout << "****************** ioVector::write *****************" << std::endl;


  std::cout << std::endl;
  std::cout << "		*********** write simple vector ***********" << std::endl;

  std::cout << "*truc1.display () " << std::endl;
  (*truc1).display();

  io.write(*truc1);


  std::cout << std::endl;

  std::cout << "		*********** read simple vector ***********" << std::endl;

  MySiconosVector *copy1   = new MySimpleVector(DENSE, 2);
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

  std::cout << std::endl;

};
