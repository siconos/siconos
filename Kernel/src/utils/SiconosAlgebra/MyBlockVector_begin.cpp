#include <vector>
#include <iostream>
#include <stdlib.h>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;
typedef vector<double, std::vector<double> > MySimpleVector;

extern "C" void mass(const unsigned int sizeOfq, const double *q, double *mass)
{
  for (unsigned i = 0; i < sizeOfq; i++)
    mass[i] = 2 * q[i];
};


int main()
{
  using namespace boost::numeric::ublas;

  typedef vector<double, std::vector<double> > MySimpleVector;
  typedef std::vector<MySimpleVector* > MyBlockVector;

  const std::vector<double> v0(4, 0);
  const std::vector<double> v1(4, 1);
  const std::vector<double> v2(8, 0);

  MySimpleVector *vect1 = new MySimpleVector(4, v1);
  MySimpleVector *vect0 = new MySimpleVector(4, v0);

  double *a = (double*)malloc(sizeof(double));
  double *b = (double*)malloc(sizeof(double));
  *a = 1;
  *b = 12;
  std::cout << "*a avant" << *a << std::endl;
  std::cout << "*b avant" << *b << std::endl;
  mass(1, a, b);
  std::cout << "*a apres" << *a << std::endl;
  std::cout << "*b apres" << *b << std::endl;
  MySimpleVector *res = new MySimpleVector(4, v0);
  std::cout << "res avant " << *res << std::endl;
  mass(4, &(*vect1)(0) , &(*res)(0));
  std::cout << "res apres " << *res << std::endl;
  int n1 = sizeof(MySimpleVector(4, v0));
  //size_t n1 = sizeof(*vect1);
  int n0 = sizeof(*vect0);
  //size_t Taille = n0 + n1;

  MySimpleVector *tmp0 = (MySimpleVector*)malloc(sizeof(*vect0) + sizeof(*vect1));
  //MySimpleVector *tmp0 = (MySimpleVector*)malloc(Taille);

  MySimpleVector *a_liberer;
  a_liberer = vect0;  // ON FIXE UN POINTEUR SUR LA MEMOIRE A LIBERER
  vect0 = tmp0;       // ON TRANSLATE VECT0 SUR LA NOUVELLE MEMOIRE ALLOUEE
  std::cout << "a_liberer " << *a_liberer << std::endl;
  std::cout << "vect0 " << *vect0 << std::endl;
  //*vect0 = *a_liberer;      //--------COPIE DU PREMIER VECTEUR
  //memcpy(&(*vect0)(0), &(*a_liberer)(0), n0);
  memcpy(vect0, a_liberer, n0);
  std::cout << "vect0 " << *vect0 << std::endl;
  std::cout << "&vect0 " << vect0 << std::endl;
  std::cout << "&a_liberer " << a_liberer << std::endl;


  //free(a_liberer);
  a_liberer = vect1;
  std::cout << "resultat " << std::endl;
  std::cout << "&vect1 " << vect1 << std::endl;
  std::cout << "&a_liberer " << a_liberer << std::endl;
  std::cout << "a_liberer " << *a_liberer << std::endl;
  MySimpleVector *tmp1 = tmp0;
  //tmp1 = tmp0 + n1;
  tmp1 = tmp0 - 1;
  vect1 = tmp1;
  std::cout << "resultat " << std::endl;

  //*vect1 = *a_liberer;
  //memcpy(&(*vect1)(0), &(*a_liberer)(0), n1);
  memcpy(vect1, a_liberer, n1);
  std::cout << "resultat " << std::endl;
  std::cout << "&vect1 " << vect1 << std::endl;
  std::cout << "&a_liberer " << a_liberer << std::endl;
  std::cout << "a_liberer " << *a_liberer << std::endl;
  //delete(a_liberer);

  MyBlockVector block ;
  std::cout << "block " << &block << std::endl;
  std::cout << "resultat " << std::endl;
  block.reserve(2);

  std::cout << "resultat " << std::endl;
  block.push_back(vect0);
  std::cout << "resultat " << std::endl;
  block.push_back(vect1);
  std::cout << "resultat " << std::endl;
  std::cout << "block0 " << *block[0] << std::endl;
  std::cout << "block1 " << *block[1] << std::endl;


  std::cout << "tmp0 " << *tmp0 << std::endl;
  std::cout << "tm1 " << *tmp1 << std::endl;

  std::cout << "&tmp0 " << tmp0 << std::endl;
  std::cout << "&tm1 " << tmp1 << std::endl;
  std::cout << "&vect0 " << vect0 << std::endl;
  std::cout << "&vect1 " << vect1 << std::endl;
  std::cout << "n1 " << n1 << std::endl;
  std::cout << "n0 " << n0 << std::endl;
  std::cout << "tmp1 - tmp0 " << (tmp1 - tmp0) << std::endl;

  MySimpleVector *res_block = new MySimpleVector(8, v2);
  std::cout << "size double " << sizeof(double) << std::endl;
  std::cout << "size int " << sizeof(unsigned int) << std::endl;
  std::cout << "size std::vector 4 " << sizeof(std::vector<double>(4)) << std::endl;
  std::cout << "&(*vect1)(0) - &(*vect0)(3) = " << &(*vect1)(0) - &(*vect0)(3) << std::endl;

  double *cour = &(*vect1)(0);
  double *cour0 = &(*vect0)(0);
  std::cout << "cour[0] = " << cour[0] << std::endl;
  std::cout << "cour[1] = " << cour[1] << std::endl;
  std::cout << "cour[2] = " << cour[2] << std::endl;
  std::cout << "cour[3] = " << cour[3] << std::endl;
  std::cout << "cour[4] = " << cour[4] << std::endl;
  std::cout << "cour[5] = " << cour[5] << std::endl;
  std::cout << "cour0[0] = " << cour0[0] << std::endl;
  std::cout << "cour0[1] = " << cour0[1] << std::endl;
  std::cout << "cour0[2] = " << cour0[2] << std::endl;
  std::cout << "cour0[3] = " << cour0[3] << std::endl;
  std::cout << "cour0[4] = " << cour0[4] << std::endl;
  std::cout << "cour0 = " << cour0 << std::endl;
  std::cout << "cour = " << cour << std::endl;
  std::cout << "cour0-cour = " << cour0 - cour << std::endl;

  std::cout << "res_block avant " << *res_block << std::endl;
  //mass(8, &(*tmp0)(0) , &(*res_block)(0) );
  //mass1(2, vect0 , res_block );
  mass(8, &(*vect0)(0) , &(*res_block)(0));
  //mass(8, &(*block[0])(0) , &(*res_block)(0) );
  std::cout << "res_block apres " << *res_block << std::endl;


  delete(res);
  delete(res_block);
  free(tmp0);
  delete(cour);
  delete(cour0);
  return 0;
}
