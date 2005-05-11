int main()
{

  // Derived1 a class derived from Base1, a class with no virtual destructor.
  Base1 *bp = new Derived1;
  delete bp; // call only the destructor of Base1 => potential bug

  // Derived2 a class derived from Base2, a class with a virtual destructor.
  Base2 *bp2 = new Derived2;
  delete bp2; // call the destructor of Derived2 and then the one of Base2
}
