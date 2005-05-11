class X();

X f()
{
  return X();  // return by value
}

void g1(X&) {} // Pass by non-const reference

void g2(const X&) {} // Pass by const reference

int main()
{
  // g1(f()); Error: const temporary created by f()
  g2(f()); // ok
}
