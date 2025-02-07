#include <doro/probability.hpp>
#include <iostream>

using namespace std;
using namespace Doro;

int main() {
  double alpha = 2.0, epsilon = 1e-15, b = 10000;
  // for (double beta = 2.0; beta < 50; beta += 0.5) {
  //   cout << loss_func(alpha,1.0 / beta, epsilon) << endl;
  // }
  auto lam = [alpha, epsilon, b](double beta) {return -loss_func(alpha, 1.0 / beta, b, epsilon);};
  cout << convex_argmax(lam, 2.0, 200.0, epsilon = 0.01) << endl;
  return 0;
}