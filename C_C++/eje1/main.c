#include <stdlib.h>
#include <stdio.h>


int main(){
  int a = 0;
  long b = 0;
  float c = 0.;
  double d = 0.;
  
  d = 3.*6.+5.;
  printf("El resultado es: %f\n\r", d);
  printf("El tamaño de int es: %ld\n\r", sizeof(a));
  printf("El tamaño de long es: %ld\n\r", sizeof(b));
  printf("El tamaño de float es: %ld\n\r", sizeof(c));
  printf("El tamaño de double es: %ld\n\r", sizeof(d));

  return 0;
}
