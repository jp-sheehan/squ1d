#include <iostream>
#include <omp.h>

int main(int argc, char *argv[])
{
  int i,j;

  int nthreads, tid;

/*  #pragma omp parallel private(nthreads,tid)
  {
    tid = omp_get_thread_num();
    std::cout << "Spaghetti Thread:  " << tid << std::endl;

    if (tid==0)
    {
      std::cout << "I am the Spaghetti King" << std::endl;
    }
  }*/

  return(0);
}
