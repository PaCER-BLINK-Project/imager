#include <limits>
#include <stdio.h>

int main()
{

    printf("Limits of double = %e - %e",-std::numeric_limits<double>::max(),std::numeric_limits<double>::max());
}