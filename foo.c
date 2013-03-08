#include <math.h>

double f(int n)
{
    double v = 0.0;
    for (int i = 0; i < n; i++)
	v = v + log(i);
    return v;
}	
	
double g(int n)
{
    return f(n);
}

int main()
{
    int n = 10000;
    int N = 10000;
    for (int i = 0; i < N; i++)
	g(n);
}

    
