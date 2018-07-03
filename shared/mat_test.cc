#include <cstdio>
#include <cstdlib>

#include "mat3.hh"

double rnd() {
	if(rand()%10==0) return 0;
	return 2*(double(rand())/RAND_MAX)-1;
}

int main() {
	double e1,e2,e3;
	mat3 Lam,q;
	sym_mat3 m(1,0,1,1,0,1);
	m.eigenvectors(e1,e2,e3,Lam);
	printf("%g %g %g\n",e1,e2,e3);
	Lam.print();puts("");
	q=Lam*mat3(e1,e2,e3)*Lam.transpose();
	q.print();
}
