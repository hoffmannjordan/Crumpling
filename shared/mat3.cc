#include "mat3.hh"
using namespace std;

void sym_mat3::eigenvalues(double &eig1,double &eig2,double &eig3) {
	double m=trace()*(1/3.0);

	// Temporarily modify matrix elements to make the computation of the
	// determinant and deviatoric modulus more straightforward
	double b11=a11-m,b22=a22-m,b33=a33-m;
	double q=0.5*(b11*(b22*b33-a23*a23)+a12*(2*a13*a23-a12*b33)-a13*a13*b22),
	       p=(1/6.0)*(b11*b11+b22*b22+b33*b33)+(1/3.0)*(a12*a12+a13*a13+a23*a23);

	// Compute intermediate values
	double sp=sqrt(p),fac=p*sp,phi;
	phi=abs(q)<abs(fac)?(1/3.0)*acos(q/fac):
			    (q<0?3.1415926535897932384626433832795/3.0:0);

	// Compute eigenvalues
	double cphi=cos(phi),sphi=sin(phi);
	eig1=m+2*sp*cphi;
	eig2=m-sp*(cphi-sqrt(3)*sphi);
	eig3=m-sp*(cphi+sqrt(3)*sphi);
}

void sym_mat3::eigenvectors(double &eig1,double &eig2,double &eig3,mat3 &Lam) {

	// Call routine to compute eigenvalues
	eigenvalues(eig1,eig2,eig3);

	// Calculate threshold to decide when eigenvalues should be treated as
	// equal
	double t=1e-8*(abs(a11)+abs(a22)+abs(a33)+2*(abs(a12)+abs(a13)+abs(a23))),nor;

	// Swap eigenvalues around to deal with the one of largest magnitude
	if(abs(eig3)>abs(eig1)) {nor=eig3;eig3=eig1;eig1=nor;}

	// Test to see if the first two eigenvalues are equal
	if(nonzero(eig1-eig2,t)) {

		// If the eigenvalues are distinct then just call the regular
		// routine to find independent eigenvalues
		compute_eigenvector(eig1,Lam.a11,Lam.a21,Lam.a31,t);
		nor=eig3;eig3=eig2;eig2=nor;
		compute_eigenvector(eig2,Lam.a12,Lam.a22,Lam.a32,t);
	} else {

		// If the eigenvalues are equal, then the matrix has rank of at
		// most one; since it is symmetric this gives it a very
		// specific form. Use this to calculate two orthogonal
		// eigenvectors
		double b11=a11-eig1;
		if(nonzero(b11,t)) {
			if(nonzero(a12,t)) {
				nor=1/sqrt(a12*a12+a13*a13);
				Lam.a11=0;Lam.a21=a13*nor;Lam.a31=-a12*nor;
				nor=1/sqrt(1+b11*b11);
				Lam.a12=nor;Lam.a22=Lam.a21*nor;Lam.a32=-Lam.a31*nor;
			} else {
				nor=1/sqrt(b11*b11+a13*a13);
				Lam.a11=a13*nor;Lam.a22=0;Lam.a32=-b11*nor;
				Lam.a12=0;Lam.a22=1;Lam.a32=0;
			}
		} else {
			Lam.a11=1;Lam.a21=0;Lam.a31=0;
			double b22=a22-eig1;
			if(nonzero(b22,t)) {
				nor=1/sqrt(b22*b22+a23*a23);
				Lam.a12=0;Lam.a22=a23*nor;Lam.a32=-b22*nor;
			} else {
				Lam.a12=0;Lam.a22=1;Lam.a32=0;
			}
		}
	}

	// Compute third eigenvector as the vector product of the first two
	Lam.a13=Lam.a21*Lam.a32-Lam.a31*Lam.a22;
	Lam.a23=Lam.a31*Lam.a12-Lam.a11*Lam.a32;
	Lam.a33=Lam.a11*Lam.a22-Lam.a21*Lam.a12;
}

void sym_mat3::compute_eigenvector(double eig,double &a,double &b,double &c,double t) {
	double b11=a11-eig,nor;

	// Split into three cases depending on whether entries in the first
	// column are zero or not
	if(nonzero(b11,t)) {
		double q=(a22-eig)-a12*(a12/b11),r=a23-a13*(a12/b11);
		if(nonzero(q,t)) {
			a=(-a12*r+a13*q);
			b=b11*r;
			c=-b11*q;
		} else {
			nor=1/sqrt(b11*b11+a12*a12);
			a=a12*nor;
			b=-b11*nor;
			c=0;
			return;
		}
	} else if(nonzero(a12*a12,t)) {
		a=a12*a23-a13*(a22-eig);
		b=a12*a13;
		c=-a12*a12;
	} else if(nonzero(a13*a13*a13,t)) {
		nor=1/sqrt(a13*a13+a23*a23);
		a=a23*nor;
		b=-a13*nor;
		c=0;
		return;
	} else {
		a=1;b=c=0;return;
	}

	// Any case with three non-zero entries needs full normalization
	nor=1/sqrt(a*a+b*b+c*c);
	a*=nor;b*=nor;c*=nor;
}

void sym_mat3::project_to_other(double v1,double v2,double v3,sym_mat3 &s) {
	double c1((a11-s.a11)*v1+(a12-s.a12)*v2+(a13-s.a13)*v3),
	       c2((a12-s.a12)*v1+(a22-s.a22)*v2+(a23-s.a23)*v3),
	       c3((a13-s.a13)*v1+(a23-s.a23)*v2+(a33-s.a33)*v3),
	       t(c1*v1+c2*v2+c3*v3);
	a11-=2*c1*v1-t*v1*v1;
	a12-=c1*v2+c2*v1-t*v1*v2;
	a13-=c1*v3+c3*v1-t*v1*v3;
	a22-=2*c2*v2-t*v2*v2;
	a23-=c2*v3+c3*v2-t*v2*v3;
	a33-=2*c3*v3-t*v3*v3;
}

void sym_mat3::project(double v1,double v2,double v3) {
	double c1(a11*v1+a12*v2+a13*v3),
	       c2(a12*v1+a22*v2+a23*v3),
	       c3(a13*v1+a23*v2+a33*v3),t(c1*v1+c2*v2+c3*v3);
	a11-=2*c1*v1-t*v1*v1;
	a12-=c1*v2+c2*v1-t*v1*v2;
	a13-=c1*v3+c3*v1-t*v1*v3;
	a22-=2*c2*v2-t*v2*v2;
	a23-=c2*v3+c3*v2-t*v2*v3;
	a33-=2*c3*v3-t*v3*v3;
}

sym_mat3 sym_mat3::trans_mult() {
	return sym_mat3(a11*a11+a12*a12+a13*a13,a11*a12+a12*a22+a13*a23,a11*a13+a12*a23+a13*a33,
			a12*a12+a22*a22+a23*a23,a12*a13+a22*a23+a23*a33,a13*a13+a23*a23+a33*a33);
}

sym_mat3 mat3::trans_mult() {
	return sym_mat3(a11*a11+a12*a12+a13*a13,a11*a21+a12*a22+a13*a23,a11*a31+a12*a32+a13*a33,
			a21*a21+a22*a22+a23*a23,a21*a31+a22*a32+a23*a33,a31*a31+a32*a32+a33*a33);
}
