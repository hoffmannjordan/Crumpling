#ifndef C3D_MAT3_HH
#define C3D_MAT3_HH

#include <cstdio>
#include <cmath>

class sym_mat3;

struct mat3 {
	double a11,a12,a13;
	double a21,a22,a23;
	double a31,a32,a33;
	mat3() {};
	mat3(double a) : a11(a), a12(0), a13(0), a21(0), a22(a), a23(0), a31(0), a32(0), a33(a) {};
	mat3(double a,double b,double c) : a11(a), a12(0), a13(0), a21(0), a22(b), a23(0), a31(0), a32(0), a33(c) {};
	mat3(double b11,double b12,double b13,double b21,double b22,double b23,double b31,double b32,double b33)
		: a11(b11), a12(b12), a13(b13), a21(b21), a22(b22), a23(b23), a31(b31), a32(b32), a33(b33) {};
	inline void set(double b11,double b12,double b13,double b21,double b22,double b23,double b31,double b32,double b33) {
		a11=b11;a12=b12;a13=b13;
		a21=b21;a22=b22;a23=b23;
		a31=b31;a32=b32;a33=b33;
	}
	inline mat3 operator+ (mat3 p) {
		return mat3(a11+p.a11,a12+p.a12,a13+p.a13,
			    a21+p.a21,a22+p.a22,a23+p.a23,
			    a31+p.a31,a32+p.a32,a33+p.a33);
	}
	inline mat3 operator- (mat3 p) {
		return mat3(a11-p.a11,a12-p.a12,a13-p.a13,
			    a21-p.a21,a22-p.a22,a23-p.a23,
			    a31-p.a31,a32-p.a32,a33-p.a33);
	}
	inline void operator*=(double p) {
		a11*=p;a12*=p;a13*=p;
		a21*=p;a22*=p;a23*=p;
		a31*=p;a32*=p;a33*=p;
	}
	inline mat3 operator* (double p) {
		return mat3(a11*p,a12*p,a13*p,
			    a21*p,a22*p,a23*p,
			    a31*p,a32*p,a33*p);
	}
	inline mat3 operator/ (double p) {
		double pinv(1/p);
		return mat3(a11*pinv,a12*pinv,a13*pinv,
			    a21*pinv,a22*pinv,a23*pinv,
			    a31*pinv,a32*pinv,a33*pinv);
	}
	inline mat3 operator* (mat3 p) {
		return mat3(a11*p.a11+a12*p.a21+a13*p.a31,a11*p.a12+a12*p.a22+a13*p.a32,a11*p.a13+a12*p.a23+a13*p.a33,
			    a21*p.a11+a22*p.a21+a23*p.a31,a21*p.a12+a22*p.a22+a23*p.a32,a21*p.a13+a22*p.a23+a23*p.a33,
			    a31*p.a11+a32*p.a21+a33*p.a31,a31*p.a12+a32*p.a22+a33*p.a32,a31*p.a13+a32*p.a23+a33*p.a33);
	}
	inline mat3 operator* (sym_mat3 p);
	inline double det() {
		return a11*(a22*a33-a23*a32)+a21*(a13*a32-a12*a33)+a31*(a12*a23-a13*a22);
	}
	inline mat3 transpose() {
		return mat3(a11,a21,a31,a12,a22,a32,a13,a23,a33);
	}
	inline mat3 inverse(double &dt) {
		dt=det();
		double idet(1/dt);
		return mat3((a22*a33-a23*a32)*idet,(a13*a32-a12*a33)*idet,(a12*a23-a13*a22)*idet,
			    (a23*a31-a21*a33)*idet,(a11*a33-a13*a31)*idet,(a13*a21-a11*a23)*idet,
			    (a21*a32-a22*a31)*idet,(a12*a31-a11*a32)*idet,(a11*a22-a12*a21)*idet);
	}
	inline double trace() {
		return a11+a22+a33;
	}
	inline void print() {
		printf("%g %g %g\n%g %g %g\n%g %g %g\n",a11,a12,a13,a21,a22,a23,a31,a32,a33);
	}
	inline double modsq() {
		return 0.5*(a11*a11+a12*a12+a13*a13+a21*a21+a22*a22+a23*a23+a31*a31+a32*a32+a33*a33);
	}
	inline double mod() {
		return sqrt(modsq());
	}
	sym_mat3 trans_mult();
};

class sym_mat3 {
	public:
		double a11,a12,a13;
		double a22,a23;
		double a33;
		sym_mat3() {};
		sym_mat3(double a_): a11(a_), a12(0) ,a13(0), a22(a_), a23(0), a33(a_) {};
		sym_mat3(double b11,double b12,double b13,double b22,double b23,double b33)
			: a11(b11), a12(b12), a13(b13), a22(b22), a23(b23), a33(b33) {};
		sym_mat3(double *b)
			: a11(*b), a12(b[1]), a13(b[2]), a22(b[3]), a23(b[4]), a33(b[5]) {};
		inline void put(double *p) {
			*(p++)=a11;*(p++)=a12;*(p++)=a13;*(p++)=a22;*(p++)=a23;*p=a33;
		}
		inline sym_mat3 operator+ (sym_mat3 p) {
			return sym_mat3(a11+p.a11,a12+p.a12,a13+p.a13,a22+p.a22,a23+p.a23,a33+p.a33);
		}
		inline sym_mat3 operator- (sym_mat3 p) {
			return sym_mat3(a11-p.a11,a12-p.a12,a13-p.a13,a22-p.a22,a23-p.a23,a33-p.a33);
		}
		inline sym_mat3 operator* (double p) {
			return sym_mat3(a11*p,a12*p,a13*p,a22*p,a23*p,a33*p);
		}
		inline sym_mat3 operator/ (double p) {
			double pinv(1/p);
			return sym_mat3(a11*pinv,a12*pinv,a13*pinv,a22*pinv,a23*pinv,a33*pinv);
		}
		inline double trace() {
			return a11+a22+a33;
		}
		inline double det() {
			return a11*(a22*a33-a23*a23)+a12*(2*a13*a23-a12*a33)-a13*a13*a22;
		}
		sym_mat3 trans_mult();

		inline void pack(double *&p) {
			*(p++)=a11;*(p++)=a12;*(p++)=a13;
			*(p++)=a22;*(p++)=a23;*(p++)=a33;
		}
		inline void print() {
			printf("%g %g %g\n* %g %g\n* * %g\n",a11,a12,a13,a22,a23,a33);
		}
		inline void clear() {
			a11=a12=a13=a22=a23=a33=0;
		}
		void project_to_other(double v1,double v2,double v3,sym_mat3 &s);
		void project(double v1,double v2,double v3);
		void eigenvalues(double &eig1,double &eig2,double &eig3);
		void eigenvectors(double &eig1,double &eig2,double &eig3,mat3 &Lam);
	private:
		void compute_eigenvector(double eig,double &a,double &b,double &c,double t);
		inline bool nonzero(double a,double threshold) {
			return std::abs(a)>threshold;
		}
};

inline mat3 operator*(const double p,mat3 m) {
	return mat3(p*m.a11,p*m.a12,p*m.a13,
		    p*m.a21,p*m.a22,p*m.a23,
		    p*m.a31,p*m.a32,p*m.a33);
}

inline sym_mat3 operator*(const double p,sym_mat3 m) {
	return sym_mat3(p*m.a11,p*m.a12,p*m.a13,
			p*m.a22,p*m.a23,p*m.a33);
}

inline mat3 mat3::operator* (sym_mat3 p) {
	return mat3(a11*p.a11+a12*p.a12+a13*p.a13,a11*p.a12+a12*p.a22+a13*p.a23,a11*p.a13+a12*p.a23+a13*p.a33,
		    a21*p.a11+a22*p.a12+a23*p.a13,a21*p.a12+a22*p.a22+a23*p.a23,a21*p.a13+a22*p.a23+a23*p.a33,
		    a31*p.a11+a32*p.a12+a33*p.a13,a31*p.a12+a32*p.a22+a33*p.a23,a31*p.a13+a32*p.a23+a33*p.a33);
}

#endif
