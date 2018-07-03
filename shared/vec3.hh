#ifndef VEC3_HH
#define VEC3_HH

#include <cmath>

struct vec3 {
	public:
		double x,y,z;
		vec3() {};
		vec3(double x_) : x(x_),y(x_),z(x_) {};
		vec3(double x_,double y_,double z_) : x(x_),y(y_),z(z_) {};
		inline vec3 operator+ (vec3 p) {
			return vec3(x+p.x,y+p.y,z+p.z);
		}
		inline vec3 operator- (vec3 p) {
			return vec3(x-p.x,y-p.y,z-p.z);
		}
		inline vec3 operator* (vec3 b) {
			return vec3(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);
		}
		inline vec3 operator* (double a) {
			return vec3(x*a,y*a,z*a);
		}
		inline vec3 operator/ (double a) {
			return vec3(x/a,y/a,z/a);
		}
		inline void operator+= (vec3 p) {
			x+=p.x;y+=p.y;z+=p.z;
		}
		inline void operator-= (vec3 p) {
			x-=p.x;y-=p.y;z-=p.z;
		}
		inline void operator*= (double a) {
			x*=a;y*=a;z*=a;
		}
		inline void operator/= (double a) {
			x/=a;y/=a;z/=a;
		}
		inline operator double() {
			return double(x*x+y*y+z*z);
		}
		inline double magnitude() {
			return sqrt(x*x+y*y+z*z);
		}
		inline void normalize() {
			double mag=x*x+y*y+z*z;
			if(mag>0) {
				mag=1./sqrt(mag);
				x*=mag;y*=mag;z*=mag;
			}
		}
		inline void add(double *p) {
			*p+=x;p[1]+=y;p[2]+=z;
		}
};

inline vec3 operator*(const double e,vec3 f) {
	return vec3(e*f.x,e*f.y,e*f.z);
}

inline vec3 operator-(vec3 f) {
	return vec3(-f.x,-f.y,-f.z);
}

inline double dot(vec3 a,vec3 b) {
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

inline double mod_sq(vec3 a) {return a.x*a.x+a.y*a.y+a.z*a.z;}

#endif
