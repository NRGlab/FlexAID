#include "flexaid.h"

double GetValueFromGaussian(double x,double max,double zero) {
	return pow( -(x-zero) * (x-(2*max-zero)) / (pow(zero-max,2.0)), 50.0 );
}

#ifdef _WIN32
 void vec_sub(float *a, const float *b, const float *c) {
#else
 void vec_sub(float *a, const float *b, const float *c) {
#endif
   a[0]=b[0]-c[0];
   a[1]=b[1]-c[1];
   a[2]=b[2]-c[2];
 }

#ifdef _WIN32
 float dot_prod(const float *v1, const float *v2) {
#else
 float dot_prod(const float *v1, const float *v2) {
#endif
   return v1[0]* v2[0] + v1[1]* v2[1] + v1[2] * v2[2];
 }

#ifdef _WIN32
 float distance2(const float *a, const float *b) {
#else
 float distance2(const float *a, const float *b){
#endif
  float delta = a[0] - b[0];
  float r2 = delta*delta;
  delta = a[1] - b[1];
  r2 += delta*delta;
  delta = a[2] - b[2];
  return r2 + delta*delta;
}

float * cross_prod(float *x1, const float *x2, const float *x3) {
  x1[0] =  x2[1]*x3[2] - x3[1]*x2[2];
  x1[1] = -x2[0]*x3[2] + x3[0]*x2[2];
  x1[2] =  x2[0]*x3[1] - x3[0]*x2[1];
  return x1;
}

float distance(const float *a, const float *b) {
  return sqrtf(distance2(a,b));
}

float angle(const float *a, const float *b, const float *c) {
  float v1[3], v2[3];
  float ab[3];
  float psin, pcos;
  vec_sub(v1,a,b);
  vec_sub(v2,b,c);
  cross_prod(ab, v1, v2);
  psin = sqrtf(dot_prod(ab, ab));
  pcos = dot_prod(v1, v2);
  return 57.2958f * atan2(psin, pcos);
}

float dihedral(const float *a1,const float *a2,const float *a3,const float *a4) {
  float r1[3], r2[3], r3[3], n1[3], n2[3];
  float psin, pcos;
  vec_sub(r1, a2, a1);
  vec_sub(r2, a3, a2);
  vec_sub(r3, a4, a3);
  
  cross_prod(n1, r1, r2);
  cross_prod(n2, r2, r3);
  
  psin = dot_prod(n1, r3) * sqrtf(dot_prod(r2, r2));
  pcos = dot_prod(n1, n2);

  return 57.2958f * atan2(psin, pcos);
}

/******************************************************************************
 * SUBROUTINE zeros calculates the zeros of a quadratic function
 * the function returns the positive zero
 ******************************************************************************/
float zero(float a, float b, float c){
	
	return (
		( -b + sqrtf( b*b - 4*a*c ) / ( 2 * a ) ) > 0 ?
		( -b + sqrtf( b*b - 4*a*c ) / ( 2 * a ) ) :
		( -b - sqrtf( b*b - 4*a*c ) / ( 2 * a ) )
		);

}

/******************************************************************************
 * SUBROUTINE distance calculates the cartesian distance between two point in 
 * n dimensions.
 ******************************************************************************/
float distance_n(float a[], float b[], int n){
  float d;
  int i;
  
  d=0.0;
  for(i=0;i<=n-1;i++){
    d +=(a[i]-b[i])*(a[i]-b[i]); 
  }
  d = sqrt(d);

  return (d);
}

/******************************************************************************
 * SUBROUTINE sqrdist calculates the square distance between two points in d=3.
 ******************************************************************************/
float sqrdist(float a[], float b[]){
  float d;
  int i;
  
  d=0.0;
  for(i=0;i<=2;i++){
    d += (a[i]-b[i])*(a[i]-b[i]); 
  }
  return (d);
}

/******************************************************************************
 * SUBROUTINE distance calculates the cartesian distance between two point in 
 * n dimensions.
 ******************************************************************************/
float dist(float a[], float b[]){
  float d;
  int i;
  
  d=0.0;
  for(i=0;i<=2;i++){
    d += (a[i]-b[i])*(a[i]-b[i]); 
  }
  d = sqrt(d);

  return (d);
}

/*******************************************************************************
 * SUBROUTINE bndang calculates the valence angle between three cosecutive atoms
 *******************************************************************************/
float bndang(float a[],float b[], float c[]){
  float absu,absv,cosa;
  int i;

  cosa=0.0;
  absu=0.0;
  absv=0.0;
  for(i=0;i<=2;i++){
    cosa += (a[i]-b[i])*(c[i]-b[i]);
    absu += (a[i]-b[i])*(a[i]-b[i]);
    absv += (c[i]-b[i])*(c[i]-b[i]);
  }
  cosa /= sqrt(absu*absv);
  cosa = (float)(acos(cosa)*180.0/PI);
  
  return(cosa);
} 

/*******************************************************************************
 * SUBROUTINE dihang calculates the torsional angle between 4 consecutive atoms
 *******************************************************************************/
float dihang(float a[],float b[], float c[], float d[]){
  float t[3],u[3],w[3],m[3],n[3],v[3];
  float absm,absn,absv,absu,costheta,theta,q;
  int i;

  for(i=0;i<=2;i++){
    t[i]=a[i]-b[i];
    u[i]=c[i]-b[i];
    w[i]=d[i]-b[i];
  }

  m[0]=t[1]*u[2]-t[2]*u[1];
  m[1]=t[2]*u[0]-t[0]*u[2];
  m[2]=t[0]*u[1]-t[1]*u[0];

  n[0]=w[1]*u[2]-w[2]*u[1];
  n[1]=w[2]*u[0]-w[0]*u[2];
  n[2]=w[0]*u[1]-w[1]*u[0];

  absm = sqrt(m[0]*m[0]+m[1]*m[1]+m[2]*m[2]);
  absn = sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);

  costheta = (m[0]*n[0]+m[1]*n[1]+m[2]*n[2])/(absm*absn);
 
  theta = acos(costheta);

  v[0]=m[1]*n[2]-m[2]*n[1];
  v[1]=m[2]*n[0]-m[0]*n[2];
  v[2]=m[0]*n[1]-m[1]*n[0];
  
  absv = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  absu = sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);

  q=(v[0]*u[0]+v[1]*u[1]+v[2]*u[2])/(absv*absu);

  theta = (float)(q*theta*180.0/PI);

  return(theta);
}
