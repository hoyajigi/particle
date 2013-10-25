#ifndef __PARTICLE_UTIL_H__
#define __PARTICLE_UTIL_H__

#define PI               3.14159265358979323846f

#include <CL/cl.h>
#include <math.h>

typedef cl_float4 float4;
typedef cl_int4 int4; 
typedef unsigned int u32;

inline float4 make_float4(float x, float y, float z, float w = 0.f) 
{
  float4 tmp;
  tmp.x = x;
  tmp.y = y;
  tmp.z = z;
  tmp.w = w;
  return tmp;
}

inline int4 make_int4(int x, int y, int z, int w = 0) 
{
  int4 tmp;
  tmp.x = x;
  tmp.y = y;
  tmp.z = z;
  tmp.w = w;
  return tmp;
}

inline float4 mul_f4(float4 a, float4 b)
{
  a.x = a.x * b.x;
  a.y = a.y * b.y;
  a.z = a.z * b.z;
  a.w = a.w * b.w;
  return a;
}

inline float4 add_f4(float4 a, float4 b)
{
  a.x = a.x + b.x;
  a.y = a.y + b.y;
  a.z = a.z + b.z;
  a.w = a.w + b.w;
  return a;
}

inline float4 sub_f4(float4 a, float4 b)
{
  a.x = a.x - b.x;
  a.y = a.y - b.y;
  a.z = a.z - b.z;
  a.w = a.w - b.w;
  return a;
}

inline float4 div_f4(float4 a, float4 b)
{
  a.x = a.x / b.x;
  a.y = a.y / b.y;
  a.z = a.z / b.z;
  a.w = a.w / b.w;
  return a;
}

typedef struct _ConstBuffer
{
  float4 m_g;
  int m_numParticles;
  float m_dt;
  float m_scale;
  float m_e;

  int4 m_nCells;
  float4 m_spaceMin;
  float m_gridScale;
}ConstBuffer;

typedef struct _ConstBufferGrid
{
	float4 m_max;
	float4 m_min;
	int4 m_nCells;
	float m_gridScale;
	u32 m_maxParticles;
} ConstBufferGrid;

float randRange(const float minV, const float maxV)
{
	float r = (rand()%10000)/10000.f;
	float range = maxV - minV;
	return (float)(minV + r*range);
}

int4 ugConvertToGridCrd(const float4 pos, float gridScale)
{
	int4 g;
	g.x = floor(pos.x*gridScale);
	g.y = floor(pos.y*gridScale);
	g.z = floor(pos.z*gridScale);
	return g;
}

int ugGridCrdToGridIdx(const int4 g, int nCellX, int nCellY, int nCellZ)
{
	return g.x+g.y*nCellX+g.z*nCellX*nCellY;
}

__inline
float dot3F4(float4 a, float4 b)
{
  float sum = 0;
	float4 a1 = make_float4(a.x, a.y, a.z,0.f);
	float4 b1 = make_float4(b.x, b.y, b.z,0.f);

  sum += a1.x * b1.x;
  sum += a1.y * b1.y;
  sum += a1.z * b1.z;
	return sum;
}

__inline
float dot3w1(const float4 point, const float4 eqn)
{
	return dot3F4(point,eqn) + eqn.w;
}

__inline
float4 calcForce1(const float4 x_i, const float4 x_j, const float4 v_i, const float4 v_j, 
				 float r_i, float r_j, float m_i, float m_j, float dt,
				 float sCoeff, float dCoeff)
{
	float4 f = make_float4(0,0,0,0);
	float4 x_ij = sub_f4(x_j, x_i);
	float dist = sqrt(dot3F4(x_ij, x_ij));


	if( dist < r_i+r_j )
	{
//		f -= sCoeff*(r_i+r_j-dist)*x_ij*div(1.f,dist);
    float4 tmp4 = mul_f4(x_ij, make_float4(1.f/dist, 1.f/dist, 1.f/dist, 1.f/dist));
    float tmp = sCoeff*(r_i+r_j-dist);
    tmp4 = mul_f4(make_float4(tmp, tmp, tmp, tmp), tmp4);
    f = sub_f4(f, tmp4);

//		f += dCoeff*(v_j - v_i);
    tmp4 = sub_f4(v_j, v_i);
    tmp4 = mul_f4(make_float4(dCoeff, dCoeff, dCoeff, dCoeff), tmp4);
    f = add_f4(f, tmp4);
	}
	return f;
}

#endif // __PARTICLE_UTIL_H__
