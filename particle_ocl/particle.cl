typedef unsigned int u32;
typedef unsigned short u16;
typedef unsigned char u8;

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

#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable

#define GET_GROUP_IDX get_group_id(0)
#define GET_LOCAL_IDX get_local_id(0)
#define GET_GLOBAL_IDX get_global_id(0)

#define MAX_IDX_PER_GRID  (1024*1024*8) // default: 6

// Vector

#define make_float4 (float4)
#define make_int4 (int4)
#define make_float2 (float2)
#define make_int2 (int2)
#define max2 max2
#define min2 min

__inline
float dot3F4(float4 a, float4 b)
{
	float4 a1 = make_float4(a.x, a.y, a.z,0.f);
	float4 b1 = make_float4(b.x, b.y, b.z,0.f);
	return dot(a1,b1);
}

__inline
float dot3w1(const float4 point, const float4 eqn)
{
	return dot3F4(point,eqn) + eqn.w;
}

__inline div(float a,float b)
{
	return native_divide(a,b);
}

__inline
float4 calcForce1(const float4 x_i, const float4 x_j, const float4 v_i, const float4 v_j, 
				 float r_i, float r_j, float m_i, float m_j, float dt,
				 float sCoeff, float dCoeff)
{
	float4 f = make_float4(0,0,0,0);
	float4 x_ij = x_j-x_i;
	float dist = native_sqrt(dot3F4(x_ij, x_ij));


	if( dist < r_i+r_j )
	{
		f -= sCoeff*(r_i+r_j-dist)*x_ij*div(1.f,dist);
		f += dCoeff*(v_j - v_i);
	}
	return f;
}

// Grid related part
int4 ugConvertToGridCrd(const float4 pos, float gridScale)
{
	int4 g;
	g.x=floor(pos.x*gridScale);
	g.y=floor(pos.y*gridScale);
	g.z=floor(pos.z*gridScale);
	return g;
}

int ugGridCrdToGridIdx(const int4 g,int nCellX,int nCellY,int nCellZ)
{
	return g.x+g.y*nCellX+g.z*nCellX*nCellY;
}

__kernel void CollideGridKernel(__global float4* posIn,
								__global float4* velIn,
								__global int* cGrid,
								__global int* cGridCounter,
								__global float4* forceOut,
								ConstBuffer cb)
{
	
}