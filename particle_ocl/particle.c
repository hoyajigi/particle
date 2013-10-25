#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "particle_util.h"
#define MAX_IDX_PER_GRID  (1024*1024*8) // default: 6

int main(int argc, char **argv) {
#ifdef DEBUG
  double checksum = .0;
#endif
  int seed = 777;

//===----------------------------------------------------------------------===//
// Data structure construction
//===----------------------------------------------------------------------===//
  srand(seed); // use time value to get real random number
  float s = 1.f;
	float density = 100.f;
  int m_numSParticles = 256;
  float m_scale = 1.f; //50.f;
  float cellSize = m_scale*2.f;
  float m_gridScale = 1.f/(cellSize);

  float4 m_min = make_float4(-s, -s, 0, 0 );
  m_min = mul_f4(m_min, make_float4(m_scale,m_scale,m_scale,m_scale));

  float4 m_max = make_float4(s, s, 0, 0 );
  m_max = mul_f4(m_max, make_float4(m_scale,m_scale,m_scale,m_scale));

  float4 extent = sub_f4(m_max, m_min);
  int nx = ((int)(extent.x/cellSize)) + 1;
  int ny = ((int)(extent.y/cellSize)) + 1;
  int nz = ((int)(extent.z/cellSize)) + 1;
  int4 nCells = make_int4(nx, ny, nz);
  int numCells = nx * ny * nz;

  float4 *m_posS = (float4*)calloc(sizeof(float4), m_numSParticles);
  float4 *m_velS= (float4*)calloc(sizeof(float4), m_numSParticles);
  float4 *m_forceS = (float4*)calloc(sizeof(float4), m_numSParticles);
  int* m_grid = (int*)calloc(sizeof(int), (numCells)*MAX_IDX_PER_GRID);
  int* m_gridCounter = (int*)calloc(sizeof(int), numCells);

  // initization
	for(int i=0; i<m_numSParticles; i++)
	{
		m_posS[i] = mul_f4(
        make_float4( randRange(-s,s), randRange(-s,s), 0 ),
        make_float4( m_scale, m_scale, m_scale, m_scale ));
		float r = m_scale;
		m_velS[i] = make_float4(0,0,0,r*r*PI*density);
		m_posS[i].w = r;
	}
	for(int i=0; i<m_numSParticles; i++)
	{
		m_forceS[i] = make_float4(0,0,0,0);
	}

  // Construct random grid
  for(int i=0; i< numCells; i++)
  {
    m_gridCounter[i] = rand()%MAX_IDX_PER_GRID;
  }

  for(int i=0; i< numCells*MAX_IDX_PER_GRID; i++)
  {
    m_grid[i] = rand()%m_numSParticles;
  }
  
	float e = 0.85f;
	float4 g = mul_f4(make_float4(0.f, -9.8f, 0.f, 0.f), 
      make_float4(0.5f, 0.5f, 0.5f, 0.5f));

	float dt = 1.f/60.f;
    dt/=1200.f;

	ConstBuffer cb;
	{
		cb.m_g = g;
		cb.m_dt = dt;
		cb.m_numParticles = m_numSParticles;
		cb.m_scale = m_scale;
		cb.m_e = e;
		cb.m_nCells = nCells;
		cb.m_spaceMin = m_min;
		cb.m_gridScale = m_gridScale;
	}

//===----------------------------------------------------------------------===//
// SS-Collide
//===----------------------------------------------------------------------===//
  float4* posIn = m_posS;
  float4* velIn = m_velS;
  int* cGrid = m_grid;
  int* cGridCounter = m_gridCounter;
  float4* forceOut = m_forceS;

  for (int idx = 0; idx < m_numSParticles; ++idx) {
    if( idx >= cb.m_numParticles ) continue;

    int4 nCells = cb.m_nCells;
    float4 spaceMin = cb.m_spaceMin;
    float gridScale = cb.m_gridScale;
    float dt = cb.m_dt;
    float e = cb.m_e;

    float4 f = make_float4(0,0,0,0);

    float4 x_i = posIn[ idx ];
    float4 v_i = velIn[ idx ];

    float sCoeff, dCoeff;
    {
      float m_i = v_i.w;
      float m_j = v_i.w;
      float m = (m_i*m_j)/(m_i+m_j);
      sCoeff = m/(dt*dt);
      dCoeff =(m*(1.f-e))/dt;
    }

    int4 iGridCrd = ugConvertToGridCrd( sub_f4(x_i, spaceMin), gridScale );

    // 1. collide particles
    int k=0;
    for(int i=-1;i<=1;i++) for(int j=-1;j<=1;j++)// for(int k=-1;k<=1;k++)
    {
      int4 gridCrd = make_int4(iGridCrd.x+i, iGridCrd.y+j, iGridCrd.z+k, 0);

      if( gridCrd.x < 0 || gridCrd.x >= nCells.x
          || gridCrd.y < 0 || gridCrd.y >= nCells.y
          || gridCrd.z < 0 || gridCrd.z >= nCells.z ) continue;

      int gridIdx = ugGridCrdToGridIdx( gridCrd, nCells.x, nCells.y, nCells.z );

      int numElem = cGridCounter[ gridIdx ];
      numElem = fmin(numElem, MAX_IDX_PER_GRID);

      for(int ie=0; ie<numElem; ie++)
      {
        int jIdx = cGrid[MAX_IDX_PER_GRID*gridIdx + ie];
        if( jIdx == idx ) continue;

        float4 x_j = posIn[jIdx];
        float4 v_j = velIn[jIdx];

        f = add_f4(f, 
            calcForce1( x_i, x_j, v_i, v_j, x_i.w, x_j.w, v_i.w, v_j.w, 
              dt, sCoeff, dCoeff ));
      }
    }

    if (idx == 0)

      // 2. collide with boundary
    {
      float sCoeff, dCoeff;
      {
        float m = v_i.w/2.f;
        sCoeff = m/(dt*dt);
        dCoeff = m/dt*(1.f-e);
      }

      float4 planes[4];
      planes[0] = make_float4(0,1,0,cb.m_scale);
      planes[1] = make_float4(-1,0,0,cb.m_scale);
      planes[2] = make_float4(1,0,0,cb.m_scale);
      planes[3] = make_float4(0,-1,0,cb.m_scale);


      for(int j=0; j<4; j++)
      {
        float4 eqn = planes[j];
        float dist = dot3w1( x_i, eqn );
        float r_i = x_i.w;
        if( dist < r_i )
        {
//          f += sCoeff*(r_i-dist)*eqn;
          float tmp = sCoeff*(r_i-dist);
          float4 tmp4 = make_float4(tmp, tmp, tmp, tmp);
          f = add_f4(f, mul_f4(tmp4, eqn));

//          f += dCoeff*(-v_i);
          tmp4 = mul_f4(make_float4(-1.0f, -1.0f, -1.0f, -1.0f), v_i);
          f = add_f4(f, mul_f4(make_float4(dCoeff, dCoeff, dCoeff, dCoeff), tmp4));
        }
      }
    }

    forceOut[ idx ] = f;
  }

//===----------------------------------------------------------------------===//
// S-Integrate
//===----------------------------------------------------------------------===//
  for (int idx = 0; idx < m_numSParticles; ++idx) {
    float4 x = m_posS[idx];
    float4 v = m_velS[idx];


    float4 m_dt4 = make_float4(cb.m_dt, cb.m_dt, cb.m_dt, cb.m_dt);
    float4 tmp = mul_f4(m_forceS[idx], m_dt4);
    float4 tmp2 = make_float4(v.w, v.w, v.w, v.w);
    tmp = div_f4(tmp, tmp2);

    v = add_f4(v, add_f4(tmp, cb.m_g));
    x = add_f4(x, mul_f4(v, m_dt4));

    m_posS[idx] = make_float4(x.x, x.y, x.z, m_posS[idx].w);
    m_velS[idx] = make_float4(v.x, v.y, v.z, m_velS[idx].w);
  }

#ifdef DEBUG
  for (int i = 0; i < 100; i++) {
    checksum += m_posS[rand()%m_numSParticles].x;
  }
  printf("Checksum: %ld (seed:%d)\n", (long)checksum, seed);
#endif

  return 0;
}
