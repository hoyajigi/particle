#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <strings.h>
#include <time.h>
#include <CL/cl.h>
#include "cl_util.h"
#include "timers.h"
#include "particle_util.h"



#define ERROR(err) fprintf(stderr, "[%s:%d] ERROR: %s\n",__FILE__,__LINE__,err);exit(EXIT_FAILURE);


#define MAX_IDX_PER_GRID  (1024*1024*8) // default: 6

int main(int argc, char **argv) {
  /* Host data structures */
  cl_platform_id   *platforms;
  cl_uint          num_platforms;
  cl_device_type   dev_type = CL_DEVICE_TYPE_DEFAULT;
  cl_device_id     dev;
  cl_context       context;
  // NOTE : You might have multiple cmd_queue but whatever
  cl_command_queue cmd_queue;
  cl_program       program;
  cl_kernel        kernel;
  // TODO : define your variables
  cl_mem           posIn,velIn,cGrid,cGridCounter,forceOut;
  cl_int           err;
  cl_uint          num_dev = 0;
  cl_event         ev_bp;
 

#ifdef DEBUG
  double checksum = .0;
#endif
  int seed = 777;






 
  // TODO : 
//  size_t lws[2]={64,4};
//  size_t gws[2]={1024,1024};


  int i;


  timer_init();

  srand(time(NULL));


  // Platform
  err = clGetPlatformIDs(0, NULL, &num_platforms);
  CHECK_ERROR(err);
  if(num_platforms == 0) {
    ERROR("No OpenCl platform");
  }
  printf("Number of platforms: %u\n",num_platforms);
  platforms = (cl_platform_id *)malloc(sizeof(cl_platform_id) * num_platforms);
  err = clGetPlatformIDs(num_platforms,platforms,NULL);
  CHECK_ERROR(err);

  //Device
  for(i=0;i<num_platforms;i++) {
    // FIXME : something wrong
    err = clGetDeviceIDs(platforms[i],dev_type,1,&dev,&num_dev);
    if(err != CL_DEVICE_NOT_FOUND) CHECK_ERROR(err);
    if(num_dev == 1) break;
  }
  if(num_dev<1) {
    ERROR("No device");
  }
  
  // Print the device name.
  size_t name_size;
  clGetDeviceInfo(dev, CL_DEVICE_NAME, 0, NULL, &name_size);
  char *dev_name = (char *)malloc(name_size + 1);
  err = clGetDeviceInfo(dev,CL_DEVICE_NAME,name_size,dev_name,NULL);
  CHECK_ERROR(err);
  printf("Device: %s\n",dev_name);
  free(dev_name);

  // Context
  context = clCreateContext(NULL, 1, &dev, NULL, NULL, &err);
  CHECK_ERROR(err);

  // Command queue
  cmd_queue = clCreateCommandQueue(context, dev, 0, &err);
  CHECK_ERROR(err);

  // Create a program
  // TODO : Get source code in your favor
  char * source_code=get_source_code("hello.cl");
  
  
  size_t source_len=strlen(source_code);
  program = clCreateProgramWithSource(context, 1, (const char **)&source_code, &source_len, &err);
  CHECK_ERROR(err);

  // Callback data for clBuildProgram
  ev_bp=clCreateUserEvent(context,&err);
  CHECK_ERROR(err);
  bp_data_t bp_data;
  bp_data.dev=dev;
  bp_data.event=&ev_bp;

  // Build the program.
  err = clBuildProgram(program, 1, &dev, NULL, build_program_callback, &bp_data);
  if (err != CL_SUCCESS) {
    // Print the build log.
    size_t log_size;
    clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
        0, NULL, &log_size);
    char *log = (char *)malloc(log_size + 1);
    clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
        log_size, log, NULL);
    fprintf(stderr,"\n");
    fprintf(stderr,"---------- BUILD LOG ----------\n");
    fprintf(stderr,"%s\n",log);
    fprintf(stderr,"-------------------------------\n");
    free(log);

    CHECK_ERROR(err);
  }
  







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
























  
  // Buffers
  // TODO: make and buffers
  /*
  clCreateBuffer(context, CL_MEM_READ_ONLY|CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
      sizeof(float) * N, A, &err);
  err=clEnqueueWriteBuffer(cmd_queue, m_array,CL_FALSE,0,size*sizeof(int),array,0,NULL,NULL);
  */
  CHECK_ERROR(err);

  clWaitForEvents(1,bp_data.event);
  
  // Kernel
  kernel = clCreateKernel(program,"blank",&err);
  CHECK_ERROR(err);

  clFinish(cmd_queue);

  // TODO : Set the arguments.
  //err=clSetKernelArg(kernel,0,sizeof(),);

  // Enqueue the kernel.
  //err=clEnqueueNDRangeKernel(cmd_queue,kernel,1,NULL,gws,lws,0,NULL,NULL);
  CHECK_ERROR(err);

  // Read the result.
  /*
  err = clEnqueueReadBuffer(cmd_queue,
      mem_C,
      CL_TRUE, 0,
      sizeof(float) * N,
      C,
      0, NULL, NULL);
  */
  CHECK_ERROR(err);

  // Release
  clReleaseEvent(ev_bp);
  //clReleaseMemObject();
  clReleaseKernel(kernel);
  clReleaseProgram(program);
  clReleaseCommandQueue(cmd_queue);
  clReleaseContext(context);
  free(platforms);

  return EXIT_SUCCESS;
}
