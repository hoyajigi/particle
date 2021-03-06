#ifndef __CL_UTIL_H__
#define __CL_UTIL_H__

#include <CL/cl.h>

typedef struct {
  cl_device_id dev;
  cl_event *event;
} bp_data_t;

/* OpenCL utility functions */
cl_device_type get_device_type();
void print_device_name(cl_device_id dev);
char *get_source_code(const char *file_name);
void print_build_log(cl_program program, cl_device_id dev);
void build_program_callback(cl_program program, void *user_data);

#define CHECK_ERROR(err)                                              \
  if (err != CL_SUCCESS) {                                            \
    fprintf(stderr, "[%s:%d] ERROR: %d\n", __FILE__, __LINE__, err);  \
    exit(EXIT_FAILURE);                                               \
  }

#endif //__CL_UTIL_H__
