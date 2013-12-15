#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include "cl_util.h"


// Get the device type from the environmental variable.
cl_device_type get_device_type() {
  char *dtype = getenv("CL_DEV_TYPE");
  if (dtype) {
    if (strcasecmp(dtype, "cpu") == 0) {
      return CL_DEVICE_TYPE_CPU;
    } else if (strcasecmp(dtype, "gpu") == 0) {
      return CL_DEVICE_TYPE_GPU;
    } else if (strcasecmp(dtype, "acc") == 0) {
      return CL_DEVICE_TYPE_ACCELERATOR;
    }
  }
  return CL_DEVICE_TYPE_DEFAULT;
}


void print_device_name(cl_device_id dev) {
  cl_int err;
  size_t name_size;

  clGetDeviceInfo(dev, CL_DEVICE_NAME, 0, NULL, &name_size);

  char *dev_name = (char *)malloc(name_size + 1);
  err = clGetDeviceInfo(dev, CL_DEVICE_NAME,
                        name_size, dev_name, NULL);
  CHECK_ERROR(err);

  printf("Device: %s\n", dev_name);

  free(dev_name);
}


char *get_source_code(const char *file_name) {
  FILE *file = fopen(file_name, "r");
  if (file == NULL) {
    fprintf(stderr, "[%s:%d] ERROR: Failed to open %s\n",
            __FILE__, __LINE__, file_name);
    exit(EXIT_FAILURE);
  }

  fseek(file, 0, SEEK_END);
  size_t length = (size_t)ftell(file);
  rewind(file);

  char *source_code = (char *)malloc(length + 1);
  if (fread(source_code, length, 1, file) != 1) {
    fclose(file);
    free(source_code);

    fprintf(stderr, "[%s:%d] ERROR: Failed to read %s\n",
            __FILE__, __LINE__, file_name);
    exit(EXIT_FAILURE);
  }
  fclose(file);

  source_code[length] = '\0';
  return source_code;
}


void print_build_log(cl_program program, cl_device_id dev) {
  size_t log_size;
  clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
                        0, NULL, &log_size);
  char *log = (char *)malloc(log_size + 1);
  clGetProgramBuildInfo(program, dev, CL_PROGRAM_BUILD_LOG,
                        log_size, log, NULL);
  fprintf(stderr, "\n---------- BUILD LOG ----------\n%s\n", log);
  fprintf(stderr, "---------------------------------------\n");
  free(log);
}


void build_program_callback(cl_program program, void *user_data) {
  cl_int err;
  cl_build_status build_status;
  bp_data_t *bp_data = (bp_data_t *)user_data;

  // Check the build status.
  err = clGetProgramBuildInfo(program, bp_data->dev,
                              CL_PROGRAM_BUILD_STATUS,
                              sizeof(cl_build_status),
                              &build_status, NULL);
  CHECK_ERROR(err);
  if (build_status != CL_BUILD_SUCCESS) {
    print_build_log(program, bp_data->dev);
    exit(EXIT_FAILURE);
  }

  // Set the event status
  err = clSetUserEventStatus(*(bp_data->event), CL_COMPLETE);
  CHECK_ERROR(err);
}
