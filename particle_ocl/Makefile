TARGET = particle
OBJS   = particle.o cl_util.o timers.o

CC       = g++
INCLUDES = #-I/home/sc56/src/AMD_Samples/include/
CFLAGS   = -std=c99 -O4 -Wall $(INCLUDES)
LDFLAGS  = -lm -lOpenCL #-L/opt/AMDAPP/lib/x86_64

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@ 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET) *.stdout *.stderr

sc:
	rm -f *.stdout *.stderr
