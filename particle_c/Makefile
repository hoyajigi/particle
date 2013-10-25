TARGET = particle
OBJS   = particle.o

CC       = g++
INCLUDES = -I.
CFLAGS   = -O2 -Wall $(INCLUDES) -DDEBUG
LDFLAGS  = -lm


all: $(TARGET)

$(TARGET): $(OBJS)
	$(CC) $(LDFLAGS) $(OBJS) -o $@ 

.c.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET) *.stdout *.stderr

sc:
	rm -f *.stdout *.stderr
