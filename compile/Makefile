

library = tackle.so
prefix = PREFIX
sharedir = ${prefix}/share
testdir = ${prefix}/test

NVCC = /usr/local/cuda/bin/nvcc
G++ = g++
gprof = -pg
gdb = -g


cpp_files = particles \
            .


CPPFLAGS = -std=c++11 -fopenmp -O3 -funroll-loops -DNDEBUG -Wall -Wno-unknown-pragmas -fPIC $(gprof) 
CUDAFLAGS = -gencode=arch=compute_20,code=compute_20 -Xcompiler \"-fopenmp\",\"-fPIC\",\"-O3\",\"-funroll-loops\",\"-DNDEBUG\",\"-Wall\",\"-Wno-unknown-pragmas\" $(gprof) 
LDFLAGS = -fPIC -fopenmp -O3 -funroll-loops -DNDEBUG -Wall -Wno-unknown-pragmas -shared -Wl,-soname,$(library)
LINKLIB = -lboost_filesystem -lboost_python -lboost_program_options -lboost_signals -lboost_iostreams -lboost_system -lpthread -lpython2.7 -lutil -ldl /usr/local/cuda/lib64/libcudart.so /usr/local/cuda/lib64/libcufft.so -Wl,-rpath,/usr/local/cuda/lib64


others := /usr/include/python2.7 /usr/local/cuda/include /usr/include
INCLUDE_HEADS_FIELS = $(addprefix -I,$(cpp_files) $(others))
vpath %.cc ${cpp_files} ${others}


%.cc.o : %.cc
	@$(G++) -c $< $(CPPFLAGS) -o $@ $(INCLUDE_HEADS_FIELS) #$(CPP_DEFINES)
	@echo "C++ compiling " $<
	
%.u.o : %.cu 
	@$(NVCC) -c $< $(CUDA_DEFINES) $(CUDAFLAGS) -o $@ $(INCLUDE_HEADS_FIELS)
	@echo "NVCC compiling " $<


filecu := $(addsuffix /*.cu,$(cpp_files))
filecc := $(addsuffix /*.cc,$(cpp_files))
cudasource := $(wildcard $(filecu))
cppsource := $(wildcard $(filecc))
cudasource_obj := $(cudasource:%.cu=%.u.o)
cppsource_obj := $(cppsource:%.cc=%.cc.o)
OBJECTS := $(cppsource_obj) $(cudasource_obj)


$(library) : $(OBJECTS)
	@$(G++) $(LDFLAGS) -o $@ $(OBJECTS) $(LINKLIB) 

clean :
	@$(RM) $(OBJECTS)
	@echo "clean .o successfully!"


