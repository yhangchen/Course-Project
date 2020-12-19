# Before making, Please module add 
# 	gcc/6.5.0
# 	cuda/9.0
# 	ffmpeg/4.1
# 	opencv/4.0.0

CXX = g++
CXXFLAGS = -std=c++11

NVCC = nvcc
NVCC_FLAGS = -std=c++11 -O3 -Wno-deprecated-gpu-targets

INCLUDES= -I/mnt/lustrefs/softwares/opencv/4.0.0/include/opencv4
LIBRARIES = -L/mnt/lustrefs/softwares/cuda/9.0/lib64 \
						-L/mnt/lustrefs/softwares/opencv/4.0.0/lib64 -lopencv_core \
						-lcudart -lopencv_imgproc -lopencv_imgcodecs -lopencv_highgui

pic2grey: main.o cpu_pic2grey.o gpu_pic2grey.o 
	$(CXX) $^ -o $@ $(LIBRARIES)

main.o: main.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $^ -o $@

cpu_pic2grey.o: pic2grey.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

gpu_pic2grey.o: pic2grey.cu
	$(NVCC) $(NVCC_FLAGS) -c $^ -o $@

clean:
	rm -f *.o pic2grey
