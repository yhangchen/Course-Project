/* 
  Foundations of Parallel and Distributed Computing, Falls 2019.
  Instructor: Prof. Chao Yang @ Peking University.
  This code shows how to convert a colorful pic to greyscale. 
*/

#include <iostream>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include "pic2grey.h"

using namespace std;
using namespace cv;

int main(){

  // Load the original image into Pin
  //Mat src_img = imread("images/cat_small.jpg", IMREAD_COLOR);
  Mat src_img = imread("images/cat_large.jpg", IMREAD_COLOR);
  int height = src_img.rows;
  int width = src_img.cols;
  cout << "height: " << height << " " << "width: " << width << endl;
  unsigned char *Pin = new unsigned char[height * width * CHANNELS];
  for(int i = 0; i < height; i++){
    const uchar *img_data = src_img.ptr<uchar>(i);
    for(int j = 0; j < width * CHANNELS; j++) {
      Pin[i * width * CHANNELS + j] = img_data[j];
    }
  }
  unsigned char* Pout = new unsigned char[height * width];

  // Process the image by CPU
  long long cpu_start_time = start_timer();
  cpu_pic2grey(Pout, Pin, width, height);
  long long cpu_time = stop_timer(cpu_start_time, (char *)"CPU");

  // Store the grey image processed by CPU 
  Mat cpu_grey_img = Mat(height, width, CV_8UC1);
  memcpy(cpu_grey_img.data, Pout, width * height);
  imwrite("cpu_grey_img.jpg", cpu_grey_img);

  // Process the image by GPU
  long long gpu_start_time = start_timer();
  gpu_pic2grey(Pout, Pin, width, height);
  long long gpu_time = stop_timer(gpu_start_time, (char *)"GPU");

  // Store the grey image processed by GPU 
  Mat gpu_grey_img = Mat(height, width, CV_8UC1);
  memcpy(gpu_grey_img.data, Pout, width * height);
  imwrite("gpu_grey_img.jpg", gpu_grey_img);
#if 0
  // Process the image by OpenCV 
  Mat cv_grey_img;
  long long cv_start_time = start_timer();
  cvtColor(src_img, cv_grey_img, COLOR_BGR2GRAY);
  long long cv_time = stop_timer(cv_start_time, (char *)"\nOpenCV");
  imwrite("cv_grey_img.jpg", cv_grey_img);
#endif
}
