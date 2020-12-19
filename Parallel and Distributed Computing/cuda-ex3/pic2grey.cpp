#include <stdio.h>
#include <math.h>
#include "pic2grey.h"

void cpu_pic2grey(unsigned char * Pout, unsigned char * Pin, int width, int height) {
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      int grey_offset = i * width + j;
      int rgb_offset = grey_offset * CHANNELS;
      unsigned char r = Pin[rgb_offset + 0]; 
      unsigned char g = Pin[rgb_offset + 1]; 
      unsigned char b = Pin[rgb_offset + 2]; 
      Pout[grey_offset] = 0.21f * r + 0.71f * g + 0.07f * b;
    } 
  }
}
