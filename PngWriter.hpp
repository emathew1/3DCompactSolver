#ifndef PNG_WRITER_HPP
#define PNG_WRITER_HPP

// ======================================================================
// PngWriter.hpp
// 
// This class can write a png image using libpng routines (i.e. png.h).
//
// Sample usage:
//
//    PngWriter png(nx,ny);
//
//    // then one or more calls to...
//    
//    png.set(i, j, red, green, blue); // 0 <= red,green,blue <= 255
// 
//    png.write("myfile.png");
//
//    // at this point you can change the image and write again...
//
//    png.set(i, j, red, green, blue);
//
//    png.write("myfile2.png");
// 
//
// History:
//  authors: Frank Ham & Phuc Quang - July 2013
// ======================================================================

// ======================================================================
// Copyright (c) 2013 Frank Ham and Phuc Quang
//
// License: MIT (http://opensource.org/licenses/MIT)
//
// use for any purpose, commercial or otherwise.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.
// ======================================================================

#include <png.h>

class PngWriter {
  
private:
  
  unsigned char (*buffer)[3]; // 0 <= r,g,b < 255 
  int nx,ny;
  
 public:
  
  PngWriter(const int width,const int height) {
    
    nx = width;
    ny = height;
    buffer = new unsigned char[nx*ny][3];

    // fill buffer with a "nice" cyan [73,175,205] -- eventually a designer should choose this ;)
    for (int i = 0; i < nx*ny; ++i) {
      buffer[i][0] = 73;
      buffer[i][1] = 175;
      buffer[i][2] = 205;
    }

  }

  ~PngWriter() {

    delete[] buffer;

  }
  
  void set(const int i,const int j,const unsigned char r,const unsigned char g,const unsigned char b) {
    // recall that for png files, the pixels are ordered from the top left, so modify 
    // this set routine's passed j so that zero is at the bottom left...
    buffer[(ny-j-1)*nx+i][0] = r;
    buffer[(ny-j-1)*nx+i][1] = g;
    buffer[(ny-j-1)*nx+i][2] = b;
  }
  
  void write(const char * filename) {

    // note: we completely skip any error handling treatment here for simplicity.
    
    FILE * fp = fopen(filename,"wb");
    if (!fp) {
      std::cout << "Warning: could not open png file: " << filename << std::endl;
      return;
    }
    
    png_structp png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING,
						  (png_voidp)NULL,NULL,NULL);
    if (!png_ptr) {
      fclose(fp);
      std::cout << "Warning: could not create png_ptr" << std::endl;
      return;
    }
    
    png_infop info_ptr = png_create_info_struct(png_ptr);
    if (!info_ptr) {
      png_destroy_write_struct(&png_ptr,(png_infopp)NULL);
      std::cout << "Warning: could not create info_ptr" << std::endl;
      return;
    }
    
    png_init_io(png_ptr, fp);

    png_set_IHDR(png_ptr, info_ptr,
		 nx, ny,             // width, height
		 8,                  // bits per pixel -- 16 does not work with blockbuster
		 PNG_COLOR_TYPE_RGB, // non-alpha options are PNG_COLOR_TYPE_RGB,PNG_COLOR_TYPE_GRAY,
		 PNG_INTERLACE_NONE,
		 PNG_COMPRESSION_TYPE_DEFAULT,
		 PNG_FILTER_TYPE_DEFAULT);

    // Some bits per pixel notes: 16 does not work with blockbuster, and there are also 
    // issues with PNG_COLOR_TYPE_GRAY interpretation, so stick to 8 and PNG_COLOR_TYPE_RGB
    // for now. Note that if you do use 16, pay attention to MSB/LSB order. Endian is 
    // flipped on my linux workstation...
    
    png_write_info(png_ptr, info_ptr);
    
    // set up row pointers to point into the raw image data in buffer...
    
    png_byte * row_pointers[ny];
    for (int i = 0; i < ny; ++i)
      row_pointers[i] = (png_byte*)(buffer + i*nx);
    
    png_write_image(png_ptr, row_pointers);

    png_write_end(png_ptr, NULL);

    png_destroy_write_struct(&png_ptr, &info_ptr);

    fclose(fp);

    // leave the image data in buffer to allow the calling process to change it and/or
    // write another image. It is deleted in the destructor.

  }
  
};

#endif
