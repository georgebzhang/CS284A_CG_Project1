#ifndef CGL_DRAWREND_H
#define CGL_DRAWREND_H

#include "CGL/CGL.h"
#include "CGL/renderer.h"
#include "CGL/color.h"
#include <vector>
#include <cstring>
#include "GLFW/glfw3.h"
#include "svg.h"

namespace CGL {

class DrawRend : public Renderer {
 public:
  DrawRend(std::vector<SVG*> svgs_):
  svgs(svgs_), current_svg(0)
  {}

  ~DrawRend( void );

  // inherited Renderer interface functions
  void init();
  void render();
  void resize( size_t w, size_t h );
  std::string name() { return "Draw"; }
  std::string info();
  void cursor_event( float x, float y );
  void scroll_event( float offset_x, float offset_y );
  void mouse_event( int key, int event, unsigned char mods );
  void keyboard_event( int key, int event, unsigned char mods );

  void set_gl(bool gl_) { gl = gl_; }

  // write current pixel buffer to disk
  void write_screenshot();

  // write only framebuffer to disk
  void write_framebuffer();

  // drawing functions
  void redraw();
  void draw_pixels();
  void draw_zoom();

  // view transform functions
  void view_init();
  void set_view(float x, float y, float span);
  void move_view(float dx, float dy, float scale);

  // rasterize a point
  void rasterize_point( float x, float y, Color color );

  // rasterize a line
  void rasterize_line( float x0, float y0,
                       float x1, float y1,
                       Color color);

  // returns 0.0 if (x,y) is on the line T = (x1-x0, y1-y0), >0.0 if the dot product between
  // V = (x-x0, y-y0) and N = perp(T) is positive, and <0.0 otherwise
  float dot_prod(float& x, float& y, float& x0, float& y0, float& x1, float& y1);

  bool point_in_tri(float& x, float& y, float& x0, float& y0, float& x1, float& y1, float& x2, float& y2);

  bool bary_point_in_tri(double a, double b, double y);

  bool bary_point_in_tri(Vector3D bc);

  Vector3D bary_coords(float cx, float cy, float x0, float y0, float x1, float y1, float x2, float y2);

  // rasterize a triangle
  void rasterize_triangle( float x0, float y0,
                           float x1, float y1,
                           float x2, float y2,
                           Color color, Triangle *tri = NULL );



private:
  // Global state variables for SVGs, pixels, and view transforms
  std::vector<SVG*> svgs; size_t current_svg;
  std::vector<Matrix3x3> svg_to_ndc;
  float view_x, view_y, view_span;

  Matrix3x3 ndc_to_screen;

  std::vector<unsigned char> framebuffer;
  size_t width, height;

  // UI state info
  float cursor_x; float cursor_y;
  bool left_clicked;
  int show_zoom;
  int sample_rate;

  PixelSampleMethod psm;
  LevelSampleMethod lsm;

  bool gl;

  // Intuitively, PixelColorStorage contains r g b intensities as unsigned char values
  typedef std::vector<unsigned char> PixelColorStorage;

  // Intuitively, a sample buffer instance is a pixel,
  // or (samples_per_side x samples_per_side) sub-pixels.
  // samples_per_side 1, 2, 3, or 4 depending on sample_rate 1, 4, 9, or 16
  struct SampleBuffer {
    std::vector<std::vector<PixelColorStorage> > sub_pixels;
    size_t samples_per_side;

    SampleBuffer(size_t sps): samples_per_side(sps) {
      clear();
    }
    
    // Fill the subpixel at i,j with the Color c
    void fill_color(int i, int j, Color c) {
      PixelColorStorage &p = sub_pixels[i][j];
      // Part 1: Overwrite PixelColorStorage p using Color c.
      //         Pay attention to different data types.
	  float temp[] = { (unsigned char)255.f*c.r, (unsigned char)255.f*c.g, (unsigned char)255.f*c.b };
	  p.assign(temp, temp + 3);
      return;
    }

	// fills all subpixels (sample_rate = 1, 4, 9, or 16) of a pixel with Color c
    void fill_pixel(Color c) {
      for (int i = 0; i < samples_per_side; ++i)
        for (int j = 0; j < samples_per_side; ++j)
          fill_color(i, j, c);
    }

    Color get_pixel_color() {
      //return Color(sub_pixels[0][0].data());
      // Part 2: Implement get_pixel_color() for supersampling.
		float r_total = 0, g_total = 0, b_total = 0;

		for (int i = 0; i < samples_per_side; ++i) {
			for (int j = 0; j < samples_per_side; ++j) {
				r_total += (float)sub_pixels[i][j][0];
				g_total += (float)sub_pixels[i][j][1];
				b_total += (float)sub_pixels[i][j][2];
			}
		}

		float r, g, b;

		// sample rate sr = samples_per_side ^ 2
		int sr = pow(samples_per_side, 2);

		r = r_total / sr / 255.f;
		g = g_total / sr / 255.f;
		b = b_total / sr / 255.f;
		
		return Color(r, g, b);
    }
    
    void clear() {
      if (sub_pixels.size() == samples_per_side) {
        for (int i = 0; i < samples_per_side; ++i)
          for (int j = 0; j < samples_per_side; ++j)
            sub_pixels[i][j].assign(3, (unsigned char)255);
        return;
      }

      sub_pixels.clear();
      PixelColorStorage white = std::vector<unsigned char>(3, 255);
      std::vector<PixelColorStorage> row;
      row.reserve(samples_per_side);
      for (int i = 0; i < samples_per_side; ++i)
        row.push_back(white);
      sub_pixels.reserve(samples_per_side);
      for (int i = 0; i < samples_per_side; ++i)
        sub_pixels.push_back(row);
    }
  };

  std::vector<std::vector<SampleBuffer> > samplebuffer;

  // This function takes the collected sub-pixel samples and
  // combines them together to fill in the framebuffer in preparation
  // for posting pixels to the screen.
  void resolve() {
    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = samplebuffer[y][x].get_pixel_color();
        for (int k = 0; k < 3; ++k) {
          framebuffer[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }
  }


};

} // namespace CGL

#endif // CGL_DRAWREND_H
