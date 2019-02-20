#ifndef CGL_TEXTURE_H
#define CGL_TEXTURE_H

#include <vector>
#include "CGL/CGL.h"
#include "CGL/color.h"
#include "CGL/vector2D.h"

#include <iostream>
using namespace std;

namespace CGL {

	// P_NEAREST: nearest-pixel sampling, P_LINEAR: bilinear sampling
typedef enum PixelSampleMethod { P_NEAREST = 0, P_LINEAR = 1 } PixelSampleMethod;
	// L_ZERO: sample from zero-th mipmap (Part 5), L_NEAREST: sample from the nearest mipmap, L_LINEAR: compute weighted sum from two samples (one from each adjacent mipmap)
typedef enum LevelSampleMethod { L_ZERO = 0, L_NEAREST = 1, L_LINEAR = 2 } LevelSampleMethod;

struct SampleParams {
  Vector2D p_uv;
  Vector2D p_dx_uv, p_dy_uv;
  PixelSampleMethod psm;
  LevelSampleMethod lsm;
};

static const int kMaxMipLevels = 14;

struct MipLevel {
	size_t width;
	size_t height;
  // RGB color values
  std::vector<unsigned char> texels;

  // Part 5 and 6
  Color get_texel(int tx, int ty) {
	  float r, g, b;

	  r = (float)texels[3 * (ty * width + tx) + 0] / 255.f;
	  g = (float)texels[3 * (ty * width + tx) + 1] / 255.f;
	  b = (float)texels[3 * (ty * width + tx) + 2] / 255.f;

	  return Color(r, g, b);
  }
};

struct Texture {
  size_t width;
  size_t height;
  std::vector<MipLevel> mipmap;

  void init(const vector<unsigned char>& pixels, const size_t& w, const size_t& h) {
    width = w; height = h;

    // A fancy C++11 feature. emplace_back constructs the element in place,
    // and in this case it uses the new {} list constructor syntax.
    mipmap.emplace_back(MipLevel{width, height, pixels});

    generate_mips();
  }

  // Generates up to kMaxMipLevels of mip maps. Level 0 contains
  // the unfiltered original pixels.
  void generate_mips(int startLevel = 0);

  Color sample(const SampleParams &sp);
  float get_level(const SampleParams &sp);

  Color sample_nearest(Vector2D uv, int level = 0);

  Color lerp(double s, Color c0, Color c1);

  Color sample_bilinear(Vector2D uv, int level = 0);

  // does bilinear sampling if int level is a whole number, otherwise does trilinear sampling
  Color sample_linear(Vector2D uv, int level = 0);
};

}

#endif // CGL_TEXTURE_H
