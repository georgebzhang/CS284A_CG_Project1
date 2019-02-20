#include "texture.h"
#include "CGL/color.h"

#include <cmath>
#include <algorithm>

namespace CGL {

Color Texture::sample(const SampleParams &sp) {
  // Parts 5 and 6: Fill this in.
  // Should return a color sampled based on the psm and lsm parameters given

	if (sp.psm == P_LINEAR)
		return sample_linear(sp.p_uv, get_level(sp));
	else if (sp.psm == P_NEAREST)
		return sample_nearest(sp.p_uv, get_level(sp));
	else
		return Color();
}

float Texture::get_level(const SampleParams &sp) {
  // Optional helper function for Parts 5 and 6
	size_t w, h;
	w = this->mipmap[0].width;
	h = this->mipmap[0].height;

	double dudx, dvdx, dudy, dvdy;
	dudx = sp.p_dx_uv[0] - sp.p_uv[0];
	dvdx = sp.p_dx_uv[1] - sp.p_uv[1];
	dudy = sp.p_dy_uv[0] - sp.p_uv[0];
	dvdy = sp.p_dy_uv[1] - sp.p_uv[1];

	dudx *= w;
	dvdx *= h;
	dudy *= w;
	dvdy *= h;

	double L = max(sqrt(pow(dudx, 2) + pow(dvdx, 2)), sqrt(pow(dudy, 2) + pow(dvdy, 2)));

	if (sp.lsm == L_ZERO)
		return 0;
	else if (sp.lsm == L_NEAREST)
		return round(log2(L));
	else if (sp.lsm == L_LINEAR)
		return log2(L);
	else
		return 0;
}

// Returns the nearest sample given a particular level and set of uv coords
Color Texture::sample_nearest(Vector2D uv, int level) {
  // Optional helper function for Parts 5 and 6
  // Feel free to ignore or create your own

	size_t w, h;
	w = this->mipmap[level].width;
	h = this->mipmap[level].height;
	double u, v;
	u = uv[0] * w;
	v = uv[1] * h;
	int tu, tv;
	tu = floor(u);
	tv = floor(v);

	return this->mipmap[level].get_texel(tu,tv);
}

Color Texture::lerp(double w, Color c0, Color c1) {
	
	float r, g, b;
	r = w * c1[0] + (1.f - w) * c0[0];
	g = w * c1[1] + (1.f - w) * c0[1];
	b = w * c1[2] + (1.f - w) * c0[2];

	return Color(r, g, b);
}

// Returns the bilinear sample given a particular level and set of uv coords
Color Texture::sample_bilinear(Vector2D uv, int level) {
  // Optional helper function for Parts 5 and 6
  // Feel free to ignore or create your own

	size_t w, h;
	w = this->mipmap[level].width;
	h = this->mipmap[level].height;
	double u, v;
	u = uv[0] * w;
	v = uv[1] * h;

	// (uc, uv) is the point in the center of the four pixels for bilinear filtering
	double uc, vc;
	uc = (abs(u - ceil(u)) < abs(u - floor(u))) ? ceil(u) : floor(u);
	vc = (abs(v - ceil(v)) < abs(v - floor(v))) ? ceil(v) : floor(v);

	Vector2D u00_uv, u01_uv, u10_uv, u11_uv;

	u00_uv = Vector2D(uc - 0.5, vc + 0.5);
	double s = u - u00_uv[0], t = -(v - u00_uv[1]);

	u00_uv = Vector2D(uc - 1, vc);
	u01_uv = Vector2D(uc - 1, vc - 1);
	u10_uv = Vector2D(uc, vc);
	u11_uv = Vector2D(uc, vc - 1);

	Color c00, c01, c10, c11, c0, c1, cf;
	c00 = this->mipmap[level].get_texel(u00_uv[0], u00_uv[1]);
	c10 = this->mipmap[level].get_texel(u10_uv[0], u10_uv[1]);
	c0 = lerp(s, c00, c10);
	c01 = this->mipmap[level].get_texel(u01_uv[0], u01_uv[1]);
	c11 = this->mipmap[level].get_texel(u11_uv[0], u11_uv[1]);
	c1 = lerp(s, c01, c11);

	return lerp(t, c0, c1);
}

// does bilinear sampling if int level is a whole number, otherwise does trilinear sampling
Color Texture::sample_linear(Vector2D uv, int level) {
	if (abs(floor(level)) == level) // if level is a whole number, do bilinear sampling
		return sample_bilinear(uv, level);
	else { // trilinear sampling
		Color c0, c1, cf;
		c0 = sample_bilinear(uv, floor(level));
		c1 = sample_bilinear(uv, ceil(level));

		double w = level - floor(level);
		
		return lerp(w, c0, c1);
	}
}

/****************************************************************************/

inline void uint8_to_float(float dst[3], unsigned char *src) {
  uint8_t *src_uint8 = (uint8_t *)src;
  dst[0] = src_uint8[0] / 255.f;
  dst[1] = src_uint8[1] / 255.f;
  dst[2] = src_uint8[2] / 255.f;
}

inline void float_to_uint8(unsigned char *dst, float src[3]) {
  uint8_t *dst_uint8 = (uint8_t *)dst;
  dst_uint8[0] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[0])));
  dst_uint8[1] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[1])));
  dst_uint8[2] = (uint8_t)(255.f * max(0.0f, min(1.0f, src[2])));
}

void Texture::generate_mips(int startLevel) {

  // make sure there's a valid texture
  if (startLevel >= mipmap.size()) {
    std::cerr << "Invalid start level";
  }

  // allocate sublevels
  int baseWidth = mipmap[startLevel].width;
  int baseHeight = mipmap[startLevel].height;
  int numSubLevels = (int)(log2f((float)max(baseWidth, baseHeight)));

  numSubLevels = min(numSubLevels, kMaxMipLevels - startLevel - 1);
  mipmap.resize(startLevel + numSubLevels + 1);

  int width = baseWidth;
  int height = baseHeight;
  for (int i = 1; i <= numSubLevels; i++) {

    MipLevel &level = mipmap[startLevel + i];

    // handle odd size texture by rounding down
    width = max(1, width / 2);
    //assert (width > 0);
    height = max(1, height / 2);
    //assert (height > 0);

    level.width = width;
    level.height = height;
    level.texels = vector<unsigned char>(3 * width * height);
  }

  // create mips
  int subLevels = numSubLevels - (startLevel + 1);
  for (int mipLevel = startLevel + 1; mipLevel < startLevel + subLevels + 1;
       mipLevel++) {

    MipLevel &prevLevel = mipmap[mipLevel - 1];
    MipLevel &currLevel = mipmap[mipLevel];

    int prevLevelPitch = prevLevel.width * 3; // 32 bit RGB
    int currLevelPitch = currLevel.width * 3; // 32 bit RGB

    unsigned char *prevLevelMem;
    unsigned char *currLevelMem;

    currLevelMem = (unsigned char *)&currLevel.texels[0];
    prevLevelMem = (unsigned char *)&prevLevel.texels[0];

    float wDecimal, wNorm, wWeight[3];
    int wSupport;
    float hDecimal, hNorm, hWeight[3];
    int hSupport;

    float result[3];
    float input[3];

    // conditional differentiates no rounding case from round down case
    if (prevLevel.width & 1) {
      wSupport = 3;
      wDecimal = 1.0f / (float)currLevel.width;
    } else {
      wSupport = 2;
      wDecimal = 0.0f;
    }

    // conditional differentiates no rounding case from round down case
    if (prevLevel.height & 1) {
      hSupport = 3;
      hDecimal = 1.0f / (float)currLevel.height;
    } else {
      hSupport = 2;
      hDecimal = 0.0f;
    }

    wNorm = 1.0f / (2.0f + wDecimal);
    hNorm = 1.0f / (2.0f + hDecimal);

    // case 1: reduction only in horizontal size (vertical size is 1)
    if (currLevel.height == prevLevel.height) {
      //assert (currLevel.height == 1);

      for (int i = 0; i < currLevel.width; i++) {
        wWeight[0] = wNorm * (1.0f - wDecimal * i);
        wWeight[1] = wNorm * 1.0f;
        wWeight[2] = wNorm * wDecimal * (i + 1);

        result[0] = result[1] = result[2] = 0.0f;

        for (int ii = 0; ii < wSupport; ii++) {
          uint8_to_float(input, prevLevelMem + 3 * (2 * i + ii));
          result[0] += wWeight[ii] * input[0];
          result[1] += wWeight[ii] * input[1];
          result[2] += wWeight[ii] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (3 * i), result);
      }

      // case 2: reduction only in vertical size (horizontal size is 1)
    } else if (currLevel.width == prevLevel.width) {
      //assert (currLevel.width == 1);

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        result[0] = result[1] = result[2] = 0.0f;
        for (int jj = 0; jj < hSupport; jj++) {
          uint8_to_float(input, prevLevelMem + prevLevelPitch * (2 * j + jj));
          result[0] += hWeight[jj] * input[0];
          result[1] += hWeight[jj] * input[1];
          result[2] += hWeight[jj] * input[2];
        }

        // convert back to format of the texture
        float_to_uint8(currLevelMem + (currLevelPitch * j), result);
      }

      // case 3: reduction in both horizontal and vertical size
    } else {

      for (int j = 0; j < currLevel.height; j++) {
        hWeight[0] = hNorm * (1.0f - hDecimal * j);
        hWeight[1] = hNorm;
        hWeight[2] = hNorm * hDecimal * (j + 1);

        for (int i = 0; i < currLevel.width; i++) {
          wWeight[0] = wNorm * (1.0f - wDecimal * i);
          wWeight[1] = wNorm * 1.0f;
          wWeight[2] = wNorm * wDecimal * (i + 1);

          result[0] = result[1] = result[2] = 0.0f;

          // convolve source image with a trapezoidal filter.
          // in the case of no rounding this is just a box filter of width 2.
          // in the general case, the support region is 3x3.
          for (int jj = 0; jj < hSupport; jj++)
            for (int ii = 0; ii < wSupport; ii++) {
              float weight = hWeight[jj] * wWeight[ii];
              uint8_to_float(input, prevLevelMem +
                                        prevLevelPitch * (2 * j + jj) +
                                        3 * (2 * i + ii));
              result[0] += weight * input[0];
              result[1] += weight * input[1];
              result[2] += weight * input[2];
            }

          // convert back to format of the texture
          float_to_uint8(currLevelMem + currLevelPitch * j + 3 * i, result);
        }
      }
    }
  }
}

}
