#pragma once

#include <algorithm>
#ifdef WINDOWS
#include <numbers>
#else
#define INT_MAX       2147483647
#endif
#include <complex>
#include <cstdint>
#include <cmath>

constexpr float DEG_TO_RAD = (float)3.14159265358979323846 / 180.0f;
constexpr float EPSILON = 0.000001f;
constexpr float SQRT3 = 1.7320508075688772935274463415059;

inline bool areEqual(float a, float b, float epsilon = EPSILON)
{
    return (b > (a - epsilon)) && (b < (a + epsilon));
}

inline float squared(float& v)
{
    return v * v;
}

inline float getMagnitude(float a, float b)
{
    return sqrtf(a * a + b * b);
}

inline double getMagnitude(double a, double b)
{
    return sqrt(a * a + b * b);
}

template <typename T>
inline T getMagnitude(std::complex<T> v)
{
    return getMagnitude(v.real(), v.imag());
}

inline float getPhase(float x, float y)
{
    return atan2f(y, x);
}

inline double getPhase(double x, double y)
{
    return atan2(y, x);
}

inline float getDistance(float x1, float y1, float x2, float y2)
{
    return sqrtf(powf(x2 - x1, 2) + powf(y2 - y1, 2));
}


inline float lerp(float a, float b, float t)
{
    return a + ((b - a) * t);
}

inline float mapRange(float value, float inpStart, float inpEnd, float outStart, float outEnd)
{
    return outStart + ((outEnd - outStart) / (inpEnd - inpStart)) * (value - inpStart);
}

inline float mapRangeClamp(float value, float inpStart, float inpEnd, float outStart, float outEnd)
{
  //  float t = std::clamp((value - inpStart) / (inpEnd - inpStart), 0.0f, 1.0f);
    float t = (value - inpStart) / (inpEnd - inpStart);
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    return outStart + t * (outEnd - outStart);
}

inline uint32_t upperPowerOf2(uint32_t v)
{
    return (uint32_t)floor(pow(2, ceil(log(v) / log(2))));
}

inline int shiftIndex(int i, int shift, int size)
{
    return (i + shift) % size;
}

inline uint32_t i32ToU32(int v)
{
    return std::max(v, 0);
}

inline int u32ToI32(uint32_t v)
{
    return std::min(v, (uint32_t)INT_MAX);
}

enum class GrayscaleType
{
    None,
    Luminance,
    Average,
    Sum,
    Maximum,
    Minimum,
    Magnitude,
    MagOverSqrt3,
    Red,
    Green,
    Blue,
    Alpha
};
constexpr uint32_t GrayscaleType_EnumSize = 12;

constexpr GrayscaleType CONV_THRESHOLD_GRAYSCALE_TYPE = GrayscaleType::Average;

inline float rgbaToGrayscale(float* rgba, GrayscaleType type)
{
    switch (type)
    {
    case GrayscaleType::Luminance:
        return (rgba[0] * 0.2126f) + (rgba[1] * 0.7152f) + (rgba[2] * 0.0722f);
        break;
    case GrayscaleType::Average:
        return (rgba[0] + rgba[1] + rgba[2]) / 3.0f;
        break;
    case GrayscaleType::Sum:
        return rgba[0] + rgba[1] + rgba[2];
        break;
    case GrayscaleType::Maximum:
        return std::max(std::max(rgba[0], rgba[1]), rgba[2]);
        break;
    case GrayscaleType::Minimum:
        return std::min(std::min(rgba[0], rgba[1]), rgba[2]);
        break;
    case GrayscaleType::Magnitude:
        return sqrtf((rgba[0] * rgba[0]) + (rgba[1] * rgba[1]) + (rgba[2] * rgba[2]));
        break;
    case GrayscaleType::MagOverSqrt3:
        return sqrtf((rgba[0] * rgba[0]) + (rgba[1] * rgba[1]) + (rgba[2] * rgba[2])) / SQRT3;
        break;
    case GrayscaleType::Red:
        return rgba[0];
        break;
    case GrayscaleType::Green:
        return rgba[1];
        break;
    case GrayscaleType::Blue:
        return rgba[2];
        break;
    case GrayscaleType::Alpha:
        return rgba[3];
        break;
    default:
        break;
    }

    return rgba[0];
}

inline float rgbToGrayscale(float* rgb, GrayscaleType type)
{
    switch (type)
    {
    case GrayscaleType::Alpha:
        return 1.0f;
        break;
    default:
        return rgbaToGrayscale(rgb, type);
        break;
    }
}

inline float transformKnee(float v)
{
    if (v <= 0.0f)
        return 0.0f;
    return powf(2.0f, v) - 1.0f;
}

inline float softThreshold(float v, float threshold, float transKnee)
{
    return transKnee > 0 ? fminf(1.0f, (v - threshold) / transKnee) : 1.0f;
}

inline bool checkBounds(int x, int y, int w, int h)
{
    if (x < 0)
        return false;

    if (x >= w)
        return false;

    if (y < 0)
        return false;

    if (y >= h)
        return false;

    return true;
}

inline void blendAddRGBA(float* colorA, uint32_t indexA, const float* colorB, uint32_t indexB, float t)
{
    colorA[indexA + 0] += colorB[indexB + 0] * t;
    colorA[indexA + 1] += colorB[indexB + 1] * t;
    colorA[indexA + 2] += colorB[indexB + 2] * t;
    colorA[indexA + 3] += colorB[indexB + 3] * t;
}

inline void blendAddRGB(float* colorA, uint32_t indexA, const float* colorB, uint32_t indexB, float t)
{
    colorA[indexA + 0] += colorB[indexB + 0] * t;
    colorA[indexA + 1] += colorB[indexB + 1] * t;
    colorA[indexA + 2] += colorB[indexB + 2] * t;
}

inline float getExposureMul(float v)
{
    return powf(2.0f, v);
}

inline double getExposureMul(double v)
{
    return pow(2.0, v);
}

float applyContrast(float v, float contrast);

void rotatePoint(float x, float y, float pivotX, float pivotY, float angle, float& outX, float& outY);
void rotatePointInPlace(float& x, float& y, float pivotX, float pivotY, float angle);

void calcFftConvPadding(
    bool powerOfTwo,
    bool square,
    uint32_t inputWidth,
    uint32_t inputHeight,
    uint32_t kernelWidth,
    uint32_t kernelHeight,
    float kernelCenterX,
    float kernelCenterY,
    uint32_t& outPaddedWidth,
    uint32_t& outPaddedHeight);

void calcDispScale(uint32_t index, uint32_t steps, float amount, float edgeOffset, float& outScale, float& outAreaMul);

float srgbToLinear_DEPRECATED(float x);
float linearToSrgb_DEPRECATED(float x);
