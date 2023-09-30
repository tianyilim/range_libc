#ifndef _RANGE_UTILS_H_INCLUDED_
#define _RANGE_UTILS_H_INCLUDED_

#include <algorithm>
#include <climits>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

// TODO add documentation

namespace utils {

// For below function
[[maybe_unused]] static unsigned long x = 123456789, y = 362436069, z = 521288629;

/// @brief  Pseudo Random number generator
unsigned long xorshf96(void);

/// @brief Returns a grayscale value by mixing 3 rgb values
/// @param[in] r color intensity value
/// @param[in] g color intensity value
/// @param[in] b color intensity value
/// @return mixed pixel value
float rgb2gray(float r, float g, float b);

/// @brief Random value within range
int randrange(int min, int max);

/// @brief Returns the 4 (or 8, if use_corners) grid cells that are neighbors of the input grid
/// coordinate (x, y).
/// @param[in] x x coordinate of grid cell
/// @param[in] y y coordinate of grid cell
/// @param[in] use_corners whether to include the diagonal neighbors
/// @return vector of neighbors
std::vector<std::pair<int, int>> outline(int x, int y, bool use_corners);

template <class key_T>
class KeyMaker {
   public:
    KeyMaker() {}
    KeyMaker(int width, int height, int theta_discretization)
    {
        y_shift = (int)std::ceil(std::log2(theta_discretization));
        x_shift = (int)std::ceil(std::log2(height)) + y_shift;
        int bitness = (int)std::ceil(std::log2(width)) + x_shift;

        if (bitness > std::log2(std::numeric_limits<key_T>::max())) {
            std::cerr << "Key bitness too large for integer packing scheme. Check "
                         "your KeyMaker template type."
                      << std::endl;
        }

        // make bit masks for unpacking the various values
        t_mask = std::pow(2, y_shift) - 1;
        y_mask = std::pow(2, x_shift) - 1 - t_mask;
        x_mask = std::pow(2, bitness) - 1 - y_mask - t_mask;
    };
    ~KeyMaker(){};
    key_T make_key(int x, int y, int t)
    {
        return ((key_T)x << x_shift) + ((key_T)y << y_shift) + (key_T)t;
    }
    std::tuple<int, int, int> unpack_key(key_T k)
    {
        return std::make_tuple((int)((k & x_mask) >> x_shift), (int)((k & y_mask) >> y_shift),
                               k & t_mask);
    }

   private:
    int y_shift;
    int x_shift;
    key_T x_mask;
    key_T y_mask;
    key_T t_mask;
};

bool has(std::string substring, std::string str);

bool has(std::string val, std::vector<std::string> vstr);

std::vector<std::string> split(std::string in, char delim);

double norminv(double q);

template <typename T, typename U>
struct is_same {
    static const bool value = false;
};

template <typename T>
struct is_same<T, T> {
    static const bool value = true;
};

template <typename T, typename U>
bool eqTypes()
{
    return is_same<T, U>::value;
}

// http://stackoverflow.com/questions/311703/algorithm-for-sampling-without-replacement
// Here's some code for sampling without replacement based on Algorithm 3.4.2S
// of Knuth's book Seminumeric Algorithms.
class NonReplacementSampler {
   public:
    NonReplacementSampler();
    ~NonReplacementSampler() {}

    void sample(int populationSize, int sampleSize, std::vector<int> &samples);

    std::uniform_real_distribution<double> rand;
    std::default_random_engine generator;
};

class FastRand {
   public:
    FastRand();
    FastRand(int n);

    float rand();

    void populate();

    int i = 0;
    int cache_size;
    float repopulate_threshold;
    // std::vector<float> cache;
    float cache[10000];
};

void serialize(std::vector<bool> &vals, std::stringstream *ss);

void serialize(std::vector<float> &vals, std::stringstream *ss);

std::string serialize(std::vector<float> &vals);
}  // namespace utils

#endif /* _RANGE_UTILS_H_INCLUDED_ */
