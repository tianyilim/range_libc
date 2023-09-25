/*
* Copyright 2017 Corey H. Walsh (corey.walsh11@gmail.com)

* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at

*     http://www.apache.org/licenses/LICENSE-2.0

* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/

#ifndef RANGE_LIB_H
#define RANGE_LIB_H

#include <time.h>
#include <unistd.h>

#include <algorithm>  // std::min
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>  // std::setw
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "RangeUtils.h"
#include "distance_transform.h"
#include "lodepng/lodepng.h"
// #define NDEBUG
#include <cassert>
#include <tuple>

#define _TRACK_LUT_SIZE 0
#define _TRACK_COLLISION_INDEXES 0

#define _EPSILON 0.00001
#define M_2PI 2 * M_PI
#define _BINARY_SEARCH_THRESHOLD 64  // if there are more than this number of elements in the lut bin, use binary search

// fast optimized version
#define _USE_CACHED_TRIG 0
#define _USE_CACHED_CONSTANTS 1
#define _USE_FAST_ROUND 0
#define _NO_INLINE 0
#define _USE_LRU_CACHE 0
#define _LRU_CACHE_SIZE 1000000

// not implemented yet -> use 16 bit integers to store zero points
#define _CDDT_SHORT_DATATYPE 1
#define _GIANT_LUT_SHORT_DATATYPE 1

#define SENSOR_MODEL_HELPERS 1

#if _USE_LRU_CACHE
#include "include/lru_cache.h"
#endif

// these defines are for yaml/JSON serialization
#define T1 "  "
#define T2 T1 T1
#define T3 T1 T1 T1
#define T4 T2 T2

#define J1 "  "
#define J2 J1 J1
#define J3 J1 J1 J1
#define J4 J2 J2

#if USE_CUDA == 1
#ifndef CHUNK_SIZE
#define CHUNK_SIZE 262144
#define CHUNK_THREADS 256
#endif
#include "CudaRangeLib.h"
#else
#define USE_CUDA 0
#endif

namespace ranges {

/// @brief Occupancy Grid
class OMap {
   public:
    typedef std::vector<std::vector<bool>> Grid_t;
    typedef std::vector<std::vector<float>> RawGrid_t;

    /// @brief Constructor from sizes
    /// @param[in] w width of occupancy grid
    /// @param[in] h height of occupancy grid
    OMap(unsigned w = 0, unsigned h = 0) : _hasError{false}, _width{w}, _height{h}, _filename{""}
    {
        // Initializes a blank occupancy grid map.
        for (int i = 0; i < w; ++i) {
            std::vector<bool> y_axis;
            for (int q = 0; q < h; ++q)
                y_axis.push_back(false);
            _OccupancyGrid.push_back(y_axis);
        }
    }

    /// @brief Constructor from an image file and belief threshold
    /// @param[in] filename path to the image file
    /// @param[in] threshold
    OMap(std::string filename, float threshold = 128) : _hasError(false), _filename(filename)
    {
        unsigned error;
        unsigned char *image;

        error = lodepng_decode32_file(&image, &_width, &_height, filename.c_str());
        if (error) {
            std::cerr << "ERROR " << error << lodepng_error_text(error) << std::endl;
            _hasError = true;
            return;
        }

        for (int i = 0; i < _width; ++i) {
            std::vector<bool> y_axis;
            for (int q = 0; q < _height; ++q)
                y_axis.push_back(false);
            _OccupancyGrid.push_back(y_axis);
        }

        for (int i = 0; i < _width; ++i) {
            std::vector<float> y_axis;
            for (int q = 0; q < _height; ++q)
                y_axis.push_back(0);
            _SignedDistFunc.push_back(y_axis);
        }

        for (int y = 0; y < _height; ++y) {
            for (int x = 0; x < _width; ++x) {
                unsigned idx = 4 * y * _width + 4 * x;
                int r = image[idx + 2];
                int g = image[idx + 1];
                int b = image[idx + 0];
                int gray = (int)utils::rgb2gray(r, g, b);
                if (gray < threshold)
                    _OccupancyGrid[x][y] = true;
                _SignedDistFunc[x][y] = gray;
            }
        }
    }

    /// @brief Return the occupancy at x, y
    inline bool isOccupied(int x, int y)
    {
        if (x < 0 || x >= _width || y < 0 || y >= _height)
            return false;
        else
            return _OccupancyGrid[x][y];
    }

    // TODO this should be moved
    /// @brief Returns the signed distance value at x, y
    inline float signedDistanceValue(int x, int y)
    {
        if (x < 0 || x >= _width || y < 0 || y >= _height)
            return -1;  // this should throw
        else
            return _SignedDistFunc[x][y];
    }

    /// @brief Save map to an image file
    /// @param[in] filename
    /// @return error code
    bool save(std::string &filename)
    {
        std::vector<unsigned char> png;
        lodepng::State state;  // optionally customize this one
        char image[_width * _height * 4];

        for (int y = 0; y < _height; ++y) {
            for (int x = 0; x < _width; ++x) {
                unsigned idx = 4 * y * _width + 4 * x;

                image[idx + 2] = (char)255;
                image[idx + 1] = (char)255;
                image[idx + 0] = (char)255;
                image[idx + 3] = (char)255;

                if (_OccupancyGrid[x][y]) {
                    image[idx + 0] = 0;
                    image[idx + 1] = 0;
                    image[idx + 2] = 0;
                }
            }
        }
        unsigned error = lodepng::encode(png, reinterpret_cast<const unsigned char *>(image), _width, _height, state);
        if (!error)
            lodepng::save_file(png, filename);
        else
            std::cerr << "Image encoder error " << error << ": " << lodepng_error_text(error) << std::endl;

        return error;
    }

    /// @brief Creates an edge map from the existing OMap
    /// @param[in] count_corners
    /// @return
    OMap makeEdgeMap(bool count_corners)
    {
        OMap edge_map = OMap(_width, _height);
        for (int x = 0; x < _width; ++x) {
            for (int y = 0; y < _height; ++y) {
                if (!isOccupied(x, y))
                    continue;

                std::vector<std::pair<int, int>> outline = utils::outline(x, y, count_corners);
                for (int i = 0; i < outline.size(); ++i) {
                    int cx, cy;

                    std::tie(cx, cy) = outline[i];
                    if (0 <= cx && 0 <= cy && cx < _width && cy < _height && !isOccupied(cx, cy)) {
                        edge_map.grid()[x][y] = true;
                        break;
                    }
                }
            }
        }
        return edge_map;
    }

    /// @brief returns memory usage in bytes
    inline int memory()
    {
        return sizeof(bool) * _width * _height;
    }

    /// @brief error accessor
    inline bool error()
    {
        return _hasError;
    }

    ///@brief modifier view into _OccupancyGrid
    Grid_t &grid() { return _OccupancyGrid; }

    /// @brief Const view into _SignedDistFunc
    const RawGrid_t &rawGrid() const { return _SignedDistFunc; }
    ///@brief modifier view into _SignedDistFunc
    RawGrid_t &rawGrid() { return _SignedDistFunc; }

    const unsigned width() const { return _width; }
    void setWidth(unsigned w) { _width = w; }

    const unsigned height() const { return _height; }
    void setHeight(unsigned h) { _height = h; }

    /// @brief filename accessor
    const std::string &filename() const { return _filename; }

    // TODO these should be added to some sort of struct with getter/setter
    float _worldScale;     ///< Real-world values
    float _worldAngle;     ///< Real-world values
    float _worldOriginX;   ///< Real-world values
    float _worldOriginY;   ///< Real-world values
    float _worldSinAngle;  ///< Real-world values
    float _worldCosAngle;  ///< Real-world values

   protected:
    // Members
    bool _hasError;  ///< error flag
    // TODO this should be removed
    unsigned _width;        ///< x axis
    unsigned _height;       ///< y axis
    std::string _filename;  ///< filename of image

    Grid_t _OccupancyGrid;      ///< Grid of boolean occupied/not-occupied values
    RawGrid_t _SignedDistFunc;  ///< Signed Distance Function
    // TODO This should be moved to the DistanceTransform subclass as it is not used
};

class DistanceTransform : public OMap {
   public:
    /// @brief Forwards to OMap constructor
    /// @param[in] w Occupancy grid width
    /// @param[in] h Occupancy grid height
    DistanceTransform(unsigned w = 0, unsigned h = 0) : OMap(w, h) {}

    ///@brief computes the distance transform of a given OMap
    ///@param[in] map Input OMap
    DistanceTransform(OMap &map) : OMap(map)
    {
        std::vector<std::size_t> grid_size({_width, _height});
        dt::MMArray<float, 2> f(grid_size.data());
        dt::MMArray<std::size_t, 2> indices(grid_size.data());

        for (std::size_t i = 0; i < _width; ++i) {
            for (std::size_t j = 0; j < _height; ++j) {
                if (map.isOccupied(i, j)) {
                    f[i][j] = 0.0f;
                }
                else {
                    f[i][j] = std::numeric_limits<float>::max();
                }
            }
        }

        dt::DistanceTransform::distanceTransformL2(f, f, indices, false);

        // allocate space in the vectors
        _SignedDistFunc.clear();
        for (int x = 0; x < _width; ++x) {
            std::vector<float> y_axis;
            for (int y = 0; y < _height; ++y) {
                y_axis.push_back(f[x][y]);
            }
            _SignedDistFunc.push_back(y_axis);
        }
    }
};

/// @brief Parent class of all range methods
class RangeMethod {
   public:
    RangeMethod(OMap m, float mr) : map(m), max_range(mr){};

    virtual ~RangeMethod(){};

    virtual float calc_range(float x, float y, float heading) = 0;

    virtual std::pair<float, float> calc_range_pair(float x, float y, float heading) { return std::make_pair(-1, -1); }

    virtual OMap *getMap() { return &map; }

    virtual void report(){};

    float maxRange() const { return max_range; }

    float memory() const { return -1; }

    // wrapper function to call calc_range repeatedly with the given array of inputs
    // and store the result to the given outputs. Useful for avoiding cython function
    // call overhead by passing it a numpy array pointer. Indexing assumes a 3xn numpy array
    // for the inputs and a 1xn numpy array of the outputs
    void numpy_calc_range(float *ins, float *outs, int num_casts)
    {
        // cache these constants on the stack for efficiency
        float inv_world_scale = 1.0 / map._worldScale;
        float _worldScale = map._worldScale;
        float _worldAngle = map._worldAngle;
        float _worldOriginX = map._worldOriginX;
        float _worldOriginY = map._worldOriginY;
        float _worldSinAngle = map._worldSinAngle;
        float _worldCosAngle = map._worldCosAngle;

        float rotation_const = -1.0 * _worldAngle - 3.0 * M_PI / 2.0;

        // avoid allocation on every loop iteration
        float x_world;
        float y_world;
        float theta_world;
        float x;
        float y;
        float temp;
        float theta;

        for (int i = 0; i < num_casts; ++i) {
            x_world = ins[i * 3];
            y_world = ins[i * 3 + 1];
            theta_world = ins[i * 3 + 2];

            x = (x_world - _worldOriginX) * inv_world_scale;
            y = (y_world - _worldOriginY) * inv_world_scale;
            temp = x;
            x = _worldCosAngle * x - _worldSinAngle * y;
            y = _worldSinAngle * temp + _worldCosAngle * y;
            theta = -theta_world + rotation_const;

            outs[i] = calc_range(y, x, theta) * _worldScale;
        }
    }

    void numpy_calc_range_angles(float *ins, float *angles, float *outs, int num_particles, int num_angles)
    {
        // cache these constants on the stack for efficiency
        float inv_world_scale = 1.0 / map._worldScale;
        float _worldScale = map._worldScale;
        float _worldAngle = map._worldAngle;
        float _worldOriginX = map._worldOriginX;
        float _worldOriginY = map._worldOriginY;
        float _worldSinAngle = map._worldSinAngle;
        float _worldCosAngle = map._worldCosAngle;
        float rotation_const = -1.0 * _worldAngle - 3.0 * M_PI / 2.0;

        // avoid allocation on every loop iteration
        float x_world;
        float y_world;
        float theta_world;
        float x;
        float y;
        float temp;
        float theta;

        for (int i = 0; i < num_particles; ++i) {
            x_world = ins[i * 3];
            y_world = ins[i * 3 + 1];
            theta_world = ins[i * 3 + 2];
            theta = -theta_world + rotation_const;

            x = (x_world - _worldOriginX) * inv_world_scale;
            y = (y_world - _worldOriginY) * inv_world_scale;
            temp = x;
            x = _worldCosAngle * x - _worldSinAngle * y;
            y = _worldSinAngle * temp + _worldCosAngle * y;

            for (int a = 0; a < num_angles; ++a)
                outs[i * num_angles + a] = calc_range(y, x, theta - angles[a]) * _worldScale;
        }
    }

    void set_sensor_model(double *table, int table_width)
    {
        // convert the sensor model from a numpy array to a vector array
        for (int i = 0; i < table_width; ++i) {
            std::vector<double> table_row;
            for (int j = 0; j < table_width; ++j)
                table_row.push_back(table[table_width * i + j]);
            sensor_model.push_back(table_row);
        }
    }

    void eval_sensor_model(float *obs, float *ranges, double *outs, int rays_per_particle, int particles)
    {
        float inv_world_scale = 1.0 / map._worldScale;
        // do no allocations in the main loop
        double weight;
        float r, d;

        for (int i = 0; i < particles; ++i) {
            weight = 1.0;
            for (int j = 0; j < rays_per_particle; ++j) {
                r = obs[j] * inv_world_scale;
                r = std::min<float>(std::max<float>(r, 0.0), (float)sensor_model.size() - 1.0);
                d = ranges[i * rays_per_particle + j] * inv_world_scale;
                d = std::min<float>(std::max<float>(d, 0.0), (float)sensor_model.size() - 1.0);
                weight *= sensor_model[(int)r][(int)d];
            }
            outs[i] = weight;
        }
    }

    // calc range for each pose, adding every angle, evaluating the sensor model
    void calc_range_repeat_angles_eval_sensor_model(float *ins, float *angles, float *obs, double *weights, int num_particles, int num_angles)
    {
        // cache these constants on the stack for efficiency
        float inv_world_scale = 1.0 / map._worldScale;
        float _worldScale = map._worldScale;
        float _worldAngle = map._worldAngle;
        float _worldOriginX = map._worldOriginX;
        float _worldOriginY = map._worldOriginY;
        float _worldSinAngle = map._worldSinAngle;
        float _worldCosAngle = map._worldCosAngle;
        float rotation_const = -1.0 * _worldAngle - 3.0 * M_PI / 2.0;

        // avoid allocation on every loop iteration
        float x_world;
        float y_world;
        float theta_world;
        float x;
        float y;
        float temp;
        float theta;

        // do no allocations in the main loop
        double weight;
        float r;
        float d;
        int i;
        int a;

        for (i = 0; i < num_particles; ++i) {
            x_world = ins[i * 3];
            y_world = ins[i * 3 + 1];
            theta_world = ins[i * 3 + 2];
            theta = -theta_world + rotation_const;

            x = (x_world - _worldOriginX) * inv_world_scale;
            y = (y_world - _worldOriginY) * inv_world_scale;
            temp = x;
            x = _worldCosAngle * x - _worldSinAngle * y;
            y = _worldSinAngle * temp + _worldCosAngle * y;

            weight = 1.0;
            for (a = 0; a < num_angles; ++a) {
                d = calc_range(y, x, theta - angles[a]);
                d = std::min<float>(std::max<float>(d, 0.0), (float)sensor_model.size() - 1.0);

                r = obs[a] * inv_world_scale;
                r = std::min<float>(std::max<float>(r, 0.0), (float)sensor_model.size() - 1.0);
                weight *= sensor_model[(int)r][(int)d];
            }
            weights[i] = weight;
        }
    }

    // this is to compute a lidar sensor model using radial (calc_range_pair) optimizations
    // this is only exact for a certain set of downsample amounts
    void calc_range_many_radial_optimized(float *ins, float *outs, int num_particles, int num_rays, float min_angle, float max_angle)
    {
        // cache these constants on the stack for efficiency
        float inv_world_scale = 1.0 / map._worldScale;
        float _worldScale = map._worldScale;
        float _worldAngle = map._worldAngle;
        float _worldOriginX = map._worldOriginX;
        float _worldOriginY = map._worldOriginY;
        float _worldSinAngle = map._worldSinAngle;
        float _worldCosAngle = map._worldCosAngle;
        float rotation_const = -1.0 * _worldAngle - 3.0 * M_PI / 2.0;

        // avoid allocation on every loop iteration
        float x_world, y_world, theta_world, x, y, temp, theta = 0.0;

        // do no allocations in the main loop
        int i, a = 0;

        float step = (max_angle - min_angle) / (num_rays - 1);
        float angle = min_angle;

        int max_pairwise_index = (float)num_rays / 3.0;
        float index_offset_float = (num_rays - 1.0) * M_PI / (max_angle - min_angle);

        // TODO: check if this index_offset_float is not very close to an integer
        // in which case throw a warning that this downsample factor is poorly compatible
        // with this radial optimization
        int index_offset = roundf(index_offset_float);
        float r, r_inv = 0.0;

        for (i = 0; i < num_particles; ++i) {
            x_world = ins[i * 3];
            y_world = ins[i * 3 + 1];
            theta_world = ins[i * 3 + 2];
            theta = -theta_world + rotation_const;

            x = (x_world - _worldOriginX) * inv_world_scale;
            y = (y_world - _worldOriginY) * inv_world_scale;
            temp = x;
            x = _worldCosAngle * x - _worldSinAngle * y;
            y = _worldSinAngle * temp + _worldCosAngle * y;

            angle = min_angle;
            for (a = 0; a <= max_pairwise_index; ++a) {
                std::tie(r, r_inv) = calc_range_pair(y, x, theta - angle);
                outs[i * num_rays + a] = r * _worldScale;
                outs[i * num_rays + a + index_offset] = r_inv * _worldScale;
                angle += step;
            }

            for (a = max_pairwise_index + 1; a < index_offset; ++a) {
                outs[i * num_rays + a] = calc_range(y, x, theta - angle) * _worldScale;
                angle += step;
            }
        }
        return;
    }

   protected:
    OMap map;  // TODO perhaps dont store an OMap and a DistanceTransform simultaneously. Or just store a DistanceTransform?
    float max_range;

    std::vector<std::vector<double>> sensor_model;
};

class BresenhamsLine : public RangeMethod {
   public:
    BresenhamsLine(OMap m, float mr) : RangeMethod(m, mr){};

    float calc_range(float x, float y, float heading)
    {
        // first check if the cell underneath the query point is occupied, if so return
        if (map.isOccupied((int)x, (int)y)) {
            return 0.0;
        }

        /*
         this defines the coordinate system such that
            ------> +x
            |
            |
            \/
            +y
          0* heading lies along the x axis, positive heading rotates towards the positive y axis
        */
        float x0 = y;
        float y0 = x;
        float x1 = y + max_range * sinf(heading);
        float y1 = x + max_range * cosf(heading);

        bool steep = false;
        if (std::abs(y1 - y0) > std::abs(x1 - x0))
            steep = true;

        if (steep) {
            float tmp = x0;
            x0 = y0;
            y0 = tmp;
            tmp = x1;
            x1 = y1;
            y1 = tmp;
        }

        float deltax = std::abs(x1 - x0);
        float deltay = std::abs(y1 - y0);

        float error = 0;
        float deltaerr = deltay;
        float _x = x0;
        float _y = y0;

        int xstep = -1;
        if (x0 < x1)
            xstep = 1;

        int ystep = -1;
        if (y0 < y1)
            ystep = 1;

        unsigned width = map.width();
        unsigned height = map.height();

        while ((int)_x != (int)(x1 + xstep)) {
            _x += xstep;
            error += deltaerr;

            if (error * 2.00 >= deltax) {
                _y += ystep;
                error -= deltax;
            }

            if (!steep) {
                if (0 <= _y && _y < width && 0 <= _x && _x < height && map.isOccupied(_y, _x)) {
                    float xd = _x - x0;
                    float yd = _y - y0;
                    return sqrtf(xd * xd + yd * yd);
                }
            }
            else {
                if (0 <= _x && _x < width && 0 <= _y && _y < height && map.isOccupied(_x, _y)) {
                    float xd = _x - x0;
                    float yd = _y - y0;
                    return sqrtf(xd * xd + yd * yd);
                }
            }
        }
        return max_range;
    }

    int memory() { return map.memory(); }
};

class RayMarchingGPU : public RangeMethod {
   public:
    RayMarchingGPU(OMap m, float mr) : RangeMethod(m, mr)
    {
        distImage = new DistanceTransform(m);
#if USE_CUDA == 1
        rmc = new RayMarchingCUDA(distImage->grid(), distImage->width, distImage->height, max_range);

        rmc->set_conversion_params(m._worldScale, m._worldAngle, m._worldOriginX, m._worldOriginY,
                                   m._worldSinAngle, m._worldCosAngle);

#else
        throw std::string("Must compile with -DWITH_CUDA=ON to use this class.");
#endif
    }
    ~RayMarchingGPU()
    {
        delete distImage;
#if USE_CUDA == 1
        delete rmc;
#else
        throw std::string("Must compile with -DWITH_CUDA=ON to use this class.");
#endif
    };

    float calc_range(float x, float y, float heading)
    {
#if USE_CUDA == 1
        std::cout << "Do not call calc_range on RayMarchingGPU, requires batched queries" << std::endl;
        return -1.0;
#else
        throw std::string("Must compile with -DWITH_CUDA=ON to use this class.");
#endif
    }

    void maybe_warn(int num_casts)
    {
#if USE_CUDA == 1
        if (!(num_casts % CHUNK_SIZE == 0) && !already_warned) {
            std::cout << "\nFor better performance, call calc_range_many with some multiple of " << CHUNK_SIZE << " queries. ";
            std::cout << "You can change the chunk size with -DCHUNK_SIZE=[integer].\n"
                      << std::endl;
            already_warned = true;
        }
#endif
    }

    void calc_range_many(float *ins, float *outs, int num_casts)
    {
#if USE_CUDA == 1
        maybe_warn(num_casts);
        int iters = std::ceil((float)num_casts / (float)CHUNK_SIZE);
        for (int i = 0; i < iters; ++i) {
            int num_in_chunk = CHUNK_SIZE;
            if (i == iters - 1)
                num_in_chunk = num_casts - i * CHUNK_SIZE;
            rmc->calc_range_many(&ins[i * CHUNK_SIZE * 3], &outs[i * CHUNK_SIZE], num_in_chunk);
        }
#else
        throw std::string("Must compile with -DWITH_CUDA=ON to use this class.");
#endif
    }

    // wrapper function to call calc_range repeatedly with the given array of inputs
    // and store the result to the given outputs. Useful for avoiding cython function
    // call overhead by passing it a numpy array pointer. Indexing assumes a 3xn numpy array
    // for the inputs and a 1xn numpy array of the outputs
    void numpy_calc_range(float *ins, float *outs, int num_casts)
    {
#if USE_CUDA == 1
        maybe_warn(num_casts);
        int iters = std::ceil((float)num_casts / (float)CHUNK_SIZE);
        for (int i = 0; i < iters; ++i) {
            int num_in_chunk = CHUNK_SIZE;
            if (i == iters - 1)
                num_in_chunk = num_casts - i * CHUNK_SIZE;
            rmc->numpy_calc_range(&ins[i * CHUNK_SIZE * 3], &outs[i * CHUNK_SIZE], num_in_chunk);
        }
#else
        throw std::string("Must compile with -DWITH_CUDA=ON to use this class.");
#endif
    }

    void numpy_calc_range_angles(float *ins, float *angles, float *outs, int num_particles, int num_angles)
    {
#if USE_CUDA == 1

        int particles_per_iter = std::ceil((float)CHUNK_SIZE / (float)num_angles);
        int iters = std::ceil((float)num_particles / (float)particles_per_iter);
        // must allways do the correct number of angles, can only split on the particles
        for (int i = 0; i < iters; ++i) {
            int num_in_chunk = particles_per_iter;
            if (i == iters - 1)
                num_in_chunk = num_particles - i * particles_per_iter;
            rmc->numpy_calc_range_angles(&ins[i * num_in_chunk * 3], angles, &outs[i * num_in_chunk * num_angles],
                                         num_in_chunk, num_angles);
        }
#else
        throw std::string("Must compile with -DWITH_CUDA=ON to use this class.");
#endif
    }

#if SENSOR_MODEL_HELPERS == 1
#if USE_CUDA == 1
    void set_sensor_model(double *table, int table_width)
    {
        // convert the sensor model from a numpy array to a vector array
        for (int i = 0; i < table_width; ++i) {
            std::vector<double> table_row;
            for (int j = 0; j < table_width; ++j)
                table_row.push_back(table[table_width * i + j]);
            sensor_model.push_back(table_row);
        }
        rmc->set_sensor_table(table, table_width);
    }
#endif

    // calc range for each pose, adding every angle, evaluating the sensor model
    void calc_range_repeat_angles_eval_sensor_model(float *ins, float *angles, float *obs, double *weights, int num_particles, int num_angles)
    {
#if USE_CUDA == 1

        int particles_per_iter = std::ceil((float)CHUNK_SIZE / (float)num_angles);
        int iters = std::ceil((float)num_particles / (float)particles_per_iter);
        // must allways do the correct number of angles, can only split on the particles
        for (int i = 0; i < iters; ++i) {
            int num_in_chunk = particles_per_iter;
            if (i == iters - 1)
                num_in_chunk = num_particles - i * particles_per_iter;
            rmc->calc_range_repeat_angles_eval_sensor_model(&ins[i * num_in_chunk * 3], angles, obs, &weights[i * num_in_chunk], num_in_chunk, num_angles);
        }
#else
        throw std::string("Must compile with -DWITH_CUDA=ON to use this class.");
#endif
    }
#endif

    int memory() { return distImage->memory(); }

   protected:
    DistanceTransform *distImage = 0;
#if USE_CUDA == 1
    RayMarchingCUDA *rmc = 0;
#endif
    bool already_warned = false;
};

class RayMarching : public RangeMethod {
   public:
    RayMarching(OMap m, float mr) : RangeMethod(m, mr)
    {
        _distImage = DistanceTransform(m);
    }

    float calc_range(const float x, const float y, const float heading)
    {
        float ray_direction_x = cosf(heading);
        float ray_direction_y = sinf(heading);

        int px, py;

        float t = 0.0;
        while (t < max_range) {
            px = x + ray_direction_x * t;
            py = y + ray_direction_y * t;

            if (px >= map.width() || px < 0 || py < 0 || py >= map.height()) {
                return max_range;
            }

            float d = _distImage.signedDistanceValue(px, py);

            if (d <= _distThreshold) {
                float xd = px - x;
                float yd = py - y;
                return sqrtf(xd * xd + yd * yd);
            }

            t += std::max<float>(d * _stepCoeff, 1.0);
        }

        return max_range;
    }

    int memory() { return _distImage.memory(); }

   protected:
    DistanceTransform _distImage;
    float _distThreshold = 0.0;
    float _stepCoeff = 0.999;
};

class CDDTCast : public RangeMethod {
   public:
    CDDTCast(OMap m, float mr, unsigned int td) : RangeMethod(m, mr), theta_discretization(td)
    {
#if _USE_CACHED_CONSTANTS
        theta_discretization_div_M_2PI = theta_discretization / M_2PI;
        M_2PI_div_theta_discretization = M_2PI / ((float)theta_discretization);
#endif

#if _USE_LRU_CACHE
        cache = cache::lru_cache<uint64_t, float>(_LRU_CACHE_SIZE, -1);
        key_maker = utils::KeyMaker<uint64_t>(m._width, m._height, theta_discretization);
#endif

        // determines the width of the projection of the map along each angle
        std::vector<int> lut_widths;
        // the angle for each theta discretization bin
        std::vector<float> angles;

        // compute useful constants and cache for later use
        for (int i = 0; i < theta_discretization; ++i) {
#if _USE_CACHED_CONSTANTS
            float angle = i * M_2PI_div_theta_discretization;
#else
            float angle = M_2PI * i / theta_discretization;
#endif
            angles.push_back(angle);

#if _USE_CACHED_TRIG == 1
            float cosfangle = cosf(angle);
            float sinfangle = sinf(angle);
            cos_values.push_back(cosfangle);
            sin_values.push_back(sinfangle);
#endif

// compute the height of the axis aligned bounding box, which will determine
// the necessary width of the lookup table for this angle
#if _USE_CACHED_TRIG == 1
            float rotated_height = std::abs(map.width() * sinfangle) + std::abs(map.height() * cosfangle);
#else
            float rotated_height = std::abs(map.width() * sinf(angle)) + std::abs(map.height() * cosf(angle));
#endif
            unsigned int lut_width = ceil(rotated_height - _EPSILON);
            lut_widths.push_back(lut_width);

/* the entire map will be rotated by the given angle. Every pixel in t hat map must be
   projected into the LUT, so we need to make sure that the index of every pixel will be
   positive when projected into LUT space. For example, here's the example with no rotation

    (0,height)  (width,height)      {
            *----------*    -----------> []
            |  a       |    -----------> [a]
            |      b   |    -----------> [b]
            |  c      d|    -----------> [c,d]
            *----------o    -----------> []
          (0,0)       (width,0)           }

     This is the case when theta = pi / 2 radians (not totally to scale)

   (-height,width) (0,width)      {
             *--------*    -----------> []
             |       d|    -----------> [d]
             |   b    |    -----------> [b]
             |        |    -----------> []
             | a   c  |    -----------> [a,c]
             *--------o    -----------> []
     (-height,0)  (0,0)                }

   Notably, the corner labeled 'o' lies on the origin no matter the rotation. Therefore,
   to ensure every LUT index is positive, we should translate the rotated map by:
             max(0, -1 * [minimum y coordinate for each corner])
*/

// this is the y-coordinate for each non-origin corner
#if _USE_CACHED_TRIG == 1
            float left_top_corner_y = map.height() * cosfangle;
            float right_top_corner_y = map.width() * sinfangle + map.height() * cosfangle;
            float right_bottom_corner_y = map.width() * sinfangle;
#else
            float left_top_corner_y = map.height() * cosf(angle);
            float right_top_corner_y = map.width() * sinf(angle) + map.height() * cosf(angle);
            float right_bottom_corner_y = map.width() * sinf(angle);
#endif

            // find the lowest corner, and determine the translation necessary to make them all positive
            float min_corner_y = std::min(left_top_corner_y, std::min(right_top_corner_y, right_bottom_corner_y));
            float lut_translation = std::max(0.0, -1.0 * min_corner_y - _EPSILON);

            lut_translations.push_back(lut_translation);
        }

        // build the empty LUT datastructure
        for (int a = 0; a < theta_discretization; ++a) {
            std::vector<std::vector<float>> projection_lut;
            for (int i = 0; i < lut_widths[a]; ++i)
            // for (int i = 0; i < 10; ++i)
            {
                std::vector<float> column;
                projection_lut.push_back(column);
            }
            compressed_lut.push_back(projection_lut);
        }

        // compute the edge map of the geometry - no ray can intersect with non-edge geometry,
        // so pruning it now will speed up LUT construction, especially for dense maps
        OMap edge_map = map.makeEdgeMap(true);
        // edge_map.save("./edge_map.png");

        // fill the LUT datastructure by projecting each occupied pixel into LUT space and storing
        // the x position in LUT space at the correct place as determined by y position and theta
        for (int x = 0; x < map.width(); ++x) {
            for (int y = 0; y < map.height(); ++y) {
                // if (map.isOccupied(x,y)) {
                if (edge_map.isOccupied(x, y)) {
                    // this (x,y) is occupied, so add it to the datastruture
                    std::pair<float, float> pixel_center = std::make_pair(x + 0.5, y + 0.5);
                    for (int a = 0; a < theta_discretization / 2.0; ++a) {
#if _USE_CACHED_TRIG == 1
                        float cosangle = cos_values[a];
                        float sinangle = sin_values[a];
#else
                        float angle = angles[a];
                        float cosangle = cosf(angle);
                        float sinangle = sinf(angle);
#endif

                        float half_lut_space_width = (std::abs(sinangle) + std::abs(cosangle)) / 2.0;

                        float lut_space_center_x = pixel_center.first * cosangle - pixel_center.second * sinangle;
                        float lut_space_center_y = (pixel_center.first * sinangle + pixel_center.second * cosangle) + lut_translations[a];

                        int upper_bin = lut_space_center_y + half_lut_space_width - _EPSILON;
                        int lower_bin = lut_space_center_y - half_lut_space_width + _EPSILON;

                        // the following is a quick hack to prevent problems in the cardinal directions
                        // where it has been known to see through walls
                        // if (std::fmod(angle, M_PI/2.0) < _EPSILON) {
                        // 	upper_bin++;
                        // 	lower_bin--;
                        // }

                        for (int i = lower_bin; i <= upper_bin; ++i)
                            compressed_lut[a][i].push_back(lut_space_center_x);

                        // std::cout << std::endl;
                        // std::cout << "angle: " << angle << std::endl;
                        // std::cout << "center: (" << pixel_center.first << ", " << pixel_center.second << ")" << std::endl;
                        // std::cout << "new center: (" << lut_space_center_x << ", " << lut_space_center_y << ")" << std::endl;
                        // std::cout << "bins:" << upper_bin << "    " << (int) lut_space_center_y << "   " << lower_bin << std::endl;
                        // std::cout << "width:" << half_lut_space_width << std::endl;
                        // std::cout << "trans" << lut_translations[a] << std::endl;
                        // std::cout << upper_bin << "   " << lower_bin << "   " << lut_translations[a] << std::endl;
                        // std::cout << lut_space_center_x << "  " << lut_space_center_y << std::endl;
                    }
                }
            }
        }

        // sort the vectors for faster lookup with binary search
        for (int a = 0; a < theta_discretization; ++a) {
            for (int i = 0; i < compressed_lut[a].size(); ++i) {
                // sort the vectors
                std::sort(compressed_lut[a][i].begin(), compressed_lut[a][i].end());

                // remove all duplicate entries, they will not change the answer
                compressed_lut[a][i].erase(unique(compressed_lut[a][i].begin(), compressed_lut[a][i].end()), compressed_lut[a][i].end());
            }
        }

#if _TRACK_LUT_SIZE == 1
        std::cout << "LUT SIZE (MB): " << lut_size() / 1000000.0 << std::endl;
#endif

#if _TRACK_COLLISION_INDEXES == 1
        for (int a = 0; a < theta_discretization; ++a) {
            std::vector<std::set<int>> projection_lut_tracker;
            for (int i = 0; i < lut_widths[a]; ++i) {
                std::set<int> collection;
                projection_lut_tracker.push_back(collection);
            }
            collision_table.push_back(projection_lut_tracker);
        }
#endif
    }

    int lut_size()
    {
        int lut_size = 0;
        // sort the vectors for faster lookup with binary search
        for (int a = 0; a < theta_discretization; ++a) {
            for (int i = 0; i < compressed_lut[a].size(); ++i) {
                lut_size += compressed_lut[a][i].size();
            }
        }
        return lut_size * sizeof(float);
    }

    int memory() { return lut_size() + map.memory() + lut_translations.size() * sizeof(float); }

    // mark all of the LUT entries that are potentially useful
    // remove all LUT entries that are not potentially useful
    void prune(float max_range)
    {
        std::vector<std::vector<std::set<int>>> local_collision_table;

        for (int a = 0; a < theta_discretization / 2.0; ++a) {
            std::vector<std::set<int>> projection_lut_tracker;
            for (int i = 0; i < compressed_lut[a].size(); ++i) {
                std::set<int> collection;
                projection_lut_tracker.push_back(collection);
            }
            local_collision_table.push_back(projection_lut_tracker);
        }

        for (int angle_index = 0; angle_index < theta_discretization / 2.0; ++angle_index) {
#if _USE_CACHED_CONSTANTS
            float angle = angle_index * M_2PI_div_theta_discretization;
#else
            float angle = M_2PI * angle_index / theta_discretization;
#endif

#if _USE_CACHED_TRIG == 1
            float cosangle = cos_values[angle_index];
            float sinangle = sin_values[angle_index];
#else
            float cosangle = cosf(angle);
            float sinangle = sinf(angle);
#endif

            // float cosangle = cos_values[angle_index];
            // float sinangle = sin_values[angle_index];
            float translation = lut_translations[angle_index];

            float lut_space_x;
            float lut_space_y;
            unsigned int lut_index;
            std::vector<float> *lut_bin;
            for (int x = 0; x < map.grid().size(); ++x) {
                float _x = 0.5 + x;
                for (int y = 0; y < map.grid()[0].size(); ++y) {
                    float _y = 0.5 + y;
                    lut_space_x = _x * cosangle - _y * sinangle;
                    lut_space_y = (_x * sinangle + _y * cosangle) + translation;
                    lut_index = (int)lut_space_y;

                    lut_bin = &compressed_lut[angle_index][lut_index];

                    // binary search for next greatest element
                    // int low = 0;
                    int high = lut_bin->size() - 1;

                    // there are no entries in this lut bin
                    if (high == -1)
                        continue;
                    if (map.isOccupied(x, y))
                        continue;

                    // the furthest entry is behind the query point
                    // if ((*lut_bin)[high] + max_range < lut_space_x) return std::make_pair(max_range, max_range);
                    if ((*lut_bin)[high] < lut_space_x && lut_space_x - (*lut_bin)[high] < max_range) {
                        local_collision_table[angle_index][lut_index].insert(high);
                        // accum += 1;
                        continue;
                    }

                    int index;
                    if (high > _BINARY_SEARCH_THRESHOLD) {
                        // once the binary search terminates, the next greatest element is indicated by 'val'
                        index = std::lower_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) - lut_bin->begin();
                    }
                    else {  // do linear search if array is very small
                        for (int i = 0; i < lut_bin->size(); ++i) {
                            if ((*lut_bin)[i] >= lut_space_x) {
                                index = i;
                                break;
                            }
                        }
                    }

                    int inverse_index = index - 1;
                    if (inverse_index == -1) {
                        local_collision_table[angle_index][lut_index].insert(index);
                        continue;
                    }
                    else {
                        local_collision_table[angle_index][lut_index].insert(index);
                        local_collision_table[angle_index][lut_index].insert(inverse_index);
                        continue;
                    }
                }
            }
        }

#if _TRACK_LUT_SIZE == 1
        std::cout << "OLD LUT SIZE (MB): " << lut_size() / 1000000.0 << std::endl;
#endif

        for (int a = 0; a < theta_discretization / 2.0; ++a) {
            for (int lut_index = 0; lut_index < compressed_lut[a].size(); ++lut_index) {
                std::vector<float> pruned_bin;

                for (int i = 0; i < compressed_lut[a][lut_index].size(); ++i) {
                    bool is_used = local_collision_table[a][lut_index].find(i) != local_collision_table[a][lut_index].end();
                    if (is_used)
                        pruned_bin.push_back(compressed_lut[a][lut_index][i]);
                }
                compressed_lut[a][lut_index] = pruned_bin;
            }
        }

#if _TRACK_LUT_SIZE == 1
        std::cout << "NEW LUT SIZE (MB): " << lut_size() / 1000000.0 << std::endl;
#endif
    }

    // takes a continuous theta space and returns the nearest theta in the discrete LUT space
    // as well as the bin index that the given theta falls into
    std::tuple<int, float, bool> discretize_theta(float theta)
    {
        theta = fmod(theta, M_2PI);
        // fmod does not wrap the angle into the positive range, so this will fix that if necessary
        if (theta < 0.0)
            theta += M_2PI;

        // exploit rotational symmetry by wrapping the theta range around to the range 0:pi
        bool is_flipped = false;
        if (theta >= M_PI) {
            is_flipped = true;
            theta -= M_PI;
        }

#if _USE_CACHED_CONSTANTS == 1
#if _USE_FAST_ROUND == 1
        int rounded = int(theta * theta_discretization_div_M_2PI + 0.5);
#else
        int rounded = (int)roundf(theta * theta_discretization_div_M_2PI);
#endif
        // this handles the special case where the theta rounds up and should wrap around
        if (rounded == theta_discretization >> 1) {
            rounded = 0;
            is_flipped = !is_flipped;
        }
        int binned = rounded % theta_discretization;
        float discrete_angle = binned * M_2PI_div_theta_discretization;
#else
#if _USE_FAST_ROUND == 1
        int rounded = int((theta * theta_discretization / M_2PI) + 0.5);
#else
        int rounded = (int)roundf(theta * theta_discretization / M_2PI);
#endif
        // this handles the special case where the theta rounds up and should wrap around
        if (rounded == theta_discretization >> 1) {
            rounded = 0;
            is_flipped = !is_flipped;
        }
        int binned = rounded % theta_discretization;
        float discrete_angle = (binned * M_2PI) / ((float)theta_discretization);
#endif
        return std::make_tuple(binned, discrete_angle, is_flipped);
    }

    float calc_range(float x, float y, float heading)
    {
#if _USE_LRU_CACHE
        int theta_key = (int)roundf(heading * theta_discretization_div_M_2PI);
        // int theta_key = angle_index;
        // if (is_flipped)
        // 	theta_key += theta_discretization/2;
        uint64_t key = key_maker.make_key(int(x), int(y), theta_key);
        float val = cache.get(key);
        if (val > 0) {
            hits += 1;
            return val;
        }
        else {
            misses += 1;
        }
// if (cache.exists(key)) { return cache.get(key); }
#endif

        int angle_index;
        float discrete_theta;
        bool is_flipped;
        std::tie(angle_index, discrete_theta, is_flipped) = discretize_theta(-1.0 * heading);

#if _USE_CACHED_TRIG == 1
        float cosangle = cos_values[angle_index];
        float sinangle = sin_values[angle_index];
#else
        float cosangle = cosf(discrete_theta);
        float sinangle = sinf(discrete_theta);
#endif

        float lut_space_x = x * cosangle - y * sinangle;
        float lut_space_y = (x * sinangle + y * cosangle) + lut_translations[angle_index];

        unsigned int lut_index = (int)lut_space_y;
        // this is to prevent segfaults
        if (lut_index < 0 || lut_index >= compressed_lut[angle_index].size())
            return max_range;
        std::vector<float> *lut_bin = &compressed_lut[angle_index][lut_index];

        // the angle is in range pi:2pi, so we must search in the opposite direction
        if (is_flipped) {
            // std::cout << "flipped" << std::endl;
            // binary search for next greatest element
            int low = 0;
            int high = lut_bin->size() - 1;

            // there are no entries in this lut bin
            if (high == -1) {
#if _USE_LRU_CACHE
                cache.put(key, max_range);
#endif
                return max_range;
            }
            // the furthest entry is behind the query point
            if ((*lut_bin)[low] > lut_space_x) {
#if _USE_LRU_CACHE
                cache.put(key, max_range);
#endif
                return max_range;
            }
            if ((*lut_bin)[high] < lut_space_x) {
                float val = lut_space_x - (*lut_bin)[high];
#if _USE_LRU_CACHE
                cache.put(key, val);
#endif
                return val;
            }

            // the query point is on top of a occupied pixel
            // this call is here rather than at the beginning, because it is apparently more efficient.
            // I presume that this has to do with the previous two return statements
            if (map.isOccupied(x, y)) {
                return 0.0;
            }

            if (high > _BINARY_SEARCH_THRESHOLD) {
                int index = std::upper_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) - lut_bin->begin();
                assert(index > 0);  // if index is 0, this will segfault. that should never happen, though.
                float val = lut_space_x - (*lut_bin)[index - 1];

#if _TRACK_COLLISION_INDEXES == 1
                collision_table[angle_index][lut_index].insert(index);
#endif

#if _USE_LRU_CACHE
                cache.put(key, val);
#endif
                return val;
            }
            else {  // do linear search if array is very small
                for (int i = high; i >= 0; --i) {
                    float obstacle_x = (*lut_bin)[i];
                    if (obstacle_x <= lut_space_x) {
#if _TRACK_COLLISION_INDEXES == 1
                        collision_table[angle_index][lut_index].insert(i);
#endif

                        float val = lut_space_x - obstacle_x;
#if _USE_LRU_CACHE
                        cache.put(key, val);
#endif
                        return val;
                    }
                }
            }
        }
        else {
            // std::cout << "not flipped" << std::endl;
            // binary search for next greatest element
            int low = 0;
            int high = lut_bin->size() - 1;

            // there are no entries in this lut bin
            if (high == -1) {
#if _USE_LRU_CACHE
                cache.put(key, max_range);
#endif
                return max_range;
            }
            // the furthest entry is behind the query point
            if ((*lut_bin)[high] < lut_space_x) {
#if _USE_LRU_CACHE
                cache.put(key, max_range);
#endif
                return max_range;
            }
            if ((*lut_bin)[low] > lut_space_x) {
                float val = (*lut_bin)[low] - lut_space_x;
#if _USE_LRU_CACHE
                cache.put(key, val);
#endif
                return val;
            }
            // the query point is on top of a occupied pixel
            // this call is here rather than at the beginning, because it is apparently more efficient.
            // I presume that this has to do with the previous two return statements
            if (map.isOccupied(x, y)) {
                return 0.0;
            }

            if (high > _BINARY_SEARCH_THRESHOLD) {
                // once the binary search terminates, the next greatest element is indicated by 'val'
                // float val = *std::lower_bound(lut_bin->begin(), lut_bin->end(), lut_space_x);
                int index = std::upper_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) - lut_bin->begin();
                float val = (*lut_bin)[index] - lut_space_x;

#if _TRACK_COLLISION_INDEXES == 1
                collision_table[angle_index][lut_index].insert(index);
#endif

#if _USE_LRU_CACHE
                cache.put(key, val);
#endif

                return val;
            }
            else {  // do linear search if array is very small
                // std::cout << "L" ;//<< std::endl;
                for (int i = 0; i < lut_bin->size(); ++i) {
                    float obstacle_x = (*lut_bin)[i];
                    if (obstacle_x >= lut_space_x) {
#if _TRACK_COLLISION_INDEXES == 1
                        collision_table[angle_index][lut_index].insert(i);
#endif

                        float val = obstacle_x - lut_space_x;

#if _USE_LRU_CACHE
                        cache.put(key, val);
#endif

                        return val;
                    }
                }
            }
        }
        // this should never occur, if it does, there's an error
        assert(0);
        return -1.0;
    }

    // returns both range for the given heading, and heading + pi/2
    // it is efficient to do both at the same time, rather than both
    // independently if they are both required
    std::pair<float, float> calc_range_pair(float x, float y, float heading)
    {
        int angle_index;
        float discrete_theta;
        bool is_flipped;
        std::tie(angle_index, discrete_theta, is_flipped) = discretize_theta(-1.0 * heading);

#if _USE_CACHED_TRIG == 1
        float cosangle = cos_values[angle_index];
        float sinangle = sin_values[angle_index];
#else
        float cosangle = cosf(discrete_theta);
        float sinangle = sinf(discrete_theta);
#endif

        float lut_space_x = x * cosangle - y * sinangle;
        float lut_space_y = (x * sinangle + y * cosangle) + lut_translations[angle_index];

        unsigned int lut_index = (int)lut_space_y;
        std::vector<float> *lut_bin = &compressed_lut[angle_index][lut_index];

        // the angle is in range pi:2pi, so we must search in the opposite direction
        if (is_flipped) {
            // std::cout << "is flipped" << std::endl;
            // binary search for next greatest element
            int low = 0;
            int high = lut_bin->size() - 1;

            // there are no entries in this lut bin
            if (high == -1)
                return std::make_pair(max_range, max_range);
            // the furthest entry is behind the query point and out of max range of the inverse query
            // if ((*lut_bin)[low] - max_range > lut_space_x) return std::make_pair(max_range, max_range);
            if ((*lut_bin)[low] > lut_space_x)
                return std::make_pair(max_range, std::min(max_range, (*lut_bin)[low] - lut_space_x));
            if ((*lut_bin)[high] < lut_space_x)
                return std::make_pair(lut_space_x - (*lut_bin)[high], max_range);
            // the query point is on top of a occupied pixel
            // this call is here rather than at the beginning, because it is apparently more efficient.
            // I presume that this has to do with the previous two return statements
            if (map.isOccupied(x, y)) {
                return std::make_pair(0.0, 0.0);
            }

            float val;
            int index;
            if (high > _BINARY_SEARCH_THRESHOLD) {
                // once the binary search terminates, the next least element is indicated by 'val'
                // float val = *std::lower_bound(lut_bin->begin(), lut_bin->end(), lut_space_x);
                index = std::upper_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) - lut_bin->begin() - 1;
                val = (*lut_bin)[index];
            }
            else {  // do linear search if array is very small
                for (int i = high; i >= 0; --i) {
                    float obstacle_x = (*lut_bin)[i];
                    if (obstacle_x <= lut_space_x) {
                        index = i;
                        val = obstacle_x;
                        break;
                    }
                }
            }

            int inverse_index = index + 1;
            if (inverse_index == lut_bin->size()) {
#if _TRACK_COLLISION_INDEXES == 1
                collision_table[angle_index][lut_index].insert(index);
#endif

                return std::make_pair(lut_space_x - val, max_range);
            }
            else {
#if _TRACK_COLLISION_INDEXES == 1
                collision_table[angle_index][lut_index].insert(index);
                collision_table[angle_index][lut_index].insert(inverse_index);
#endif

                return std::make_pair(lut_space_x - val, (*lut_bin)[inverse_index] - lut_space_x);
            }
        }
        else {
            // std::cout << "flipped" << std::endl;
            // binary search for next greatest element
            // int low = 0;
            int high = lut_bin->size() - 1;

            // there are no entries in this lut bin
            if (high == -1)
                return std::make_pair(max_range, max_range);
            // the furthest entry is behind the query point
            // if ((*lut_bin)[high] + max_range < lut_space_x) return std::make_pair(max_range, max_range);
            if ((*lut_bin)[high] < lut_space_x)
                return std::make_pair(max_range, std::min(max_range, lut_space_x - (*lut_bin)[high]));
            // TODO might need another early return case here
            // return std::make_pair(max_range, std::min(max_range, lut_space_x - (*lut_bin)[high]));
            // the query point is on top of a occupied pixel
            // this call is here rather than at the beginning, because it is apparently more efficient.
            // I presume that this has to do with the previous two return statements
            if (map.isOccupied(x, y)) {
                std::make_pair(0.0, 0.0);
            }

            float val;
            int index;
            if (high > _BINARY_SEARCH_THRESHOLD) {
                // once the binary search terminates, the next greatest element is indicated by 'val'
                // float val = *std::lower_bound(lut_bin->begin(), lut_bin->end(), lut_space_x);
                index = std::lower_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) - lut_bin->begin();
                val = (*lut_bin)[index];
            }
            else {  // do linear search if array is very small
                // std::cout << "L" ;//<< std::endl;
                for (int i = 0; i < lut_bin->size(); ++i) {
                    float obstacle_x = (*lut_bin)[i];
                    if (obstacle_x >= lut_space_x) {
                        val = obstacle_x;
                        index = i;
                        break;
                    }
                }
            }

            int inverse_index = index - 1;
            if (inverse_index == -1) {
#if _TRACK_COLLISION_INDEXES == 1
                collision_table[angle_index][lut_index].insert(index);
#endif

                return std::make_pair(val - lut_space_x, max_range);
            }
            else {
#if _TRACK_COLLISION_INDEXES == 1
                collision_table[angle_index][lut_index].insert(index);
                collision_table[angle_index][lut_index].insert(inverse_index);
#endif
                return std::make_pair(val - lut_space_x, lut_space_x - (*lut_bin)[inverse_index]);
            }
        }
    }

    // this works ok, but yaml deserialization is REALLY slow (at least in Python)
    void serializeYaml(std::stringstream *ss)
    {
        // (*ss) << std::fixed;
        (*ss) << std::setprecision(7);

        (*ss) << "cddt:" << std::endl;
        (*ss) << T1 << "theta_discretization: " << theta_discretization << std::endl;
        (*ss) << T1 << "lut_translations: ";
        utils::serialize(lut_translations, ss);
        (*ss) << std::endl;
        (*ss) << T1 << "max_range: " << max_range << std::endl;
        (*ss) << T1 << "map: " << std::endl;
        (*ss) << T2 << "# note: map data is width and then height (width is number of rows) transposed from expectation:" << std::endl;
        (*ss) << T2 << "path: " << map.filename() << std::endl;
        (*ss) << T2 << "width: " << map.width() << std::endl;
        (*ss) << T2 << "height: " << map.height() << std::endl;
        (*ss) << T2 << "data: " << std::endl;
        for (int i = 0; i < map.width(); ++i) {
            (*ss) << T3 << "- ";
            utils::serialize(map.grid()[i], ss);
            (*ss) << std::endl;
        }
        (*ss) << T1 << "compressed_lut: " << std::endl;
        for (int i = 0; i < compressed_lut.size(); ++i) {
#if _USE_CACHED_CONSTANTS
            float angle = i * M_2PI_div_theta_discretization;
#else
            float angle = M_2PI * i / theta_discretization;
#endif

            (*ss) << T2 << "- slice: " << std::endl;
            (*ss) << T3 << "theta: " << angle << std::endl;
            (*ss) << T3 << "zeros: " << std::endl;

            for (int j = 0; j < compressed_lut[i].size(); ++j) {
                (*ss) << T4 << "- ";
                utils::serialize(compressed_lut[i][j], ss);
                (*ss) << std::endl;
            }
        }
    }

    // directly generates JSON
    void serializeJson(std::stringstream *ss)
    {
        // (*ss) << std::fixed;
        (*ss) << std::setprecision(7);

        (*ss) << "{\"cddt\": {" << std::endl;
        (*ss) << J1 << "\"theta_discretization\": " << theta_discretization << "," << std::endl;
        (*ss) << J1 << "\"lut_translations\": ";
        utils::serialize(lut_translations, ss);
        (*ss) << "," << std::endl;
        (*ss) << J1 << "\"max_range\":" << max_range << "," << std::endl;
        (*ss) << J1 << "\"map\": {" << std::endl;
        // (*ss) << J2 << "# note: map data is width and then height (width is number of rows) transposed from expectation:"  << std::endl;
        (*ss) << J2 << "\"path\": \"" << map.filename() << "\"," << std::endl;
        (*ss) << J2 << "\"width\": " << map.width() << "," << std::endl;
        (*ss) << J2 << "\"height\": " << map.height() << "," << std::endl;

        (*ss) << J2 << "\"data\": [";  // utils::serialize(map.grid()[0], ss);
        for (int i = 0; i < map.width(); ++i) {
            if (i > 0)
                (*ss) << ",";
            utils::serialize(map.grid()[i], ss);
        }
        (*ss) << "]," << std::endl;
        (*ss) << J1 << "}," << std::endl;
        (*ss) << J1 << "\"compressed_lut\": [" << std::endl;
        for (int i = 0; i < compressed_lut.size(); ++i) {
#if _USE_CACHED_CONSTANTS
            float angle = i * M_2PI_div_theta_discretization;
#else
            float angle = M_2PI * i / theta_discretization;
#endif

            (*ss) << J2 << "{" << std::endl;
            (*ss) << J3 << "\"theta\": " << angle << "," << std::endl;
            (*ss) << J3 << "\"zeros\": [";

            for (int j = 0; j < compressed_lut[i].size(); ++j) {
                if (j > 0)
                    (*ss) << ",";
                utils::serialize(compressed_lut[i][j], ss);
            }
            (*ss) << "]" << std::endl;
            if (i == compressed_lut.size() - 1)
                (*ss) << J2 << "}" << std::endl;
            else
                (*ss) << J2 << "}," << std::endl;
        }
        (*ss) << J1 << "]" << std::endl;
        (*ss) << "}}" << std::endl;
    }

    void report()
    {
#if _USE_LRU_CACHE
        std::cout << "cache hits: " << hits << "  cache misses: " << misses << std::endl;
#endif
    }
    // protected:
    unsigned int theta_discretization;

    // compressed_lut[theta][offset] -> list of obstacle positions
    std::vector<std::vector<std::vector<float>>> compressed_lut;
    // cached list of y translations necessary to project points into lut space
    std::vector<float> lut_translations;

#if _USE_CACHED_TRIG == 1
    std::vector<float> cos_values;
    std::vector<float> sin_values;
#endif

#if _USE_CACHED_CONSTANTS == 1
    float theta_discretization_div_M_2PI;
    float M_2PI_div_theta_discretization;
#endif

#if _TRACK_COLLISION_INDEXES == 1
    std::vector<std::vector<std::set<int>>> collision_table;
#endif

#if _USE_LRU_CACHE
    cache::lru_cache<uint64_t, float> cache;
    utils::KeyMaker<uint64_t> key_maker;
    int hits = 0;
    int misses = 0;
#endif
};

class GiantLUTCast : public RangeMethod {
   public:
#if _GIANT_LUT_SHORT_DATATYPE
    typedef uint16_t lut_t;
#else
    typedef float lut_t;
#endif

    /// @brief Returns a LUT index based on the discretized inputs
    inline size_t getLutIdx(size_t x, size_t y, size_t theta)
    {
        return theta +
               y * _thetaDiscretization +
               x * _thetaDiscretization * _height;
    }

    /// @brief Constructor.
    /// @param[in] m Input Occupancy Grid Map
    /// @param[in] mr Max range
    /// @param[in] td theta discretization
    GiantLUTCast(OMap m, float mr, int td) : _height{m.height()},
                                             _width{m.width()},
                                             _thetaDiscretization{td},
                                             RangeMethod(m, mr),
                                             _thetaDiscretization_div_M_2PI{(float)_thetaDiscretization / M_2PI},
                                             _M_2PI_div_thetaDiscretization{M_2PI / (float)_thetaDiscretization},
                                             _maxDivLimits{max_range / std::numeric_limits<uint16_t>::max()},
                                             _limitsDivMax{std::numeric_limits<uint16_t>::max() / max_range}
    {
        // To initalize the LUT, use Ray Marching
        RayMarching seed_cast = RayMarching(m, mr);

        _GiantLUT.resize(_height * _width * _thetaDiscretization);
        for (int x = 0; x < _width; ++x) {
            for (int y = 0; y < _height; ++y) {
                for (int i = 0; i < _thetaDiscretization; ++i) {
                    float angle = i * _M_2PI_div_thetaDiscretization;
                    float r = seed_cast.calc_range(x, y, angle);

#if _GIANT_LUT_SHORT_DATATYPE
                    r = std::min(max_range, r);
                    uint16_t val = r * _limitsDivMax;
                    _GiantLUT[getLutIdx(x, y, i)] = val;
#else
                    _GiantLUT[getLutIdx(x, y, i)] = r;
#endif
                }
            }
        }
    }

    int lut_size()
    {
        return map.width() * map.height() * _thetaDiscretization * sizeof(lut_t);
    }

    int memory() { return lut_size(); }

    /// @brief takes a continuous theta space and returns the nearest theta in the discrete LUT space
    /// as well as the bin index that the given theta falls into
    /// @param[in] theta the angle
    /// @return the discretized value
    int discretize_theta(float theta)
    {
        theta = fmod(theta, M_2PI);
        // fmod does not wrap the angle into the positive range, so this will fix that if necessary
        if (theta < 0.0)
            theta += M_2PI;

#if _USE_FAST_ROUND == 1
        int rounded = int(theta * _thetaDiscretization_div_M_2PI + 0.5);
#else
        int rounded = (int)roundf(theta * _thetaDiscretization_div_M_2PI);
#endif
        int binned = rounded % _thetaDiscretization;

        return binned;
    }

    float calc_range(float x, float y, float heading)
    {
        if (x < 0 || x >= map.width() || y < 0 || y >= map.height()) {
            return max_range;
        }

        size_t idx = getLutIdx((size_t)x, (size_t)y, discretize_theta(heading));
#if _GIANT_LUT_SHORT_DATATYPE
        return _GiantLUT[idx] * _maxDivLimits;
#else
        return _GiantLUT[idx];
#endif
    }

    /// @brief Returns a slice of the LUT based on a theta value
    /// @param[in] theta the slice value
    /// @return a DistanceTransform unique_ptr
    std::unique_ptr<DistanceTransform> get_slice(float theta)
    {
        std::unique_ptr<DistanceTransform> slice = std::make_unique<DistanceTransform>(_width, _height);
        const int dtheta = discretize_theta(theta);

        size_t idx = dtheta;  // iterator index

        for (int x = 0; x < _width; ++x) {
            for (int y = 0; y < _height; ++y) {
                slice->grid()[x][y] = _GiantLUT[idx];
                idx += _thetaDiscretization;
            }
            idx += _thetaDiscretization * _height;
        }

        return slice;
    }

   protected:
    const int _height, _width;                   ///< dimensions of omap
    const int _thetaDiscretization;              ///< how many values of theta are inside
    const float _thetaDiscretization_div_M_2PI;  ///< theta_discretization / 2pi
    const float _M_2PI_div_thetaDiscretization;  ///< 2pi / theta_discretization
    const float _maxDivLimits;                   ///< Max range / uint16_t max
    const float _limitsDivMax;                   ///< uint16_t max / Max range

    std::vector<lut_t> _GiantLUT;  ///< actual LUT
};

}  // namespace ranges

namespace benchmark {
template <class range_T>
class Benchmark {
   public:
    Benchmark(range_T rm) : range(rm)
    {
        map = range.getMap();
    };
    ~Benchmark(){};

    void set_log(std::stringstream *ss) { log = ss; }

    int memory() { return range.memory(); }

    void grid_sample(int step_size, int num_rays, int samples)
    {
        float coeff = (2.0 * M_PI) / num_rays;
        double t_accum = 0;
        float num_cast = 0;

        volatile clock_t t;
        t = clock();

        if (log)
            (*log) << "x,y,theta,time" << std::endl;
        if (log)
            (*log) << std::fixed;
        if (log)
            (*log) << std::setprecision(9);

        for (int i = 0; i < num_rays; ++i) {
            float angle = i * coeff;
            for (int x = 0; x < map->width(); x += step_size) {
                for (int y = 0; y < map->height(); y += step_size) {
                    auto start_time = std::chrono::high_resolution_clock::now();
                    for (int j = 0; j < samples; ++j) {
                        volatile float r = range.calc_range(x, y, angle);
                    }

                    auto end_time = std::chrono::high_resolution_clock::now();

                    num_cast += samples;
                    // std::cout << (end_time - start_time).count() << std::endl;

                    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

                    t_accum += time_span.count();

                    if (log)
                        (*log) << x << "," << y << "," << angle << "," << time_span.count() << std::endl;
                }
            }
        }

        std::cout << "finished grid sample after: " << (((float)(clock() - t)) / CLOCKS_PER_SEC) << " sec" << std::endl;
        std::cout << " -avg time per ray: " << (t_accum / (float)num_cast) << " sec" << std::endl;
        std::cout << " -rays cast: " << num_cast << std::endl;
        std::cout << " -total time: " << t_accum << " sec" << std::endl;
    }

    void grid_sample2(int step_size, int num_rays, int samples)
    {
        float coeff = (2.0 * M_PI) / num_rays;
        double t_accum = 0;

        if (log)
            (*log) << "x,y,theta,time" << std::endl;
        if (log)
            (*log) << std::fixed;
        if (log)
            (*log) << std::setprecision(9);

        int num_samples = num_grid_samples(step_size, num_rays, samples, map->width(), map->height());
        float *samps = new float[num_samples * 3];
        float *outs = new float[num_samples];
        get_grid_samples(samps, step_size, num_rays, samples, map->width(), map->height());

        auto start_time = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_samples; ++i)
            outs[i] = range.calc_range(samps[i * 3], samps[i * 3 + 1], samps[i * 3 + 2]);
        auto end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);
        t_accum += time_span.count();

        // print first few outputs for sanity checking
        for (int i = 0; i < 10; ++i)
            std::cout << outs[i] << std::endl;

        std::cout << "finished grid sample after: " << (float)t_accum << " sec" << std::endl;
        std::cout << " -avg time per ray: " << (t_accum / (float)num_samples) << " sec" << std::endl;
        std::cout << " -rays cast: " << num_samples << std::endl;
        std::cout << " -total time: " << t_accum << " sec" << std::endl;
    }

    static int num_grid_samples(int step_size, int num_rays, int samples, int map_width, int map_height)
    {
        int num_samples = 0;
        for (int i = 0; i < num_rays; ++i)
            for (int x = 0; x < map_width; x += step_size)
                for (int y = 0; y < map_height; y += step_size)
                    for (int j = 0; j < samples; ++j)
                        num_samples++;
        return num_samples;
    }

    static void get_grid_samples(float *queries, int step_size, int num_rays, int samples, int map_width, int map_height)
    {
        float coeff = (2.0 * M_PI) / num_rays;
        double t_accum = 0;
        int ind = 0;
        for (int i = 0; i < num_rays; ++i) {
            float angle = i * coeff;
            for (int x = 0; x < map_width; x += step_size) {
                for (int y = 0; y < map_height; y += step_size) {
                    for (int j = 0; j < samples; ++j) {
                        queries[ind * 3] = (float)x;
                        queries[ind * 3 + 1] = (float)y;
                        queries[ind * 3 + 2] = angle;
                        ind++;
                    }
                }
            }
        }
    }

    void random_sample(int num_samples)
    {
        std::default_random_engine generator;
        generator.seed(clock());
        std::uniform_real_distribution<float> randx = std::uniform_real_distribution<float>(1.0, map->width() - 1.0);
        std::uniform_real_distribution<float> randy = std::uniform_real_distribution<float>(1.0, map->height() - 1.0);
        std::uniform_real_distribution<float> randt = std::uniform_real_distribution<float>(0.0, M_2PI);

        double t_accum = 0;
        for (int i = 0; i < num_samples; ++i) {
            float x = randx(generator);
            float y = randy(generator);
            float angle = randt(generator);

            auto start_time = std::chrono::high_resolution_clock::now();
            volatile float r = range.calc_range(x, y, angle);
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

            t_accum += time_span.count();
            if (log)
                (*log) << x << "," << y << "," << angle << "," << time_span.count() << std::endl;
        }

        std::cout << "finished random sample after: " << t_accum << " sec" << std::endl;
        std::cout << " -avg time per ray: " << (t_accum / (float)num_samples) << " sec" << std::endl;
        std::cout << " -rays cast: " << num_samples << std::endl;
    }

    static void get_random_samples(float *queries, int num_samples, int map_width, int map_height)
    {
        std::default_random_engine generator;
        generator.seed(std::chrono::duration_cast<std::chrono::nanoseconds>(
                           std::chrono::system_clock::now().time_since_epoch())
                           .count());
        std::uniform_real_distribution<float> randx = std::uniform_real_distribution<float>(1.0, map_width - 1.0);
        std::uniform_real_distribution<float> randy = std::uniform_real_distribution<float>(1.0, map_height - 1.0);
        std::uniform_real_distribution<float> randt = std::uniform_real_distribution<float>(0.0, M_2PI);

        for (int i = 0; i < num_samples; ++i) {
            queries[3 * i] = randx(generator);
            queries[3 * i + 1] = randy(generator);
            queries[3 * i + 2] = randt(generator);
        }
    }

    ranges::OMap *getMap() { return range.getMap(); }

   protected:
    range_T range;
    ranges::OMap *map;
    std::stringstream *log = NULL;
};
}  // namespace benchmark

#endif
