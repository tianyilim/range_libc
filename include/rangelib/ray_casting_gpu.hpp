#ifndef RANGELIB_RAY_MARCHING_GPU_HPP_
#define RANGELIB_RAY_MARCHING_GPU_HPP_

#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"

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

    int memory() { return distImage->memory(); }

   protected:
    DistanceTransform *distImage = 0;
#if USE_CUDA == 1
    RayMarchingCUDA *rmc = 0;
#endif
    bool already_warned = false;
};

}  // namespace ranges

#endif