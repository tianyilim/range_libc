/*
* Copyright 2017 Corey H. Walsh (corey.walsh11@gmail.com)
  Modified by Tianyi Lim (0.tianyi.lim@gmail.com)

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

#include "rangelib/bresenhams.hpp"
#include "rangelib/cddt.hpp"
#include "rangelib/lookup_table.hpp"
#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"
#include "rangelib/ray_casting.hpp"
#include "rangelib/ray_casting_gpu.hpp"

namespace benchmark {
template <class range_T>
class Benchmark {
   public:
    /// @brief benchmark constructor.
    /// @param rm: m_Range method
    Benchmark(range_T rm) : m_Range(rm)
    {
        m_OMap = m_Range.getMap();
    };

    ~Benchmark(){};

    /// @brief Set logfile to be a stringstream pointer.
    /// @param ss Stringstream pointer
    void set_log(std::stringstream *ss) { m_Log = ss; }

    /// @brief Query the memory requirement of given range method.
    int memory() { return m_Range.memory(); }

    void grid_sample(int step_size, int num_rays, int samples)
    {
        float coeff = (2.0 * M_PI) / num_rays;
        double t_accum = 0;
        float num_cast = 0;

        volatile clock_t t;
        t = clock();

        if (m_Log)
            (*m_Log) << "x,y,theta,time" << std::endl;
        if (m_Log)
            (*m_Log) << std::fixed;
        if (m_Log)
            (*m_Log) << std::setprecision(9);

        for (int i = 0; i < num_rays; ++i) {
            float angle = i * coeff;
            for (int x = 0; x < m_OMap->width(); x += step_size) {
                for (int y = 0; y < m_OMap->height(); y += step_size) {
                    auto start_time = std::chrono::high_resolution_clock::now();
                    for (int j = 0; j < samples; ++j) {
                        volatile float r = m_Range.calc_range(x, y, angle);
                    }

                    auto end_time = std::chrono::high_resolution_clock::now();

                    num_cast += samples;
                    // std::cout << (end_time - start_time).count() << std::endl;

                    std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

                    t_accum += time_span.count();

                    if (m_Log)
                        (*m_Log) << x << "," << y << "," << angle << "," << time_span.count() << std::endl;
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

        if (m_Log)
            (*m_Log) << "x,y,theta,time" << std::endl;
        if (m_Log)
            (*m_Log) << std::fixed;
        if (m_Log)
            (*m_Log) << std::setprecision(9);

        int num_samples = num_grid_samples(step_size, num_rays, samples, m_OMap->width(), m_OMap->height());
        float *samps = new float[num_samples * 3];
        float *outs = new float[num_samples];
        get_grid_samples(samps, step_size, num_rays, samples, m_OMap->width(), m_OMap->height());

        auto start_time = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < num_samples; ++i)
            outs[i] = m_Range.calc_range(samps[i * 3], samps[i * 3 + 1], samps[i * 3 + 2]);
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
        std::uniform_real_distribution<float> randx = std::uniform_real_distribution<float>(1.0, m_OMap->width() - 1.0);
        std::uniform_real_distribution<float> randy = std::uniform_real_distribution<float>(1.0, m_OMap->height() - 1.0);
        std::uniform_real_distribution<float> randt = std::uniform_real_distribution<float>(0.0, M_2PI);

        double t_accum = 0;
        for (int i = 0; i < num_samples; ++i) {
            float x = randx(generator);
            float y = randy(generator);
            float angle = randt(generator);

            auto start_time = std::chrono::high_resolution_clock::now();
            volatile float r = m_Range.calc_range(x, y, angle);
            auto end_time = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(end_time - start_time);

            t_accum += time_span.count();
            if (m_Log)
                (*m_Log) << x << "," << y << "," << angle << "," << time_span.count() << std::endl;
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

    ranges::OMap *getMap() const { return m_Range.getMap(); }

   protected:
    range_T m_Range;                  ///< Range method
    ranges::OMap *m_OMap;             ///< Map
    std::stringstream *m_Log = NULL;  ///< Logging output
};
}  // namespace benchmark

#endif
