#ifndef RANGELIB_RANGE_METHOD_HPP_
#define RANGELIB_RANGE_METHOD_HPP_

#include "rangelib/omap.hpp"

#define _EPSILON 0.00001
#define M_2PI 2 * M_PI

// TODO fix with updated API

namespace ranges {
/// @brief Parent class of all range methods
class RangeMethod {
   public:
    RangeMethod(OMap m, float mr) : map(m), max_range(mr){};

    virtual ~RangeMethod(){};

    virtual float calc_range(float x, float y, float heading) = 0;

    virtual std::pair<float, float> calc_range_pair(float x, float y, float heading)
    {
        return std::make_pair(-1, -1);
    }

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

    void numpy_calc_range_angles(float *ins, float *angles, float *outs, int num_particles,
                                 int num_angles)
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
            for (int j = 0; j < table_width; ++j) table_row.push_back(table[table_width * i + j]);
            sensor_model.push_back(table_row);
        }
    }

    void eval_sensor_model(float *obs, float *ranges, double *outs, int rays_per_particle,
                           int particles)
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
    void calc_range_repeat_angles_eval_sensor_model(float *ins, float *angles, float *obs,
                                                    double *weights, int num_particles,
                                                    int num_angles)
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
    void calc_range_many_radial_optimized(float *ins, float *outs, int num_particles, int num_rays,
                                          float min_angle, float max_angle)
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
    OMap map;  // TODO perhaps dont store an OMap and a DistanceTransform simultaneously. Or just
               // store a DistanceTransform?
    float max_range;

    std::vector<std::vector<double>> sensor_model;
};

}  // namespace ranges

#endif