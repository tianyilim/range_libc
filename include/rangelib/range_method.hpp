#ifndef RANGELIB_RANGE_METHOD_HPP_
#define RANGELIB_RANGE_METHOD_HPP_

#include "rangelib/omap.hpp"

#define _EPSILON 0.00001
#define M_2PI 2 * M_PI

namespace ranges {

using Pose2Df_t = std::tuple<float, float, float>;

/// @brief Parent class of all range methods
class RangeMethod {
   public:
    /// @brief Constructor.
    /// @param[in] m Input map
    /// @param[in] mr Max range of sensor
    RangeMethod(DistanceTransform map, float maxRange)
        : _distTransform(map),
          _maxRange(maxRange),
          _worldScale{_distTransform.worldValues().worldScale},
          _invWorldScale{1.0f / _worldScale},
          _worldAngle{_distTransform.worldValues().worldAngle},
          _worldOriginX{_distTransform.worldValues().worldOriginX},
          _worldOriginY{_distTransform.worldValues().worldOriginY},
          _worldSinAngle{_distTransform.worldValues().worldSinAngle},
          _worldCosAngle{_distTransform.worldValues().worldCosAngle},
          _mapOriginX{_worldOriginX * _worldCosAngle - _worldOriginY * _worldSinAngle},
          _mapOriginY{(_worldOriginX * _worldSinAngle + _worldOriginY * _worldCosAngle) +
                      _distTransform.height() * _worldScale}
    {
    }

    virtual ~RangeMethod(){};

    /// @brief Copy Constructor
    /// @param r Other
    RangeMethod(const RangeMethod &r)
        : _distTransform(r._distTransform),
          _maxRange(r._maxRange),
          _worldScale{r._worldScale},
          _invWorldScale{r._invWorldScale},
          _worldAngle{r._worldAngle},
          _worldOriginX{r._worldOriginX},
          _worldOriginY{r._worldOriginY},
          _worldSinAngle{r._worldSinAngle},
          _worldCosAngle{r._worldCosAngle},
          _mapOriginX{r._mapOriginX},
          _mapOriginY{r._mapOriginY},
          _sensorModel{r._sensorModel}
    {
    }

    /// @brief Move constructor
    RangeMethod(RangeMethod &&r)
        : _maxRange(r._maxRange),
          _worldScale{r._worldScale},
          _invWorldScale{r._invWorldScale},
          _worldAngle{r._worldAngle},
          _worldOriginX{r._worldOriginX},
          _worldOriginY{r._worldOriginY},
          _worldSinAngle{r._worldSinAngle},
          _worldCosAngle{r._worldCosAngle},
          _mapOriginX{r._mapOriginX},
          _mapOriginY{r._mapOriginY}
    {
        _distTransform = std::move(r._distTransform);
        _sensorModel = std::move(r._sensorModel);
    }

    /// @brief Base function to give the range measurement at this pose.
    virtual float calc_range(float x, float y, float heading) const = 0;

    /// @brief Overridden in CDDT.
    virtual std::pair<float, float> calc_range_pair([[maybe_unused]] float x,
                                                    [[maybe_unused]] float y,
                                                    [[maybe_unused]] float heading) const
    {
        return std::make_pair(-1, -1);
    }

    /// @brief OMap Accessor
    const OMap &getMap() const { return _distTransform; }

    virtual void report(){};

    float maxRange() const { return _maxRange; }

    virtual int memory() const = 0;

    /// @brief Wrapper function to call calc_range repeatedly. Indexing assumes a 3xn numpy array
    /// for the inputs and a 1xn numpy array of the outputs.
    /// @param[in] ins (3 x num_casts) array of poses.
    /// @param[out] outs (1 x num_casts) array of output ranges
    /// @param[in] num_casts
    void numpy_calc_range(float *ins, float *outs, int num_casts)
    {
        // avoid allocation on every loop iteration
        float x, y, theta;

        for (int i = 0; i < num_casts; ++i) {
            std::tie(x, y, theta) = mapCoordinates(ins[i * 3], ins[i * 3 + 1], ins[i * 3 + 2]);
            outs[i] = calc_range(x, y, theta) * _worldScale;
        }
    }

    /// @brief Wrapper function to call calc_range repeatedly.
    /// @param[in] ins (3 x num_particles) array of poses (of sensor)
    /// @param[in] angles (1 x num_angles) array of query angles
    /// @param[out] outs (1 x num_angles*num_particles) array of ranges. Organized with
    /// num_particles as major index.
    /// @param[in] num_particles
    /// @param[in] num_angles
    void numpy_calc_range_angles(float *ins, float *angles, float *outs, int num_particles,
                                 int num_angles)
    {
        // avoid allocation on every loop iteration
        float x, y, theta;

        for (int i = 0; i < num_particles; ++i) {
            std::tie(x, y, theta) = mapCoordinates(ins[i * 3], ins[i * 3 + 1], ins[i * 3 + 2]);

            for (int a = 0; a < num_angles; ++a) {
                outs[i * num_angles + a] = calc_range(x, y, theta - angles[a]) * _worldScale;
            }
        }
    }

    /// @brief Sets a pre-calculated sensor model from numpy. The sensor model is the likelihood
    /// score of observing a particular measurement.
    /// @param[in] table
    /// @param[in] table_width
    void set_sensor_model(double *table, int table_width)
    {
        // TODO convert into 1d vector
        // convert the sensor model from a numpy array to a vector array
        for (int i = 0; i < table_width; ++i) {
            std::vector<double> table_row;
            for (int j = 0; j < table_width; ++j) {
                table_row.push_back(table[table_width * i + j]);
            }
            _sensorModel.push_back(table_row);
        }
    }

    /// @brief Evaluating the (discretized) sensor model is equivalent to a 2D array lookup.
    /// @param[in] obs
    /// @param[in] ranges
    /// @param[in] outs
    /// @param[in] rays_per_particle
    /// @param[in] particles
    void eval_sensor_model(float *obs, float *ranges, double *outs, int rays_per_particle,
                           int num_particles) const
    {
        double weight;
        float r;  // observed ranges
        float d;  // expected values

        for (int i = 0; i < num_particles; ++i) {
            weight = 1.0;
            for (int j = 0; j < rays_per_particle; ++j) {
                // Convert into map scale
                r = obs[j] * _invWorldScale;
                d = ranges[i * rays_per_particle + j] * _invWorldScale;

                // Clamp
                r = std::min<float>(std::max<float>(r, 0.0), (float)_sensorModel.size() - 1.0);
                d = std::min<float>(std::max<float>(d, 0.0), (float)_sensorModel.size() - 1.0);

                // Discretize and evaluate
                weight *= _sensorModel[(int)r][(int)d];
            }
            outs[i] = weight;  // weight each particle
        }
    }

    /// @brief Function that calculates expected ranges for each particle and directly outputs
    /// particle weights
    /// @param ins
    /// @param angles
    /// @param obs
    /// @param weights
    /// @param num_particles
    /// @param num_angles
    void calc_range_repeat_angles_eval_sensor_model(float *ins, float *angles, float *obs,
                                                    double *weights, int num_particles,
                                                    int num_angles)
    {
        float x, y, theta;
        double weight;
        float r;  // observed ranges
        float d;  // expected values

        for (int particleIdx = 0; particleIdx < num_particles; ++particleIdx) {
            std::tie(x, y, theta) = mapCoordinates(ins[particleIdx * 3], ins[particleIdx * 3 + 1],
                                                   ins[particleIdx * 3 + 2]);

            weight = 1.0;
            for (int angleIdx = 0; angleIdx < num_angles; ++angleIdx) {
                // Convert into map scale
                r = obs[angleIdx] * _invWorldScale;
                d = calc_range(x, y, theta - angles[angleIdx]);

                // Clamp
                r = std::min<float>(std::max<float>(r, 0.0), (float)_sensorModel.size() - 1.0);
                d = std::min<float>(std::max<float>(d, 0.0), (float)_sensorModel.size() - 1.0);

                // Discretize and evaluate
                weight *= _sensorModel[(int)r][(int)d];
            }
            weights[particleIdx] = weight;
        }
    }

    // this is to compute a lidar sensor model using radial (calc_range_pair) optimizations
    // this is only exact for a certain set of downsample amounts
    [[deprecated("Currently unused, waiting for CDDT refactoring")]] void
    calc_range_many_radial_optimized(float *ins, float *outs, int num_particles, int num_rays,
                                     float min_angle, float max_angle)
    {
        float x, y, theta = 0.0;

        float step = (max_angle - min_angle) / (num_rays - 1);
        float angle = min_angle;

        int max_pairwise_index = (float)num_rays / 3.0;
        float index_offset_float = (num_rays - 1.0) * M_PI / (max_angle - min_angle);

        // TODO: check if this index_offset_float is not very close to an integer
        // in which case throw a warning that this downsample factor is poorly compatible
        // with this radial optimization
        int index_offset = roundf(index_offset_float);
        float r, r_inv = 0.0;

        for (int i = 0; i < num_particles; ++i) {
            std::tie(x, y, theta) = mapCoordinates(ins[i * 3], ins[i * 3 + 1], ins[i * 3 + 2]);

            angle = min_angle;
            for (int a = 0; a <= max_pairwise_index; ++a) {
                std::tie(r, r_inv) = calc_range_pair(y, x, theta - angle);
                outs[i * num_rays + a] = r * _worldScale;
                outs[i * num_rays + a + index_offset] = r_inv * _worldScale;
                angle += step;
            }

            for (int a = max_pairwise_index + 1; a < index_offset; ++a) {
                outs[i * num_rays + a] = calc_range(x, y, theta - angle) * _worldScale;
                angle += step;
            }
        }
        return;
    }

    /// @brief Convert world coordinates into map-discretized values
    inline Pose2Df_t mapCoordinates(float x_world, float y_world, float theta_world)
    {
        float thetaMap = -(theta_world - _worldAngle);
        float xMap =
            ((x_world * _worldCosAngle - y_world * _worldSinAngle) - _mapOriginX) * _invWorldScale;
        float yMap =
            (_mapOriginY - (x_world * _worldSinAngle + y_world * _worldCosAngle)) * _invWorldScale;

        return {xMap, yMap, thetaMap};
    }

   protected:
    DistanceTransform _distTransform;
    const float _maxRange;

    /// @brief map-specific parameters
    const float _worldScale;     ///< Resolution of map (m/px)
    const float _invWorldScale;  ///< Resolution (px/m)
    const float _worldAngle;     ///< Yaw of map origin relative to image frame
    const float _worldOriginX;   ///< x coords of bottom left px of map in world frame
    const float _worldOriginY;   ///< y coords of bottom left px of map in world frame
    const float _worldSinAngle;  ///< sin(_worldAngle)
    const float _worldCosAngle;  ///< cos(_worldAngle)
    const float _mapOriginX;  ///< x coords of top left px of map (rangelibc origin) in world frame
    const float _mapOriginY;  ///< y coords of top left px of map (rangelibc origin) in world frame

    std::vector<std::vector<double>> _sensorModel;
};

}  // namespace ranges

#endif