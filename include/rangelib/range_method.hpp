#ifndef RANGELIB_RANGE_METHOD_HPP_
#define RANGELIB_RANGE_METHOD_HPP_

#include <cstddef>
#include <stdexcept>

#include "eigen3/Eigen/Dense"
#include "rangelib/omap.hpp"

#define _EPSILON 0.00001
#define M_2PI 2 * M_PI

namespace ranges {

using Pose2Df_t = std::tuple<float, float, float>;
using VectorFloat = Eigen::VectorXf;
using VectorDouble = Eigen::VectorXd;
using PosesFloat = Eigen::Matrix<float, 3, Eigen::Dynamic, Eigen::RowMajor>;
using PosesDouble = Eigen::Matrix<double, 3, Eigen::Dynamic, Eigen::RowMajor>;

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

    /// @brief Wrapper function to call calc_range repeatedly.

    /// @param[in] ins (3 x N) array of poses.
    /// @return (1xN) array of ranges.
    VectorFloat batchCalcRange(const Eigen::Ref<const PosesFloat> poses) const
    {
        if (poses.rows() != 3) {
            throw std::invalid_argument("Poses input must be 3xN");
        }

        float x, y, theta;
        VectorFloat outs(poses.cols());

        for (int i = 0; i < poses.cols(); ++i) {
            const auto &currPose = poses.col(i);
            std::tie(x, y, theta) = mapCoordinates(currPose(0), currPose(1), currPose(2));
            outs(i) = calc_range(x, y, theta) * _worldScale;
        }

        return outs;
    }

    ///@brief Wrapper function to call calc_range repeatedly on a range of poses and angles.
    /// Returns in pose-major order, ie. for input `[P1 P2], [A1 A2 A3]`, the output will be ordered
    /// like `[P1A1 P1A2 P1A3 P2A1 P2A2 P2A3]`.
    ///
    ///@param poses (3 x N) array of poses.
    ///@param angles (1 x M) array of query angles.
    ///@return VectorFloat (1xNM) array of ranges in "pose-major" order.
    VectorFloat batchCalcRangeAngles(const Eigen::Ref<const PosesFloat> poses,
                                     const Eigen::Ref<const VectorFloat> angles) const
    {
        if (poses.rows() != 3) {
            throw std::invalid_argument("Poses input must be 3xN");
        }

        float x, y, theta;
        const auto numPoses = poses.cols();
        const auto numAngles = angles.size();
        VectorFloat outs(numPoses * numAngles);

        for (int i = 0; i < numPoses; ++i) {
            const auto &currPose = poses.col(i);
            std::tie(x, y, theta) = mapCoordinates(currPose(0), currPose(1), currPose(2));

            for (int j = 0; j < numAngles; ++j) {
                const auto currAngleOffset = angles(j);
                outs(i * numAngles + j) = calc_range(x, y, theta - currAngleOffset) * _worldScale;
            }
        }

        return outs;
    }

    ///@brief Sets a pre-calculated sensor model from numpy. The sensor model is the likelihood
    /// score of observing a particular measurement.
    ///
    ///@throws std::invalid argument if input is not square
    ///@param[in] sensorModel a NxN matrix
    void setSensorModel(
        const Eigen::Ref<const Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic>> sensorModel)
    {
        if (sensorModel.rows() != sensorModel.cols()) {
            throw std::invalid_argument(
                "Sensor Model must have the same number of rows and columns!");
        }

        _sensorModel = sensorModel;
    }

    ///@brief Evaluates how likely an expected range is, based on the observed range and the sensor
    /// model. Equivalent to a 2d array lookup.
    ///
    ///@param observedRanges (1 x N) array of ranges observed by a sensor.
    ///@param expectedRanges (1 x N) array of ranges computed by rangelib.
    ///@throws std::invalid argument if input arrays do not have same size
    ///@throws std::runtime_error if sensor model has not been initialized
    ///@return a likelihood as a float.
    float evalSensorModel(const Eigen::Ref<const VectorFloat> observedRanges,
                          const Eigen::Ref<const VectorFloat> expectedRanges) const
    {
        if (observedRanges.size() != expectedRanges.size()) {
            throw std::invalid_argument("Observations and Ranges must have same dimensions");
        }

        if (_sensorModel.rows() == 0) {
            throw std::runtime_error("The sensor model has size 0. Has it been initialized?");
        }

        const auto numRanges = observedRanges.size();

        float weight = 1.0;
        float observedRange, expectedRange;
        for (int i = 0; i < numRanges; ++i) {
            // Convert to map scale
            observedRange = observedRanges(i) * _invWorldScale;
            expectedRange = expectedRanges(i) * _invWorldScale;

            // Clamp to size of sensor model array
            observedRange = std::min(std::max(observedRange, 0.0f), float(_sensorModel.rows() - 1));
            expectedRange = std::min(std::max(expectedRange, 0.0f), float(_sensorModel.cols() - 1));

            weight *= _sensorModel(int(observedRange), int(expectedRange));
        }

        return weight;
    }

    ///@brief Evaluates how likely a set of expected ranges corresponding multiple particles are,
    /// based on the observed ranges and the sensor model.
    ///
    ///@param observedRanges (1 x N) array of ranges observed by a sensor.
    ///@param expectedRangesPerParticle (1 x N*numParticles) "pose-major" array of ranges computed
    /// by rangelib. This is usually the output of batchCalcRangeAngles.
    ///@throws std::invalid argument if expectedRangesPerParticle.size is not cleanly divisible by
    /// observedRanges.size, i.e. numParticles is not obtainable.
    ///@return (1 x numParticles) vector of particle weights
    VectorFloat batchEvalSensorModel(
        const Eigen::Ref<const VectorFloat> observedRanges,
        const Eigen::Ref<const VectorFloat> expectedRangesPerParticle) const
    {
        const size_t numObservations = observedRanges.size();
        if (expectedRangesPerParticle.size() % numObservations != 0) {
            throw std::invalid_argument(
                "Input Expected ranges must be evenly divisible by number of observations!");
        }
        const size_t numParticles = expectedRangesPerParticle.size() / numObservations;

        VectorFloat outs(numParticles);

        for (size_t particleIdx = 0; particleIdx < numParticles; ++particleIdx) {
            const auto startIdx = numObservations * particleIdx;
            const auto weight = evalSensorModel(
                observedRanges, expectedRangesPerParticle.segment(startIdx, numObservations));

            outs(particleIdx) = weight;
        }

        return outs;
    }

    ///@brief Function that calculates expected ranges for each particle and directly outputs
    /// particle weights
    ///
    ///@param poses (3 x N) array of poses.
    ///@param angles (1 x M) array of query angles.
    ///@param observedRanges (1 x M) array of ranges observed by a sensor.
    ///@throws std::invalid_argument if size of angles != size of observedRanges
    ///@return 1xN array of weights of each particle
    VectorFloat batchCalcParticleWeights(const Eigen::Ref<const PosesFloat> poses,
                                         const Eigen::Ref<const VectorFloat> angles,
                                         const Eigen::Ref<const VectorFloat> observedRanges) const
    {
        if (angles.size() != observedRanges.size()) {
            throw std::invalid_argument(
                "Number of query angles must correspond to number of observed ranges!");
        }

        VectorFloat ranges = batchCalcRangeAngles(poses, angles);
        VectorFloat vals = batchEvalSensorModel(observedRanges, ranges);

        return vals;
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
    inline Pose2Df_t mapCoordinates(float x_world, float y_world, float theta_world) const
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

    Eigen::MatrixXf _sensorModel;
};

}  // namespace ranges

#endif