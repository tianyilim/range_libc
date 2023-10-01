#ifndef RANGELIB_RAY_MARCHING_HPP_
#define RANGELIB_RAY_MARCHING_HPP_

#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"

namespace ranges {

class RayMarching : public RangeMethod {
   public:
    RayMarching(OMap m, float mr, float distThreshold = 0.0, float stepCoeff = 0.999)
        : RangeMethod(m, mr), _distThreshold{distThreshold}, _stepCoeff{stepCoeff}
    {
    }

    float calc_range(const float x, const float y, const float heading) const override
    {
        // Catch illegal array access
        if (x < 0 || x >= _distTransform.width() || y < 0 || y >= _distTransform.height()) {
            return _maxRange;
        }

        float ray_direction_x = cosf(heading);
        float ray_direction_y = sinf(heading);

        int px, py;

        float t = 0.0;
        while (t < _maxRange) {
            px = x + ray_direction_x * t;
            py = y + ray_direction_y * t;

            if (px >= (int)_distTransform.width() || px < 0 || py < 0 ||
                py >= (int)_distTransform.height()) {
                return _maxRange;
            }

            float d = _distTransform.signedDistanceValue(px, py);

            if (d <= _distThreshold) {
                float xd = px - x;
                float yd = py - y;
                return sqrtf(xd * xd + yd * yd);
            }

            t += std::max<float>(d * _stepCoeff, 1.0);
        }

        return _maxRange;
    }

    int memory() const override { return _distTransform.memory(); }

   protected:
    const float _distThreshold;
    const float _stepCoeff;
};

}  // namespace ranges

#endif