#ifndef RANGELIB_RAY_MARCHING_HPP_
#define RANGELIB_RAY_MARCHING_HPP_

#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"

namespace ranges {

class RayMarching : public RangeMethod {
   public:
    RayMarching(OMap m, float mr) : RangeMethod(m, mr) { _distImage = DistanceTransform(m); }

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

}  // namespace ranges

#endif