#ifdef RANGELIB_CDDT_CAST_HPP_  // disables this file for now
#define RANGELIB_CDDT_CAST_HPP_

#include "rangelib/omap.hpp"
#include "rangelib/range_method.hpp"

#if _USE_LRU_CACHE
#include "lru_cache.h"
#endif

#define _TRACK_LUT_SIZE 0
#define _TRACK_COLLISION_INDEXES 0

#define _BINARY_SEARCH_THRESHOLD \
    64  // if there are more than this number of elements in the lut bin, use binary search

// fast optimized version
#define _USE_CACHED_TRIG 0
#define _USE_CACHED_CONSTANTS 1
#define _USE_FAST_ROUND 0
#define _NO_INLINE 0
#define _USE_LRU_CACHE 0
#define _LRU_CACHE_SIZE 1000000

// these defines are for yaml/JSON serialization
#define T1 "  "
#define T2 T1 T1
#define T3 T1 T1 T1
#define T4 T2 T2

#define J1 "  "
#define J2 J1 J1
#define J3 J1 J1 J1
#define J4 J2 J2

#define _CDDT_SHORT_DATATYPE 1

namespace ranges {

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
            float rotated_height =
                std::abs(map.width() * sinfangle) + std::abs(map.height() * cosfangle);
#else
            float rotated_height =
                std::abs(map.width() * sinf(angle)) + std::abs(map.height() * cosf(angle));
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

            // find the lowest corner, and determine the translation necessary to make them all
            // positive
            float min_corner_y =
                std::min(left_top_corner_y, std::min(right_top_corner_y, right_bottom_corner_y));
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

                        float half_lut_space_width =
                            (std::abs(sinangle) + std::abs(cosangle)) / 2.0;

                        float lut_space_center_x =
                            pixel_center.first * cosangle - pixel_center.second * sinangle;
                        float lut_space_center_y =
                            (pixel_center.first * sinangle + pixel_center.second * cosangle) +
                            lut_translations[a];

                        int upper_bin = lut_space_center_y + half_lut_space_width - _EPSILON;
                        int lower_bin = lut_space_center_y - half_lut_space_width + _EPSILON;

                        // the following is a quick hack to prevent problems in the cardinal
                        // directions where it has been known to see through walls if
                        // (std::fmod(angle, M_PI/2.0) < _EPSILON) { 	upper_bin++; 	lower_bin--;
                        // }

                        for (int i = lower_bin; i <= upper_bin; ++i)
                            compressed_lut[a][i].push_back(lut_space_center_x);

                        // std::cout << std::endl;
                        // std::cout << "angle: " << angle << std::endl;
                        // std::cout << "center: (" << pixel_center.first << ", " <<
                        // pixel_center.second << ")" << std::endl; std::cout << "new center: (" <<
                        // lut_space_center_x << ", " << lut_space_center_y << ")" << std::endl;
                        // std::cout << "bins:" << upper_bin << "    " << (int) lut_space_center_y
                        // << "   " << lower_bin << std::endl; std::cout << "width:" <<
                        // half_lut_space_width << std::endl; std::cout << "trans" <<
                        // lut_translations[a] << std::endl; std::cout << upper_bin << "   " <<
                        // lower_bin << "   " << lut_translations[a] << std::endl; std::cout <<
                        // lut_space_center_x << "  " << lut_space_center_y << std::endl;
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
                compressed_lut[a][i].erase(
                    unique(compressed_lut[a][i].begin(), compressed_lut[a][i].end()),
                    compressed_lut[a][i].end());
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
    void prune(float _maxRange)
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
                    if (high == -1) continue;
                    if (map.isOccupied(x, y)) continue;

                    // the furthest entry is behind the query point
                    // if ((*lut_bin)[high] + _maxRange < lut_space_x) return
                    // std::make_pair(_maxRange, _maxRange);
                    if ((*lut_bin)[high] < lut_space_x &&
                        lut_space_x - (*lut_bin)[high] < _maxRange) {
                        local_collision_table[angle_index][lut_index].insert(high);
                        // accum += 1;
                        continue;
                    }

                    int index;
                    if (high > _BINARY_SEARCH_THRESHOLD) {
                        // once the binary search terminates, the next greatest element is indicated
                        // by 'val'
                        index = std::lower_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) -
                                lut_bin->begin();
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
                    bool is_used = local_collision_table[a][lut_index].find(i) !=
                                   local_collision_table[a][lut_index].end();
                    if (is_used) pruned_bin.push_back(compressed_lut[a][lut_index][i]);
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
        if (theta < 0.0) theta += M_2PI;

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
        if (lut_index < 0 || lut_index >= compressed_lut[angle_index].size()) return _maxRange;
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
                cache.put(key, _maxRange);
#endif
                return _maxRange;
            }
            // the furthest entry is behind the query point
            if ((*lut_bin)[low] > lut_space_x) {
#if _USE_LRU_CACHE
                cache.put(key, _maxRange);
#endif
                return _maxRange;
            }
            if ((*lut_bin)[high] < lut_space_x) {
                float val = lut_space_x - (*lut_bin)[high];
#if _USE_LRU_CACHE
                cache.put(key, val);
#endif
                return val;
            }

            // the query point is on top of a occupied pixel
            // this call is here rather than at the beginning, because it is apparently more
            // efficient. I presume that this has to do with the previous two return statements
            if (map.isOccupied(x, y)) {
                return 0.0;
            }

            if (high > _BINARY_SEARCH_THRESHOLD) {
                int index = std::upper_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) -
                            lut_bin->begin();
                assert(index >
                       0);  // if index is 0, this will segfault. that should never happen, though.
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
                cache.put(key, _maxRange);
#endif
                return _maxRange;
            }
            // the furthest entry is behind the query point
            if ((*lut_bin)[high] < lut_space_x) {
#if _USE_LRU_CACHE
                cache.put(key, _maxRange);
#endif
                return _maxRange;
            }
            if ((*lut_bin)[low] > lut_space_x) {
                float val = (*lut_bin)[low] - lut_space_x;
#if _USE_LRU_CACHE
                cache.put(key, val);
#endif
                return val;
            }
            // the query point is on top of a occupied pixel
            // this call is here rather than at the beginning, because it is apparently more
            // efficient. I presume that this has to do with the previous two return statements
            if (map.isOccupied(x, y)) {
                return 0.0;
            }

            if (high > _BINARY_SEARCH_THRESHOLD) {
                // once the binary search terminates, the next greatest element is indicated by
                // 'val' float val = *std::lower_bound(lut_bin->begin(), lut_bin->end(),
                // lut_space_x);
                int index = std::upper_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) -
                            lut_bin->begin();
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
            if (high == -1) return std::make_pair(_maxRange, _maxRange);
            // the furthest entry is behind the query point and out of max range of the inverse
            // query if ((*lut_bin)[low] - _maxRange > lut_space_x) return std::make_pair(_maxRange,
            // _maxRange);
            if ((*lut_bin)[low] > lut_space_x)
                return std::make_pair(_maxRange,
                                      std::min(_maxRange, (*lut_bin)[low] - lut_space_x));
            if ((*lut_bin)[high] < lut_space_x)
                return std::make_pair(lut_space_x - (*lut_bin)[high], _maxRange);
            // the query point is on top of a occupied pixel
            // this call is here rather than at the beginning, because it is apparently more
            // efficient. I presume that this has to do with the previous two return statements
            if (map.isOccupied(x, y)) {
                return std::make_pair(0.0, 0.0);
            }

            float val;
            int index;
            if (high > _BINARY_SEARCH_THRESHOLD) {
                // once the binary search terminates, the next least element is indicated by 'val'
                // float val = *std::lower_bound(lut_bin->begin(), lut_bin->end(), lut_space_x);
                index = std::upper_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) -
                        lut_bin->begin() - 1;
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

                return std::make_pair(lut_space_x - val, _maxRange);
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
            if (high == -1) return std::make_pair(_maxRange, _maxRange);
            // the furthest entry is behind the query point
            // if ((*lut_bin)[high] + _maxRange < lut_space_x) return std::make_pair(_maxRange,
            // _maxRange);
            if ((*lut_bin)[high] < lut_space_x)
                return std::make_pair(_maxRange,
                                      std::min(_maxRange, lut_space_x - (*lut_bin)[high]));
            // TODO might need another early return case here
            // return std::make_pair(_maxRange, std::min(_maxRange, lut_space_x -
            // (*lut_bin)[high])); the query point is on top of a occupied pixel this call is here
            // rather than at the beginning, because it is apparently more efficient. I presume that
            // this has to do with the previous two return statements
            if (map.isOccupied(x, y)) {
                std::make_pair(0.0, 0.0);
            }

            float val;
            int index;
            if (high > _BINARY_SEARCH_THRESHOLD) {
                // once the binary search terminates, the next greatest element is indicated by
                // 'val' float val = *std::lower_bound(lut_bin->begin(), lut_bin->end(),
                // lut_space_x);
                index = std::lower_bound(lut_bin->begin(), lut_bin->end(), lut_space_x) -
                        lut_bin->begin();
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

                return std::make_pair(val - lut_space_x, _maxRange);
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
        (*ss) << T1 << "_maxRange: " << _maxRange << std::endl;
        (*ss) << T1 << "map: " << std::endl;
        (*ss) << T2
              << "# note: map data is width and then height (width is number of rows) transposed "
                 "from expectation:"
              << std::endl;
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
        (*ss) << J1 << "\"_maxRange\":" << _maxRange << "," << std::endl;
        (*ss) << J1 << "\"map\": {" << std::endl;
        // (*ss) << J2 << "# note: map data is width and then height (width is number of rows)
        // transposed from expectation:"  << std::endl;
        (*ss) << J2 << "\"path\": \"" << map.filename() << "\"," << std::endl;
        (*ss) << J2 << "\"width\": " << map.width() << "," << std::endl;
        (*ss) << J2 << "\"height\": " << map.height() << "," << std::endl;

        (*ss) << J2 << "\"data\": [";  // utils::serialize(map.grid()[0], ss);
        for (int i = 0; i < map.width(); ++i) {
            if (i > 0) (*ss) << ",";
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
                if (j > 0) (*ss) << ",";
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
}  // namespace ranges

#endif