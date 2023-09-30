#include "RangeUtils.h"

namespace utils {

unsigned long xorshf96(void)
{  // period 2^96-1
    // return std::rand() / RAND_MAX;
    unsigned long t;
    x ^= x << 16;
    x ^= x >> 5;
    x ^= x << 1;

    t = x;
    x = y;
    y = z;
    z = t ^ x ^ y;

    return z;
}

float rgb2gray(float r, float g, float b) { return 0.229 * r + 0.587 * g + 0.114 * b; }

int randrange(int min, int max) { return min + (rand() % (int)(max - min + 1)); }

std::vector<std::pair<int, int>> outline(int x, int y, bool use_corners)
{
    std::vector<std::pair<int, int>> corners;

    corners.push_back(std::make_pair(x + 1, y));
    corners.push_back(std::make_pair(x - 1, y));
    corners.push_back(std::make_pair(x, y + 1));
    corners.push_back(std::make_pair(x, y - 1));

    if (use_corners) {
        corners.push_back(std::make_pair(x + 1, y + 1));
        corners.push_back(std::make_pair(x - 1, y + 1));
        corners.push_back(std::make_pair(x + 1, y - 1));
        corners.push_back(std::make_pair(x - 1, y - 1));
    }

    return corners;
}

bool has(std::string substring, std::string str)
{
    return str.find(substring) != std::string::npos;
}

bool has(std::string val, std::vector<std::string> vstr)
{
    return std::find(vstr.begin(), vstr.end(), val) != vstr.end();
}

std::vector<std::string> split(std::string in, char delim)
{
    std::vector<std::string> result;
    std::stringstream ss(in);
    while (ss.good()) {
        std::string substr;
        std::getline(ss, substr, delim);
        result.push_back(substr);
    }
    return result;
}

double norminv(double q)
{
    if (q == .5) return 0;

    q = 1.0 - q;

    double p = (q > 0.0 && q < 0.5) ? q : (1.0 - q);
    double t = sqrt(log(1.0 / pow(p, 2.0)));

    double c0 = 2.515517;
    double c1 = 0.802853;
    double c2 = 0.010328;

    double d1 = 1.432788;
    double d2 = 0.189269;
    double d3 = 0.001308;

    double x =
        t - (c0 + c1 * t + c2 * pow(t, 2.0)) / (1.0 + d1 * t + d2 * pow(t, 2.0) + d3 * pow(t, 3.0));

    if (q > .5) x *= -1.0;

    return x;
}

NonReplacementSampler::NonReplacementSampler()
{
    rand = std::uniform_real_distribution<double>(0.0, 1.0);
    generator.seed(clock());
}

void NonReplacementSampler::sample(int populationSize, int sampleSize, std::vector<int> &samples)
{
    int t = 0;  // total input records dealt with
    int m = 0;  // number of items selected so far
    double u;

    while (m < sampleSize) {
        // u = rand(generator);
        u = std::rand() / (float)RAND_MAX;
        if ((populationSize - t) * u >= sampleSize - m)
            t++;
        else {
            samples.push_back(t);
            t++;
            m++;
        }
    }
    // samples.push_back(1);
}

FastRand::FastRand() : FastRand(10000){};
FastRand::FastRand(int n) : cache_size(n)
{
    populate();
    repopulate_threshold = 1.0 / cache_size;
}

float FastRand::rand()
{
    // return std::rand() / (float)RAND_MAX;
    // float v = cache[i];
    if (i++ > cache_size - 1) i = 0;
    // if (v < repopulate_threshold) populate();
    return cache[i];
}

void FastRand::populate()
{
    // cache.empty();
    // for (int i = 0; i < cache_size; ++i) cache.push_back(std::rand() /
    // (float)RAND_MAX);
    for (int i = 0; i < cache_size; ++i) cache[i] = std::rand() / (float)RAND_MAX;
}

void serialize(std::vector<bool> &vals, std::stringstream *ss)
{
    if (vals.size() == 0) {
        (*ss) << "[]";
        return;
    }
    (*ss) << "[" << vals[0];
    for (size_t i = 1; i < vals.size(); ++i) {
        (*ss) << "," << vals[i];
    }
    (*ss) << "]";
}

void serialize(std::vector<float> &vals, std::stringstream *ss)
{
    if (vals.size() == 0) {
        (*ss) << "[]";
        return;
    }
    (*ss) << "[" << vals[0];
    for (size_t i = 1; i < vals.size(); ++i) {
        (*ss) << "," << vals[i];
    }
    (*ss) << "]";
}

std::string serialize(std::vector<float> &vals)
{
    std::stringstream ss;
    serialize(vals, &ss);
    return ss.str();
}

}  // namespace utils