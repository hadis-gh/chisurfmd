#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <random>
#include <algorithm>
#include "lettuce/IntersectionLineCircle.h"
#include "lettuce/Vec.h"
#include "lettuce/Circle.h"
#include "lettuce/CircleDistribution.h"

template<typename T>
std::pair<Circle<T>, T> findStopPoint(Circle<T> &startCircle, Vec<T> &direction, const Circle<T> &closestCircle) {
    const T a = direction.abs2();
    const T b = static_cast<T>(2.0) * direction * (startCircle.c - closestCircle.c);
    const T c = (startCircle.c - closestCircle.c).abs2() - 4 * startCircle.r * startCircle.r;

    const T discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return {{NAN, NAN}, NAN};
    }

    const T t1 = (-b + std::sqrt(discriminant)) / (2 * a);
    const T t2 = (-b - std::sqrt(discriminant)) / (2 * a);

    T t;
    if (t1 >= 0 && t2 >= 0) {
        t = std::min(t1, t2);
    } else if (t1 >= 0) {
        t = t1;
    } else if (t2 >= 0) {
        t = t2;
    } else {
        return {{NAN, NAN}, NAN};
    }
    return {{startCircle.c + t * direction}, t};
}

template<typename T>
Vec<T> shootToCenter(const Circle<T> &startCircle, const T &areaWidth) {
    return -startCircle.c + areaWidth / static_cast<T>(2);
}

template<typename T>
void writeCircles(T begin, T end, const std::string& fname) {
    std::ofstream output_file(fname);
    for (auto c = begin; c != end; ++c)
        output_file << c->c << ", " << c->r << "\n";
}

const float r = 1.0;
const float width = 50.0;

int main() {
    Circle<float> newCircle;
    Vec<float> direction;

    const Circle<float> target = {{{width / 2, width / 2}}, r};
    std::vector<Circle<float>> finalCircles;
    finalCircles.push_back(target);

    for (int j = 0; j < 10; ++j) {
        newCircle = startCircleRandom(r, width);
        direction = shootToCenter(newCircle, width);

        bool stuck = false;
        for (auto &c : finalCircles) {
            auto stoppedCircle = findStopPoint(newCircle, direction, c).first;
            if (!std::isnan(stoppedCircle.c[0])) {
                newCircle = stoppedCircle;
                stuck = true;
                break;
            }
        }
        if (stuck) {
            finalCircles.push_back(newCircle);
        }
    }

    writeCircles(finalCircles.begin(), finalCircles.end(), "circlecDLAtest.txt");

    return 0;
}
