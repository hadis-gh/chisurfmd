#include <iostream>
#include <lettuce/Intersection_line_circle.h>


int main() {
    const Vec<float> point ({0.0, 0.0});
    const Vec<float> direction ({1.0, 0.0});     //= {1.0, 1.0}

    const Circle<float> c1 = {{{3.0, -1.0}}, 1.0};
    const Circle<float> c2 = {{{0.0, 3.0}}, 1.0};
    const Circle<float> c3 = {{{0.0, 3.0}}, 1.0};

    const auto distance = calculateIntersection(point, direction, c1);
    if (std::isnan(distance)) {
        std::cout << "No intersection.\n";
    } else {
        std::cout << "Distance: " << distance << "\n";
        std::cout << "Intersection point: " << (point + distance*direction) << "\n";
    }

    return 0;
}