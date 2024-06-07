#include <iostream>
#include <cmath>
#include <vector>

std::vector<double> calculate_intersection(double px, double py, 
                                           double dx, double dy,
                                           double cx, double cy, double r) 
{   
    double a = dx * dx + dy * dy;
    double b = 2 * (dx * (px - cx) + dy * (py - cy));
    double c = (px - cx) * (px - cx) + (py - cy) * (py - cy) - r * r;

    double discriminant = b * b - 4 * a * c;

    if (discriminant < 0) {
        return {};
    }
    double t1 = (-b + std::sqrt(discriminant)) / (2 * a);
    double t2 = (-b - std::sqrt(discriminant)) / (2 * a);

    double t;
    if (t1 >= 0 && t2 >= 0) {
        t = std::min(t1, t2);
    } else if (t1 >= 0) {
        t = t1;
    } else if (t2 >= 0) {
        t = t2;
    } else {
        return {};
    }

    double ix = px + t * dx;
    double iy = py + t * dy;

    double distance = t * std::sqrt(dx * dx + dy * dy);

    return {ix, iy, distance};
}

int main() {
    double px = 0.0, py = 0.0;
    double dx = 1.0, dy = 1.0;
    double cx = 3.0, cy = 3.0, r = 1.0;

    auto result = calculate_intersection(px, py, dx, dy, cx, cy, r);

    if (!result.empty()) {
        std::cout << "Intersection point: (x= " << result[0] << ", y= " << result[1] << ")\n";
        std::cout << "Distance: " << result[2] << "\n";
    } else {
        std::cout << "No intersection.\n";
    }

    return 0;
}