#include <iostream>
#include <array>

class Vec2D {
    std::array<double, 2> elements;

public:

    Vec2D() : elements() {}
    Vec2D(double x, double y) : elements{x, y} {}

    double& operator[](int index) {
        return elements[index];
    }

    double operator[](int index) const {
        return elements[index];
    }

    Vec2D operator+(const Vec2D& other) const {
        return Vec2D(elements[0] + other.elements[0], elements[1] + other.elements[1]);
    }

    Vec2D operator-(const Vec2D& other) const {
        return Vec2D(elements[0] - other.elements[0], elements[1] - other.elements[1]);
    }

    Vec2D operator-() const {
        return Vec2D(-elements[0], -elements[1]);
    }

    double operator* (const Vec2D& other) const {
        return (elements[0] * other.elements[0] + elements[1] * other.elements[1]);
    }

    double operator^(const Vec2D& other) const {
        return (elements[0] * other.elements[1] - elements[1] * other.elements[0]);
    }

    Vec2D operator* (const double& number)const{
        return Vec2D (number*elements[0], number*elements[1]);
    }

    Vec2D& operator*= (const double& number){
        elements[0] *= number;
        elements[1] *= number;
        return *this;
    }

    friend Vec2D operator* (const double& number, const Vec2D& vec);

};

Vec2D operator* (const double& number, const Vec2D& vec){
    return Vec2D (number*vec.elements[0], number*vec.elements[1]);
}

std::ostream& operator<<(std::ostream& COUT, const Vec2D& ref) {
    COUT << "(" << ref[0] << " , " << ref[1] << " )";
    return COUT;
}

int main() {
    Vec2D v1(1, 0);
    Vec2D v2(1, -1);

    std::cout << -v1 << std::endl;
    std::cout << v1 + v2 << std::endl;
    std::cout << v1 - v2 << std::endl;
    std::cout << v1 * v2 << std::endl;
    std::cout << (v1 ^ v2) << std::endl;
    std::cout << v1 * 2 << std::endl;
    std::cout << 2 * v1 << std::endl;
    

    return 0;
}