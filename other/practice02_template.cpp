#include <iostream>
#include <array>

template<typename T , int D = 2> 
class Vec2D {
    std::array<T, D> elements;

public:

    Vec2D() : elements() {}
    Vec2D(T x, T y) : elements{x, y} {}

    T& operator[](int index) {
        return elements[index];
    }

    const T& operator[](int index) const {
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

    auto operator* (const Vec2D& other) const {
        return (elements[0] * other.elements[0] + elements[1] * other.elements[1]);
    }

    auto operator^(const Vec2D& other) const {
        return (elements[0] * other.elements[1] - elements[1] * other.elements[0]);
    }

    Vec2D operator* (const double& number)const{
        return Vec2D (number*elements[0], number*elements[1]);
    }

    friend Vec2D operator* (const double& number, const Vec2D<T>& vec);

};

template <typename T>
Vec2D<T> operator* (const double& number, const Vec2D<T>& vec){
    return Vec2D<T> (number*vec.elements[0], number*vec.elements[1]);
}
// can't I define * inside the class some how I could calculate 2 * v1 ? I mean before elements of vector

std::ostream& operator<<(std::ostream& COUT, const Vec2D<double>& ref) {
    COUT << "(" << ref[0] << " , " << ref[1] << " )";
    return COUT;
}

int main() {
    Vec2D<double> v1(1, 0);     // why this Vec2D v1(1, 0); does not work and we shouold use <>
    Vec2D<double> v2(1, -1);

    std::cout << -v1 << std::endl;
    std::cout << (v1 + (v2 * 4)) << std::endl;
    std::cout << v1 - v2 << std::endl;
    std::cout << v1 * v2 << std::endl;
    std::cout << (v1 ^ v2) << std::endl;
    std::cout << v1 * 2 << std::endl;
    std::cout << 2 * v1 << std::endl;
    
    return 0;
}