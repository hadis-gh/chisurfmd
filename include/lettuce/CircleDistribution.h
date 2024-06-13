#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <random>
#include "Circle.h"

template <typename T>
bool has_overlap (const Circle<T> &c1,const Circle<T> &c2){
    T distance2 = (c1.c - c2.c).abs2();
    const auto radiuses (c1.r+c2.r);
    return distance2 < radiuses*radiuses;
}

template <typename T>
std::vector <Circle<T>> distCircles (const int &circles_number, const T &L, const T &radius, std::mt19937 &gen){
    std::vector <Circle<T>> circles;
    std::uniform_real_distribution<> dis(radius, L-radius);
    
    while (circles.size() < circles_number){
        Circle<T> new_circle;
        new_circle.c[0] = dis(gen);
        new_circle.c[1] = dis(gen);
        new_circle.r = radius;

        bool overlap = false;
        
        for (auto circle: circles){
            if(has_overlap(new_circle, circle)){
                overlap = true;
                break;
            }
        }
        if (!overlap){
            circles.push_back(new_circle);
        }
    }
    return circles;
}