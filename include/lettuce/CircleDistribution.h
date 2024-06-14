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
    std::uniform_real_distribution<> dis(0, L);
    
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

template<typename T>
std::vector<Circle<T>> distCirclesPBC (const int &circles_number, const T &L, const T &radius, std::mt19937 &gen){
    std::vector<Circle<T>> repCircles = distCircles (circles_number, L, radius, gen);
    
    for (const auto &circle : repCircles) {
        Vec<T> c = circle.c;
        T r = circle.r;
        // left, right, bottom, up
        if (c[0] - r < 0)
            repCircles.emplace_back(Vec<T>({c[0] + L, c[1]}), r);
        if (c[0] + r >= L) 
            repCircles.emplace_back(Vec<T>({c[0] - L, c[1]}), r);
        if (c[1] - r < 0) 
            repCircles.emplace_back(Vec<T>({c[0], c[1] + L}), r);
        if (c[1] + r >= L)
            repCircles.emplace_back(Vec<T>({c[0], c[1] - L}), r);
        //corners
        if (c[0] - r < 0 && c[1] - r < 0)
            repCircles.emplace_back(Vec<T>({c[0] + L, c[1] + L}), r);
        if (c[0] + r >= L && c[1] - r < 0) 
            repCircles.emplace_back(Vec<T>({c[0] - L, c[1] + L}), r);
        if (c[0] - r < 0 && c[1] + r >= L)
            repCircles.emplace_back(Vec<T>({c[0] + L, c[1] - L}), r);
        if (c[0] + r >= L && c[1] + r >= L)
            repCircles.emplace_back(Vec<T>({c[0] - L, c[1] - L}), r);
    }
    return repCircles;
}