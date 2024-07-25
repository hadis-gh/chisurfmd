#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <cmath>
#include <algorithm>
#include "Vec.h"
#include "Circle.h"

// template<typename T>
// std::pair<Particle<T>, T> findStopPoint (const Particle<T> &startParticle, const T &radius, Vec<T> &direction, const Particle<T> &closestParticle){

//     const T a = direction.abs2();
//     const T b = static_cast<T>(2.0) * direction*(startParticle.r - closestParticle.r);
//     const T c = (startParticle.r - startParticle.r).abs2() - 4 * radius * radius;

//     const T discriminant = b * b - 4 * a * c;
//     if (discriminant < 0) {
//         return {{NAN, NAN}, NAN};
//     }

//     const T t1 = (-b + std::sqrt(discriminant)) / (2 * a);
//     const T t2 = (-b - std::sqrt(discriminant)) / (2 * a);

//     T t;
//     if (t1 >= 0 && t2 >= 0) {
//         t = std::min(t1, t2);
//     } else if (t1 >= 0) {
//         t = t1;
//     } else if (t2 >= 0) {
//         t = t2;
//     } else {
//         return {{NAN, NAN}, NAN};
//     }
//     return {{startParticle.r + t * direction}, t};
// }

// template<typename T>
// Particle<T> findStopPointAll (const Particle<T> &startParticle, const T& radius, Vec<T> &direction, const std::vector<Particle<T>> &particles){
//     std::vector<std::pair<Particle<T>, T>> pairs;
//     for (const auto &p: particles){
//         pairs.push_back(findStopPoint(startParticle, radius, direction, p));
//     }

//     auto it = std::min_element(pairs.begin(), pairs.end(),
//         [](const std::pair<Particle<T>, T>& a, const std::pair<Particle<T>, T>& b) {
//             if (!std::isnan(a.second) && !std::isnan(b.second)) {
//                 return a.second < b.second;
//             }
//             else if (std::isnan(a.second)) {
//                 return false;
//             }
//             else {
//                 return true;
//             }
//         });

//     return it->first;
// }