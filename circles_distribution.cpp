#include <iostream>
#include <vector>
#include <array>
#include <cmath>
#include <cstdlib>
#include <ctime>

struct circle {
    double x;
    double y;
};

double random_double(const double &min, const double &max) {
    srand((time(0)));
    return min + rand() % int(max) +1;
}

bool has_overlap (const circle &new_circle,const circle &old_circle, const int &radius){
    double distance = std::sqrt ((new_circle.x -old_circle.x)*(new_circle.x -old_circle.x) + (new_circle.y -old_circle.y)*(new_circle.y -old_circle.y));
    if (distance < 2*double(radius))
        return true;
    return false;
}

std::vector <circle> dist_circles (const int &circles_number, const double &L, const double &radius){
    std::vector <circle> circles;
    
    while (circles.size() < circles_number){
    circle new_circle;
    new_circle.x = random_double(radius, L-radius);
    new_circle.y = random_double(radius, L-radius);
    bool overlap = false;
    
    for (auto circle: circles){
        if(!has_overlap(new_circle, circle, radius)){
            circles.push_back(new_circle);
        }else{
            overlap = true;
            break;
        }
    }
    }
    return circles;
}

int main() {
    int circles_number = 20;
    double length = 100.0;
    double radius = 1.0;

    std::vector <circle> circles = dist_circles (circles_number, length, radius);
    for (auto c: circles)
        std::cout << c.x << ", "<< c.y <<"\n";
    return 0;
}