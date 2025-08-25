#pragma once

template<typename T>
struct Species
{
    T mass;
    T momentOfInertia;
    T radius;

    Species(T mass, T momentOfInertia, T radius)
        : mass(mass), momentOfInertia(momentOfInertia), radius(radius)
    {}
};