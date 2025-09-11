#pragma once
#include <vector>
#include <iostream>
#include <cmath>

#include "config.h"

class Vector3 {
    private:
        real data_[3];
    
    public:
        real& x();
        real& y();
        real& z();
        real* data();

        const real& x() const;
        const real& y() const;
        const real& z() const;
        const real* data() const;

        // Constructors
        Vector3();
        Vector3(real x, real y, real z);
        explicit Vector3(const std::vector<real>& data);

        // Core operations
        Vector3 operator+(const Vector3& other) const;
        Vector3 operator-(const Vector3& other) const;
        Vector3 operator*(real scalar) const;         // Scalar multiplication
        Vector3 operator/(real scalar) const;         // Scalar division
        friend Vector3 operator*(real scalar, const Vector3& Vector3); // Commutation
        Vector3& operator+=(const Vector3& other);      // In-place addition
        Vector3& operator-=(const Vector3& other);      // In-place subtraction
        Vector3& operator*=(real scalar);   
        Vector3& operator/=(real scalar);   

        real dot(const Vector3& other) const;  // Element-wise multiplication
        Vector3 cross(const Vector3& other) const;  // Element-wise multiplication
 
        real mod() const; // module of the vector
        real mod2() const; // square module

        Vector3 norm() const; // Normalized vector

        // Periodic boundaries conditions, using a box centralized at the origin
        // In-place wrapping
        void wrap_PBC(const Vector3& box_size);

        // Supplementary PBC functions
        Vector3 PBCWrap(const Vector3& box_size) const;
        // Minimum image distance
        Vector3 PBCDisplacement(const Vector3& other, const Vector3& box_size) const;

        // Helper functions
        void print() const;
};