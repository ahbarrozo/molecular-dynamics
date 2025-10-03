#include "vector3.h"

Vector3::Vector3() : data_{0.0f, 0.0f, 0.0f} {}  // Creates an empty vector at origin
Vector3::Vector3(real x, real y, real z) : data_{x, y, z} {}
Vector3::Vector3(const std::vector<real>& data) {
    if (data.size() != 3) {
        throw std::invalid_argument("Vector3 requires exactly 3 elements");
    }
    data_[0] = data[0];
    data_[1] = data[1];
    data_[2] = data[2];
}

real& Vector3::x() { return data_[0]; }
real& Vector3::y() { return data_[1]; }
real& Vector3::z() { return data_[2]; }
real* Vector3::data() { return data_; }

const real& Vector3::x() const { return data_[0]; }
const real& Vector3::y() const { return data_[1]; }
const real& Vector3::z() const { return data_[2]; }
const real* Vector3::data() const { return data_; }


// Addition operator
Vector3 Vector3::operator+(const Vector3& other) const {
    return Vector3(
        data_[0] + other.data_[0],
        data_[1] + other.data_[1],
        data_[2] + other.data_[2]
    );
}

// In-place addition 
Vector3& Vector3::operator+=(const Vector3& other) {
    data_[0] += other.data_[0];
    data_[1] += other.data_[1];
    data_[2] += other.data_[2];
    return *this;
}

// Subtraction operator
Vector3 Vector3::operator-(const Vector3& other) const {
    return Vector3(
        data_[0] - other.data_[0],
        data_[1] - other.data_[1],
        data_[2] - other.data_[2]
    );
}

// In-place subtraction
Vector3& Vector3::operator-=(const Vector3& other) {
    data_[0] -= other.data_[0];
    data_[1] -= other.data_[1];
    data_[2] -= other.data_[2];
    return *this;
}

// Scalar product operator
Vector3 Vector3::operator*(real scalar) const {
    return Vector3(
        scalar * data_[0],
        scalar * data_[1],
        scalar * data_[2]
    );
}

// Scalar division operator
Vector3 Vector3::operator/(real scalar) const {
    if (std::abs(scalar) < 1e-10f) {
        throw std::runtime_error("Division by zero");
    }

    return Vector3(
        data_[0] / scalar,
        data_[1] / scalar,
        data_[2] / scalar
    );
}

// In-place scalar product
Vector3& Vector3::operator*=(real scalar) {
    data_[0] *= scalar;
    data_[1] *= scalar;
    data_[2] *= scalar;
    return *this;
}

// In-place scalar division
Vector3& Vector3::operator/=(real scalar) {
    data_[0] /= scalar;
    data_[1] /= scalar;
    data_[2] /= scalar;
    return *this;
}

// Commuted scalar multiplication: scalar * tensor (friend function)
Vector3 operator*(real scalar, const Vector3& vector3) {
    return vector3 * scalar;  // Delegate to the member function
}

// Commuted scalar multiplication: scalar * tensor (friend function)
Vector3 operator/(real scalar, const Vector3& vector3) {
    return scalar / vector3;  // Delegate to the member function
}

real Vector3::dot(const Vector3& other) const {
    return (data_[0] * other.x() + data_[1] * other.y() + data_[2] * other.z());
}

Vector3 Vector3::cross(const Vector3& other) const {
    return Vector3(
        data_[1] * other.data_[2] - data_[2] * other.data_[1],
        data_[2] * other.data_[0] - data_[0] * other.data_[2],
        data_[0] * other.data_[1] - data_[1] * other.data_[0]
    );
}

real Vector3::mod() const {
    return sqrt(dot(*this));  // Reuse dot product
}

real Vector3::mod2() const {
    return dot(*this);  // Avoid sqrt when possible!
}

Vector3 Vector3::norm() const {
    real module = mod();
    if (module > 0.0f) {
        return (*this) * (1.0f / module);
    }
    return *this;  // or throw exception for zero vector
}

// Wrapped vector in a PBC box centralized at 0
void Vector3::wrap_PBC(const Vector3& box_size) {
        data_[0] -= box_size.x() * round(data_[0] / box_size.x());
        data_[1] -= box_size.y() * round(data_[1] / box_size.y());
        data_[2] -= box_size.z() * round(data_[2] / box_size.z());
}

// Wrapped vector in a PBC box centralized at 0
Vector3 Vector3::PBCWrap(const Vector3& box_size) const {
    return Vector3(
        data_[0] - box_size.x() * round(data_[0] / box_size.x()),
        data_[1] - box_size.y() * round(data_[1] / box_size.y()),
        data_[2] - box_size.z() * round(data_[2] / box_size.z())
    );
}

// Minimum image distance
Vector3 Vector3::PBCDisplacement(const Vector3& other, const Vector3& box_size) const {
    Vector3 diff = other - *this;
    return diff.PBCWrap(box_size);
}

void Vector3::print() const {
    std::cout << "Vector3 coordinates: ( " ;

    for (size_t i = 0; i < 3; ++i) {
        std::cout << data_[i] << " ";
    }

    std::cout << ")" << std::endl;
}