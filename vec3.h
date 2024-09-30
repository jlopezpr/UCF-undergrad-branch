#ifndef VEC3_H
#define VEC3_H

#include <cmath>
#include <iostream>

struct alignas(32) vec3 {
    double x, y, z;

    vec3(double x = 0, double y = 0, double z = 0) : x(x), y(y), z(z) {}

     // Operador de igualdad
    bool operator==(const vec3& other) const {
        return x == other.x && y == other.y && z == other.z;
    }

     // Operador de desigualdad
    bool operator!=(const vec3& other) const {
        return !(*this == other);
    }

    double& operator[](int index) {
        return *(&x + index);
    }

    const double& operator[](int index) const {
        return *(&x + index);
    }

    vec3 operator+(const vec3& other) const {
        return vec3(x + other.x, y + other.y, z + other.z);
    }

    vec3 operator-(const vec3& other) const {
        return vec3(x - other.x, y - other.y, z - other.z);
    }

    vec3 operator*(double scalar) const {
        return vec3(x * scalar, y * scalar, z * scalar);
    }

    // Agrega multiplicación por un escalar al principio
    friend vec3 operator*(double scalar, const vec3& vec) {
        return vec * scalar;
    }

    double dot(const vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    double magnitude() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    double squaremag() const {
        return x * x + y * y + z * z;
    }

    vec3 normalize() const {
        double mag = magnitude();
        if (mag > 0) {
            return vec3(x / mag, y / mag, z / mag);
        }
        return *this;
    }

    vec3 cross(const vec3& other) const {
        return vec3(
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        );
    }

    vec3 hadda(const vec3& other) const {
        return vec3(
            x * other.x,
            y * other.y,
            z * other.z
        );
    }

    friend std::ostream& operator<<(std::ostream& os, const vec3& vec) {
        os << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
        return os;
    }
};

#endif // VEC3_H

