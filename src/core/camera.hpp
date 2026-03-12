#pragma once
#include "geodesic.hpp"
#include <cmath>
#include <array>

// Camera in spacetime with position, orientation, and field of view
struct Camera {
    // Position in spacetime coordinates (t, r, θ, φ) or (t, x, y, z)
    std::array<double, 4> position;
    
    // Orientation: spherical angles (theta, phi) for look direction
    double look_theta;  // Polar angle (0 to π)
    double look_phi;    // Azimuthal angle (0 to 2π)
    
    // Camera parameters
    double fov;         // Field of view in radians (e.g., π/3 for 60 degrees)
    int width;          // Image width in pixels
    int height;         // Image height in pixels
    
    // Velocity (for moving camera, used in Doppler calculations)
    std::array<double, 4> velocity;
    
    // Constructor
    Camera(const std::array<double, 4>& pos, 
           double theta, double phi,
           double field_of_view,
           int w, int h)
        : position(pos)
        , look_theta(theta)
        , look_phi(phi)
        , fov(field_of_view)
        , width(w)
        , height(h)
        , velocity({1.0, 0.0, 0.0, 0.0})  // Stationary by default
    {}
    
    // Get camera basis vectors in local frame
    // Returns: (right, up, forward) unit vectors
    struct Basis {
        std::array<double, 3> right;
        std::array<double, 3> up;
        std::array<double, 3> forward;
    };
    
    Basis get_basis() const {
        Basis basis;
        
        // Forward direction (look direction)
        double sin_theta = std::sin(look_theta);
        double cos_theta = std::cos(look_theta);
        double sin_phi = std::sin(look_phi);
        double cos_phi = std::cos(look_phi);
        
        basis.forward = {
            sin_theta * cos_phi,
            sin_theta * sin_phi,
            cos_theta
        };
        
        // Up direction (perpendicular to forward, points toward north pole)
        basis.up = {
            cos_theta * cos_phi,
            cos_theta * sin_phi,
            -sin_theta
        };
        
        // Right direction (cross product: forward × up)
        basis.right = {
            -sin_phi,
            cos_phi,
            0.0
        };
        
        return basis;
    }
    
    // Convert pixel coordinates to ray direction
    // Returns direction in camera's local frame (unit vector)
    std::array<double, 3> pixel_to_direction(int px, int py) const {
        // Normalized device coordinates: [-1, 1] × [-1, 1]
        // (0, 0) is top-left, center is (width/2, height/2)
        double aspect = static_cast<double>(width) / height;
        double u = (2.0 * px / width - 1.0) * aspect;
        double v = (1.0 - 2.0 * py / height);  // Flip y axis
        
        // Scale by FOV
        double tan_half_fov = std::tan(fov / 2.0);
        u *= tan_half_fov;
        v *= tan_half_fov;
        
        // Direction in camera space (right, up, forward)
        Basis basis = get_basis();
        
        std::array<double, 3> dir = {
            u * basis.right[0] + v * basis.up[0] + basis.forward[0],
            u * basis.right[1] + v * basis.up[1] + basis.forward[1],
            u * basis.right[2] + v * basis.up[2] + basis.forward[2]
        };
        
        // Normalize
        double len = std::sqrt(dir[0]*dir[0] + dir[1]*dir[1] + dir[2]*dir[2]);
        if (len > 1e-10) {
            dir[0] /= len;
            dir[1] /= len;
            dir[2] /= len;
        }
        
        return dir;
    }
};

// Ray generator: converts camera + pixel → initial geodesic state
struct RayGenerator {
    const Camera& camera;
    const Metric& metric;
    
    RayGenerator(const Camera& cam, const Metric& met)
        : camera(cam), metric(met) {}
    
    // Generate initial geodesic state for a given pixel
    // This is coordinate-system specific (needs to be specialized per metric)
    GeodesicState generate_ray(int px, int py) const {
        GeodesicState state;
        
        // Get ray direction in camera's local frame
        std::array<double, 3> dir = camera.pixel_to_direction(px, py);
        
        // Initial position = camera position
        state.x = camera.position;
        
        // Initial 4-velocity (needs coordinate transformation)
        // For now, simple flat-space approximation
        // TODO: Proper transformation using metric at camera position
        //This prrtion not to be used 
        // Normalize to null geodesic (light ray): g_μν u^μ u^ν = 0
        // For flat space: u^0 = 1, u^i = direction
        state.u = {1.0, dir[0], dir[1], dir[2]};
        
        return state;
    }
};

// Schwarzschild-specific ray generator
struct SchwarzschildRayGenerator {
    const Camera& camera;
    const Metric& metric;
    
    SchwarzschildRayGenerator(const Camera& cam, const Metric& met)
        : camera(cam), metric(met) {}
    
    GeodesicState generate_ray(int px, int py) const {
        GeodesicState state;
        
        // Get ray direction in camera's local frame
        std::array<double, 3> dir = camera.pixel_to_direction(px, py);
        
        // Camera position in Schwarzschild coordinates (t, r, θ, φ)
        state.x = camera.position;
        double r = state.x[1];
        double theta = state.x[2];
        
        // Schwarzschild metric factor
        double f = 1.0 - 2.0 / r;  // Assumes M=1
        
        if (f <= 0.0) {
            // Inside horizon, return zero velocity
            state.u = {0.0, 0.0, 0.0, 0.0};
            return state;
        }
        
        // Transform direction to Schwarzschild coordinate velocities
        // This is simplified: assumes camera is stationary in Schwarzschild coordinates
        // and direction is given in local orthonormal frame
        
        // For light ray: E = 1 (energy), then construct 4-velocity
        double E = 1.0;
        state.u[0] = E / f;  // u^t
        
        // Spatial components (approximate, needs proper basis transformation)
        // For distant observer: u^r ≈ dir[0], u^θ/r ≈ dir[1], u^φ/(r sinθ) ≈ dir[2]
        state.u[1] = dir[0] * E;  // u^r
        state.u[2] = dir[1] * E / r;  // u^θ
        state.u[3] = dir[2] * E / (r * std::sin(theta));  // u^φ
        
        return state;
    }
};