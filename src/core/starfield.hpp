#pragma once
#include <cmath>
#include <random>
#include <vector>

// Simple procedural star field
struct StarField {
    struct Star {
        double theta;      // Polar angle (0 to π)
        double phi;        // Azimuthal angle (0 to 2π)
        double brightness; // 0 to 1
    };
    
    std::vector<Star> stars;
    
    // Generate random star field
    StarField(int num_stars = 5000, unsigned int seed = 42) {
        std::mt19937 rng(seed);
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        
        stars.reserve(num_stars);
        
        for (int i = 0; i < num_stars; ++i) {
            Star s;
            // Uniform distribution on sphere
            s.theta = std::acos(2.0 * uniform(rng) - 1.0);
            s.phi = 2.0 * M_PI * uniform(rng);
            s.brightness = std::pow(uniform(rng), 2.0); // Bias toward dimmer stars
            stars.push_back(s);
        }
    }
    
    // Get color at direction (theta, phi) on celestial sphere
    struct Color {
        double r, g, b;
    };
    
    Color sample(double theta, double phi) const {
        // Background sky color (very dark blue)
        Color bg = {0.01, 0.01, 0.05};
        
        // Check if we're near any star
        const double star_radius = 0.01; // Angular size in radians
        
        for (const auto& star : stars) {
            // Angular distance to star
            double cos_dist = std::sin(theta) * std::sin(star.theta) * 
                             std::cos(phi - star.phi) +
                             std::cos(theta) * std::cos(star.theta);
            double angular_dist = std::acos(std::clamp(cos_dist, -1.0, 1.0));
            
            if (angular_dist < star_radius) {
                // Smooth falloff
                double falloff = 1.0 - angular_dist / star_radius;
                double intensity = star.brightness * falloff * falloff;
                
                // Star color (slightly warm white)
                return {
                    0.9 * intensity,
                    0.85 * intensity,
                    0.7 * intensity
                };
            }
        }
        
        return bg;
    }
};

// Simple accretion disk model
struct AccretionDisk {
    double inner_radius;  // In Schwarzschild radii
    double outer_radius;
    double thickness;     // Vertical thickness
    
    AccretionDisk(double r_inner = 3.0, double r_outer = 10.0, double thick = 0.3)
        : inner_radius(r_inner), outer_radius(r_outer), thickness(thick) {}
    
    // Check if a point (r, θ) intersects the disk
    // Returns: {intersects, emission_color}
    struct Emission {
        bool hit;
        double r, g, b;
    };
    
    Emission sample(double r, double theta) const {
        // Check if in equatorial plane
        double z = r * std::cos(theta);
        if (std::abs(z) > thickness) {
            return {false, 0, 0, 0};
        }
        
        // Check if in disk radial range
        double r_cyl = r * std::sin(theta); // Cylindrical radius
        if (r_cyl < inner_radius || r_cyl > outer_radius) {
            return {false, 0, 0, 0};
        }
        
        // Temperature decreases with radius (simple model)
        // Hotter near inner edge (bluer), cooler at outer edge (redder)
        double temp_factor = (outer_radius - r_cyl) / (outer_radius - inner_radius);
        
        // Black body approximation
        double red   = 0.5 + 0.5 * (1.0 - temp_factor);
        double green = 0.4 + 0.4 * temp_factor;
        double blue  = 0.3 + 0.7 * temp_factor;
        
        // Brightness decreases with radius
        double brightness = temp_factor * 2.0;
        
        return {
            true,
            red * brightness,
            green * brightness,
            blue * brightness
        };
    }
};