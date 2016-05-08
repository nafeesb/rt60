//
//  main.cpp
//  rt60
//
//  Created by Nafees Bin Zafar on 5/8/16.
//  Copyright (c) 2016 Nafees Bin Zafar. All rights reserved.
//

#include <iostream>
#include <vector>
#include <list>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>

using namespace std;

class Vec3
{
public:
    typedef float value_type;
    typedef value_type value_t;
    value_type data[3];

    /// @name ctors
    /// @{
    Vec3() :
        data{0,0,0}
    {}
    
    Vec3( const value_t x, const value_t y, const value_t z) :
        data{x,y,z}
    {}
    
    Vec3( const float v ) :
        data{v,v,v}
    {}
    /// @}
    
    /// @name accessors
    /// @{
    value_t x() const { return data[0]; }
    value_t y() const { return data[1]; }
    value_t z() const { return data[2]; }
    value_t r() const { return data[0]; }
    value_t g() const { return data[1]; }
    value_t b() const { return data[2]; }
    
    value_t& operator[](int comp) { return data[comp]; }
    const value_t& operator[](int comp) const { return data[comp]; }
    /// @}
    
    /// @name vector arithmetic
    /// @{
    Vec3& operator*=(const value_t mult)
    {
        data[0]*=mult; data[1]*=mult; data[2]*=mult;
        return *this;
    }

    Vec3& operator/=(const value_t denom)
    {
        assert( 0.f != denom );
        data[0]/=denom; data[1]/=denom; data[2]/=denom;
        return *this;
    }
    /// @}
    
    /// @name vector algebra
    /// @{
    value_t dot( const Vec3& b ) const { return x()*b.x() + y()*b.y() + z()*b.z(); }
    value_t length() const { return std::sqrt( dot(*this) ); }
    value_t length2() const { return dot(*this); }
    void
    normalize()
    {
        const value_t len2 = length2();
        if (0.f==len2)
            return;
        
        *this /= std::sqrt(len2);
    }
    
    static Vec3
    unit_vector( const Vec3& src )
    {
        Vec3 unit = src;
        unit.normalize();
        return unit;
    }
    /// @}
};

Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3(a.data[0] + b.data[0],
                                                           a.data[1] + b.data[1],
                                                           a.data[2] + b.data[2]); }

Vec3 operator*(const Vec3& a, const Vec3& b) { return Vec3(a.data[0] * b.data[0],
                                                           a.data[1] * b.data[1],
                                                           a.data[2] * b.data[2]); }

Vec3 operator*(const Vec3& a, const float t) { return Vec3(a.data[0]*t,
                                                           a.data[1]*t,
                                                           a.data[2]*t); }

Vec3 operator/(const Vec3& a, const float t) { return Vec3(a.data[0]/t,
                                                           a.data[1]/t,
                                                           a.data[2]/t); }

class Ray
{
public:
    Ray() {}
    Ray( const Vec3& _orig, const Vec3& _dir ) :
        m_origin( _orig ),
        m_direction( _dir )
    {}
    
    void
    setRay( const Vec3& _orig, const Vec3& _dir )
    {
        m_origin = _orig;
        m_direction = _dir;
    }
    
    const Vec3&
    direction() const { return m_direction; }
    
    const Vec3&
    origin() const { return m_origin; }

    Vec3
    point( const float t ) const
    {
        return m_direction * t + m_origin;
    }
    
private:
    Vec3 m_origin;
    Vec3 m_direction;
};


/// image data structure
class Framebuffer
{
public:
    Framebuffer(int w, int h) :
        m_width(w),
        m_height(h),
        m_buffer(w*h)
    {}
    
    const Vec3&
    operator()( int x, int y ) const
    {
        return m_buffer[m_width*y+x];
    }
    
    Vec3&
    operator()( int x, int y )
    {
        return m_buffer[m_width*y+x];
    }
    
    /// write a PPM file
    void
    write( const string& fname )
    {
        ofstream imgf( fname );
        imgf << "P3" << endl
             << m_width << " " << m_height << "\n"
             << "255\n" /* max color */;

        for (int r=m_height-1; r>=0; --r) {
            for (int c=0; c<m_width; ++c) {
                const Vec3& color = operator()(c,r);
                imgf << int(255.0 * color[0]) << " "
                     << int(255.0 * color[1]) << " "
                     << int(255.0 * color[2]) << "\n";
            }
        }
    }
    
    void
    test()
    {
        Framebuffer fb(m_width,m_height);
        Vec3 rgb;
        for (int r=0; r<m_height; ++r) {
            for (int c=0; c<m_width; ++c) {
                rgb[0] = float(r)/float(m_height);
                rgb[1] = float(c)/float(m_width);
                rgb[2] = 0.f;
                fb(c,r) = rgb;
            }
        }
        write( "test.ppm" );
    }
    
private:
    const int m_width;
    const int m_height;
    vector<Vec3> m_buffer;
};


/// Base class for a geometric entity in R3
class Geometry
{
public:
    virtual bool
    intersects( const Ray& ray ) const = 0;
};


/// Sphere
class Sphere : public Geometry
{
public:
    bool intersects( const Ray& ray ) const
    {
        
        return false;
    }
    
private:
};


/// Some physical material
class Material
{
public:
    
};

/// Something that can be drawn: geometry + material
class Drawable
{
public:
    
};

/// Data structure for the scene
class World
{
public:
    bool intersects( const Ray& ray )
    {
        /* for each object
            - test intersection
            - save iterator to closest object
           return closest
         */
        bool hit = false;
        
        return hit;
    }
private:
    std::list<Drawable> m_objects;
};


/// Camera class to project between worldspace, and filmback UV coordinates
class Camera
{
public:
    Vec3 corner;  ///< lower left corner
    Vec3 offset;  ///< offset to the right top corner

    /// Convert UV coords into world positions.  UVs must be in the [0,1] range.
    Vec3
    project( float u, float v ) const
    {
        Vec3 pixpos(u,v,0);
        return corner + offset * pixpos;
    }
};

/// sample the world with @a ray, and produce a color.
Vec3
sample( const World& scene, const Ray& ray )
{
    Vec3 color;
    Vec3 unit_v = Vec3::unit_vector(ray.direction());
    const float a = 0.5 * (unit_v.y()+1);
    color = (1.0-a) * Vec3(1,1,1) + a * Vec3(0.5,0.7,1.0);
    
    return color;
}

int main(int argc, const char * argv[]) {
    /*
     dev plan:
     - write test image (DONE)
     - render quad
     */
    const int W=600;
    const int H=300;
    Framebuffer frame(W,H);
    
    World scene;
    
    /* for each pixel:
        - create ray
        - sample ray
     */
    const Vec3 origin(0,0,0);
    Camera camera = { Vec3(-2,-1,-1), Vec3(4,2,0) };
    
    for (int j=0; j<H; ++j) {
        for (int i=0; i<W; ++i) { // for each pixel
            const float U = float(i)/float(W);
            const float V = float(j)/float(H);
            Ray ray( origin, camera.project(U, V) );
            const Vec3 color = sample( scene, ray );
            frame(i,j) = color;
        }
    }

    frame.write("test.ppm");
    
    return 0;
}
