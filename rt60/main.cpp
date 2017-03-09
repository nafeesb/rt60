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
#include <limits>
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
    Vec3& operator+=(const Vec3& rhs)
    {
        data[0] += rhs[0];
        data[1] += rhs[1];
        data[2] += rhs[2];

        return *this;
    }

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

    void
    reverse()
    {
        data[0] = -data[0];
        data[1] = -data[1];
        data[2] = -data[2];
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


/// @name Binary operations on vectors
/// @{

/// addition
Vec3 operator+(const Vec3& a, const Vec3& b) { return Vec3(a.data[0] + b.data[0],
                                                           a.data[1] + b.data[1],
                                                           a.data[2] + b.data[2]); }

/// subtraction
Vec3 operator-(const Vec3& a, const Vec3& b) { return Vec3(a.data[0] - b.data[0],
                                                           a.data[1] - b.data[1],
                                                           a.data[2] - b.data[2]); }


/// multiplication
Vec3 operator*(const Vec3& a, const Vec3& b) { return Vec3(a.data[0] * b.data[0],
                                                           a.data[1] * b.data[1],
                                                           a.data[2] * b.data[2]); }

/// scale
Vec3 operator*(const Vec3& a, const float t) { return Vec3(a.data[0]*t,
                                                           a.data[1]*t,
                                                           a.data[2]*t); }

/// scale by division
Vec3 operator/(const Vec3& a, const float t) { return a*(1.f/t); }

/// negation
Vec3 operator-(const Vec3& a) { return Vec3(-a.data[0],
                                            -a.data[1],
                                            -a.data[2]); }

/// Scalar product
Vec3::value_type dot( const Vec3& a, const Vec3& b ) { return a.dot( b ); }

/// @}


/// Represents a ray. Defined by an origin, and a normalized direction.
class Ray
{
public:
    Ray() {}
    Ray( const Vec3& _orig, const Vec3& _dir ) :
        m_origin( _orig ),
        m_direction( _dir )
    { m_direction.normalize(); }

    void
    setRay( const Vec3& _orig, const Vec3& _dir )
    {
        m_origin = _orig;
        m_direction = _dir;
        m_direction.normalize();
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


/// Image data structure
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


// forward declaration
class Drawable;
class Geometry;


/// Hit point
class Intersection
{
public:
    float t; ///< Ray parameter at intersection
    Vec3  p; ///< Intersection position
    Vec3  n; ///< Normal
    const Geometry* obj;
};


/// Some physical material
class Material
{
public:

};


/// Base class for a geometric entity in R3
class Geometry
{
public:
    /// dtor
    virtual ~Geometry() {};

    /// generic intersection method
    virtual bool
    intersects( const Ray& ray,
                const float t_min, const float t_max,
                Intersection& hit  ) const = 0;
};


/// Sphere
class Sphere : public Geometry
{
public:
    Sphere( const Vec3& _pos, const float _rad ) :
        m_pos(_pos),
        m_radius(_rad),
        m_radius_squared( _rad * _rad )
    {}

    virtual bool
    intersects( const Ray& ray, const float t_min, const float t_max, Intersection& hit ) const
    {
        // use quadratic equation to compute ray-sphere intersection points
        const Vec3 rc = ray.origin() - m_pos;
        const float a = dot( ray.direction(), ray.direction() );
        const float b = 2.f * dot( ray.direction(), rc );
        const float c = dot( rc, rc ) - m_radius_squared;

        const float discriminant = b*b - 4.f * a * c;
        if (discriminant>0) {
            // real solutions
            float t = (-b - sqrt(discriminant)) / (2.f*a);
            if( (t>t_min) && (t<t_max) ) {
                hit.t = t;
                hit.p = ray.point(t);
                hit.n = (hit.p-m_pos)*(1.f/m_radius);
                hit.obj = this;
                return true;
            }


            t = (-b + sqrt(discriminant)) / (2.f*a);
            if( (t>t_min) && (t<t_max) ) {
                hit.t = t;
                hit.p = ray.point(t);
                hit.n = (hit.p-m_pos)*(1.f/m_radius);
                hit.obj = this;
                return true;
            }

        }

        return false;
    }

    /// returns a random point inside a unit sphere
    static Vec3
    random_point_inside()
    {
        // use rejection sampling

        Vec3 pos;
        do {
            pos = Vec3(drand48(), drand48(), drand48()) * 2.f - Vec3(1);
        } while( dot(pos,pos) >= 1.f );

        return pos;
    }

private:
    Vec3             m_pos;
    Vec3::value_type m_radius;
    Vec3::value_type m_radius_squared;
};


/// Data structure for the scene
class World
{
public:

    ~World()
    {
        // delete list elements in a single pass
        while (not m_scene.empty()) {
            delete m_scene.front();
            m_scene.pop_front();
        }
    }

    bool intersects( const Ray& ray, const float t_min, const float t_max, Intersection& nearest ) const
    {
        /* for each object
            - test intersection
            - save iterator to closest object
           return closest
         */
        bool hit_something = false;
        Intersection testpt;
        nearest.t = std::numeric_limits<float>::max();

        for (const auto& obj: m_scene) {
            bool hit = obj->intersects(ray, t_min, nearest.t, testpt);
            if (true==hit) {
                nearest = testpt;
                hit_something = true;
            }
        }

        return hit_something;
    }

    void
    addDrawable( Geometry* obj )
    {
        m_scene.push_back(obj);
    }

private:
    typedef std::list<Geometry*> scene_t;
    scene_t m_scene;
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


/// sample the @a world with @a ray, and produce a color.
Vec3
sample( const World& scene, const Ray& ray )
{
    Vec3 color;
    Intersection nearest;

    bool hit = scene.intersects(ray, 0.f, numeric_limits<float>::max(), nearest);

    if (hit) {
        // lambertian
        Vec3 bounce = nearest.n + Sphere::random_point_inside();
        bounce.normalize();
        color = Vec3(0.3,0.3,0.6) * sample( scene, Ray(nearest.p, bounce) );
    }
    else {
        Vec3 unit_v = Vec3::unit_vector(ray.direction());
        const float a = 0.5 * (unit_v.y()+1);
        color = (1.0-a) * Vec3(1,1,1) + a * Vec3(0.5,0.7,1.0);
    }
    return color;
}


int
main(int argc, const char * argv[])
{
    srand48(0xDEADBEEF);

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
    const int num_samples = 16;
    const Vec3 origin(0,0,0);
    Camera camera = { Vec3(-2,-1,-1), Vec3(4,2,0) };

    // create a scene
    scene.addDrawable( new Sphere( Vec3(0,0,-1), 0.5 ));
    scene.addDrawable( new Sphere( Vec3(0,-50.5,-1), 50 ));

    for (int j=0; j<H; ++j) {
        for (int i=0; i<W; ++i) { // for each pixel

            Vec3 color;
            for (int s=0; s<num_samples; ++s) {
                const float U = float(i+drand48())/float(W);
                const float V = float(j+drand48())/float(H);
                Ray ray( origin, camera.project(U, V) );
                color += sample( scene, ray );
            }

            color /= num_samples;

            frame(i,j) = color;
        }
    }

    frame.write("scene.ppm");

    return 0;
}
