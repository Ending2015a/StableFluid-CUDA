#include <iostream>
#include <vector>
#include <fstream>

#include <cmath>
#include <cassert>
#include <cstring>


#include "lodepng/lodepng.h"


template<typename T>
T clamp(T x, T mn, T mx)
{
    return x<mn?mn:x>mx?mx:x;
}

// ===== vec2 =====
struct vec2
{
    double x, y;
};

// ===== Grid =====
class Grid
{
public:
    std::vector<vec2> particle_list;
    int *grid_state;

    int nx;
    int ny;
    double dx;
    double dy;

    Grid(int nx, int ny, double dx, double dy, 
            std::vector<vec2> &particles)
        : nx(nx), ny(ny), dx(dx), dy(dy), particle_list(particles)
    { grid_state = new int[nx*ny]{}; }

    ~Grid(){ delete[] grid_state; }

    void update()
    {
        memset((char*)grid_state, 0, nx*ny*sizeof(int));
        for(auto& p:particle_list)
        {
            int ix = (int)(p.x/dx);
            int iy = (int)(p.y/dy);
            if(ix < 0 || iy < 0 || ix > nx-1 || iy > ny-1)
                continue;
            grid_state[iy*nx+ix] = 1;
        }
    }

    int &operator()(int x, int y)
    {
        return grid_state[y*nx+x];
    }
};

// ===== Color =====
struct Color
{
    double r, g, b, a;

    void fill(double n){ r = g = b = a = n; }
    void fill(Color &n){ r = n.r, g = n.g, b = n.b, a = n.a; }
};

// ===== Image =====
class Image
{
public:
    unsigned width;
    unsigned height;
    Color *raw_image;
    Color **image;

    Image(unsigned w, unsigned h):width(w), height(h)
    { 
        raw_image = new Color[w*h]{};
        image = new Color*[h]{};
        for(unsigned i=0;i<h;++i)
        {
            image[i] = raw_image + i*w;
        }
    }
    ~Image(){ delete[] raw_image; delete[] image;}

    Color& operator()(unsigned w, unsigned h){ return image[h][w]; }

    void save(const char * filename)
    {
        unsigned char *raw = new unsigned char[width*height*4]{};

        for(int i=0;i<height;++i)
        {
            for(int j=0;j<width;++j)
            {
                unsigned p = (height-1-i)*width*4+j*4;
                raw[p] = (unsigned char)(image[i][j].r*255.);
                raw[p+1] = (unsigned char)(image[i][j].g*255.);
                raw[p+2] = (unsigned char)(image[i][j].b*255.);
                raw[p+3] = 255;
            }
        }

        unsigned err = lodepng_encode32_file(filename, (unsigned char*)raw, width, height);
        if(err)
            std::cout << "[lodepng] error " << err << ": " << lodepng_error_text(err);

        delete[] raw;
    }
};

// ===== Canvas =====
class Canvas
{
public:

    unsigned width;
    unsigned height;
    Image image;

    Canvas(unsigned w, unsigned h) : width(w), height(h), image(w, h){}

    void drawLine(int x1, int y1, int x2, int y2, Color &c)
    {
        double dx = x2-x1;
        double dy = y2-y1;

        double dist = sqrt(dx*dx+dy*dy);

        for(int i=0;i<dist;++i)
        {
            int xi = x1 + dx/dist * i;
            int yi = y1 + dy/dist * i;

            image(clamp(xi, 0, (int)width-1), clamp(yi, 0, (int)height-1)).fill(c);
        }
    }

    void drawPoint(int x, int y, int size, Color &c)
    {
        int x1 = x-size/2;
        int x2 = x+size/2;
        int y1 = y-size/2;
        int y2 = y+size/2;
        fillRect(x1, y1, x2, y2, c, c);
    }

    void fillRect(int x1, int y1, int x2, int y2, Color &margin, Color &c)
    {
        x1 = clamp(x1, 0, (int)width-1);
        x2 = clamp(x2, 0, (int)width-1);
        y1 = clamp(y1, 0, (int)height-1);
        y2 = clamp(y2, 0, (int)height-1);
        for(int i=y1;i<=y2;i++)
        {
            for(int j=x1;j<=x2;j++)
            {
                if(i == y1||i==y2||j==x1||j==x2)
                    image(j, i).fill(margin);
                else
                    image(j, i).fill(c);
            }
        }
    }

    void fill(Color &c)
    {
        for(unsigned y=0;y<height;++y)
        {
            for(unsigned x=0;x<width;++x)
            {
                image(x, y).fill(c);
            }
        }
    }
};


// ===== Renderer =====
class Renderer
{
public:
    Color white;
    Color gray;
    Color blue1;
    Color blue2;

    Canvas canvas;
    Grid grid;

    int particle_size;
    int pivot_x;
    int pivot_y;
    int len_x;
    int len_y;

    double margin_x;
    double margin_y;
    
    int nx;
    int ny;
    double dx;
    double dy;

    Renderer(unsigned w, unsigned h, 
            int nx, int ny, double dx, double dy,
            std::vector<vec2> &particle_list) 
        : canvas(w, h), grid(nx, ny, dx, dy, particle_list),
          nx(nx), ny(ny), dx(dx), dy(dy), particle_size(2),
          white({1., 1., 1., 1.}),
          gray({.5, .5, .5, 1.}),
          blue1({0., 0., 1., 1.}),
          blue2({0.5, 0.5, 1., 1.}),
          margin_x(0.9), margin_y(0.9)
    { grid.update(); canvas.fill(white); }

    void analyzeGrid()
    {
        double w = (double)canvas.width*0.9;
        double h = (double)canvas.height*0.9;
        double c_ratio = w/h;

        double g_ratio = ((nx+2)*dx)/((ny+2)*dy);
        if(c_ratio > g_ratio)
        {
            len_y = (int)(canvas.height*margin_y);
            len_x = (int)(len_y * g_ratio);
            pivot_y = (int)(canvas.height* (1.-margin_x)/2. );
            pivot_x = (canvas.width-len_x)/2;
        }
        else
        {
            len_x = (int)(canvas.width * margin_x);
            len_y = (int)(len_x/g_ratio);
            pivot_x = (int)(canvas.width * (1.-margin_y)/2.);
            pivot_y = (canvas.height-len_y)/2;
        }
        
    }

    void drawGrid()
    {
        // grid cell length
        double ddx = (double)len_x/(double)(nx+2);
        double ddy = (double)len_y/(double)(ny+2);

        // top bottom boundary
        for(int i=0;i<nx+2;++i)
        {
            canvas.fillRect(pivot_x+ddx*i, pivot_y, pivot_x+ddx*(i+1), pivot_y+ddy, white, gray);
            canvas.fillRect(pivot_x+ddx*i, pivot_y+ddy*(ny+1), pivot_x+ddx*(i+1), pivot_y+ddy*(ny+2), white, gray);
        }

        for(int i=0;i<ny+2;++i)
        {
            canvas.fillRect(pivot_x, pivot_y+ddy*i, pivot_x+ddx, pivot_y+ddy*(i+1), white, gray);
            canvas.fillRect(pivot_x+ddx*(nx+1), pivot_y+ddy*i, pivot_x+ddx*(nx+2), pivot_y+ddy*(i+1), white, gray);
        }

        Color fillc;

        for(int i=0;i<ny;++i)
        {
            for(int j=0;j<nx;++j)
            {
                if(grid(j, i) == 1)  //FLUID
                    fillc = blue2;
                else   //AIR
                    fillc = white;

                int off_x = pivot_x + ddx * (j+1);
                int off_y = pivot_y + ddy * (i+1);

                canvas.fillRect(off_x, off_y, off_x+ddx, off_y+ddy, white, fillc);
            }
        }

    }

    void drawParticles()
    {
        double ddx = (double)len_x/(double)(nx+2);
        double ddy = (double)len_y/(double)(ny+2);


        for(auto& p:grid.particle_list)
        {
            int ix = (p.x/dx+1)*ddx + pivot_x;
            int iy = (p.y/dy+1)*ddy + pivot_y;
            canvas.drawPoint(ix, iy, particle_size, blue1);
        }
    }

    void render()
    {
        analyzeGrid();
        drawGrid();
        drawParticles();
    }

    void save(const char *filename){ canvas.image.save(filename); }
};



// ./exe input output width height
int main(int argc, char **argv)
{
    assert(argc == 5);

    unsigned width = atoi(argv[3]);
    unsigned height = atoi(argv[4]);

    std::ifstream fin(argv[1], std::ios::binary);
    int nx;
    int ny;
    double dx;
    double dy;
    int num;

    fin.read((char*)&nx, sizeof(int));
    fin.read((char*)&ny, sizeof(int));
    fin.read((char*)&dx, sizeof(double));
    fin.read((char*)&dy, sizeof(double));
    fin.read((char*)&num, sizeof(int));
    

    std::vector<vec2> particles;
    double x, y;

    for(int i=0;i<num;++i)
    {
        fin.read((char*)&x, sizeof(double));
        fin.read((char*)&y, sizeof(double));
        particles.push_back({x, y});
    }

    fin.close();


    Renderer renderer(width, height, nx, ny, dx, dy, particles);

    renderer.render();
    renderer.save(argv[2]);

    return 0;
}
