#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <random>
#include <algorithm>
#include <string>

#include <cmath>
#include <cassert>
#include <cstring>
#include <ctime>
#include <cstdlib>

#include <lodepng.h>
#include <cmdline.h>

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



// ===== Themes =====
class Theme
{
public:
    bool draw_grid;

    virtual void render(Canvas &canvas, const char *filename) = 0;
    virtual void initialize() = 0;
    virtual void plot_arrow(double u, double v, 
                            double off_x, double off_y, 
                            double len_x, double len_y, 
                            double magni,
                            double angle,
                            Canvas &canvas) = 0;
    virtual void plot_boundary(double x1, double y1, 
                            double x2, double y2,
                            Canvas &canvas) = 0;
    virtual void plot_water(double x1, double y1,
                            double x2, double y2,
                            Canvas &canvas) = 0;
    virtual void plot_air(double x1, double y1,
                            double x2, double y2,
                            Canvas &canvas) = 0;
    virtual void plot_marker(double x1, double y1,
                            Canvas &canvas) = 0;
    virtual void plot_background(Canvas &canvas) = 0;
    virtual Color get_color_palette(double t) = 0;
};

class Normal : public Theme
{
public:
    int marker_size;
    Color black;
    Color white;
    Color gray;
    Color water;
    Color marker;
    Color background;

    virtual void render(Canvas &canvas, const char *filename) override
    {
        Image& image = canvas.image;
        image.save(filename);
    }

    virtual void initialize() override
    {
        black = Color({0., 0., 0., 1.});
        white = Color({1., 1., 1., 1.});
        marker = Color({0., 0.18, 0.83, 1.});
        water = Color({0.40, 0.64, 0.98, 1.});
        gray = Color({0.5, 0.5, 0.5, 1.});

        background = white;

        marker_size = 2;
    }

    virtual void plot_arrow(double u, double v,
                            double off_x, double off_y,
                            double len_x, double len_y,
                            double magni,
                            double angle, 
                            Canvas &canvas)
    {
        double unorm = u / magni * len_x;
        double vnorm = v / magni * len_y;

        double x1 = off_x - unorm;
        double y1 = off_y - vnorm;
        double x2 = off_x + unorm;
        double y2 = off_y + vnorm;
        
        Color cc = get_color_palette(magni);

        canvas.drawLine(x1, y1, x2, y2, cc);
        
        angle = angle/180. * 3.14159265358979323846;

        {
            double ax = ((x1-x2) * cos(angle) - (y1-y2) * sin(angle))/2.;
            double ay = ((x1-x2) * sin(angle) + (y1-y2) * cos(angle))/2.;
            canvas.drawLine(ax+x2, ay+y2, x2, y2, cc);
        }

        {
            double ax = ((x1-x2) * cos(-angle) - (y1-y2) * sin(-angle))/2.;
            double ay = ((x1-x2) * sin(-angle) + (y1-y2) * cos(-angle))/2.;
            canvas.drawLine(ax+x2, ay+y2, x2, y2, cc);
        }    
    }

    virtual void plot_boundary(double x1, double y1, 
                            double x2, double y2,
                            Canvas &canvas) override
    {
        if(draw_grid)
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, white, gray);
        else
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, gray, gray);
    }

    virtual void plot_water(double x1, double y1,
                            double x2, double y2,
                            Canvas &canvas) override
    {
        if(draw_grid)
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, white, water);
        else    
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, water, water);
    }

    virtual void plot_air(double x1, double y1,
                            double x2, double y2,
                            Canvas &canvas) override
    {
        canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, white, white);
    }

    virtual void plot_background(Canvas &canvas) override
    {
        canvas.fill(background);
    }

    virtual void plot_marker(double x1, double y1, Canvas &canvas) override
    {
        canvas.drawPoint(x1, y1, marker_size, marker);
    }


    virtual Color get_color_palette(double t)
    {
        return black;
    }
};



class Nightmare : public Theme
{
public:

    int marker_size;

    Color background;
    Color black;
    Color gray;
    Color gray2;
    Color water;
    Color marker;

    virtual void render(Canvas &canvas, const char *filename) override
    {
        Image &image = canvas.image;
        image.save(filename);
    }

    virtual void initialize() override
    {
        background = Color({0., 0., 0., 1.});
        black = background;
        gray = Color({0.2, 0.2, 0.2, 1.0});
        gray2 = Color({0.15, 0.15, 0.15, 1.0});
        water = Color({0.64, 0.30, 0.58, 1.});
        marker = Color({0.94, 0.67, 0.87, 1.});

        marker_size = 2;
    }

    virtual void plot_arrow(double u, double v,  //normalized arrow vector
                            double off_x, double off_y,  //center position
                            double len_x, double len_y, //arrow length
                            double magni,  //magitude //0~1
                            double angle, //arrow angle
                            Canvas &canvas
                            ) override
    {
        double sq = sqrt(u*u+v*v);
        double unorm = u/sq * len_x;
        double vnorm = v/sq * len_y;

        double x1 = off_x - unorm;
        double y1 = off_y - vnorm;
        double x2 = off_x + unorm;
        double y2 = off_y + vnorm;
        
        Color cc = get_color_palette(sq/magni);

        canvas.drawLine(x1, y1, x2, y2, cc);
        
        angle = angle/180. * 3.14159265358979323846;

        {
            double ax = ((x1-x2) * cos(angle) - (y1-y2) * sin(angle))/2.;
            double ay = ((x1-x2) * sin(angle) + (y1-y2) * cos(angle))/2.;
            canvas.drawLine(ax+x2, ay+y2, x2, y2, cc);
        }

        {
            double ax = ((x1-x2) * cos(-angle) - (y1-y2) * sin(-angle))/2.;
            double ay = ((x1-x2) * sin(-angle) + (y1-y2) * cos(-angle))/2.;
            canvas.drawLine(ax+x2, ay+y2, x2, y2, cc);
        }
    }

    virtual void plot_boundary(double x1, double y1, 
                            double x2, double y2,
                            Canvas &canvas) override
    {
        canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, gray2, gray2);
    }

    virtual void plot_water(double x1, double y1,
                            double x2, double y2,
                            Canvas &canvas) override
    {
        canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, water, water);
    }

    virtual void plot_air(double x1, double y1,
                            double x2, double y2,
                            Canvas &canvas) override
    {
        if(draw_grid)
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, gray2, background);
        else
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, background, background);
    }

    virtual void plot_background(Canvas &canvas) override
    {
        canvas.fill(background);
    }

    virtual void plot_marker(double x1, double y1, Canvas &canvas) override
    {
        canvas.drawPoint(x1, y1, marker_size, marker);
    }

    virtual Color get_color_palette(double t) override
    {
        t = (1.-t)*(-0.8);
        Color g({0., 0., 0., 1.});
        Color a({0.5, 0.5, 0.5, 0.});
        Color b({0.5, 0.5, 0.5, 0.});
        Color c({1., 1., 1., 0.});
        Color d({0.0, 0.33, 0.67, 0.});

        g.r = clamp(a.r + b.r * cos(2.*3.14159625358979323846*(c.r*t + d.r)), 0., 1.);
        g.g = clamp(a.g + b.g * cos(2.*3.14159625358979323846*(c.g*t + d.g)), 0., 1.);
        g.b = clamp(a.b + b.b * cos(2.*3.14159625358979323846*(c.b*t + d.b)), 0., 1.);

        return g;
    }
};


class Paper : public Theme
{
public:
    int marker_size;
    Color black;
    Color white;
    Color gray;
    Color water;
    Color marker;
    Color background;

    double *noise=0;

    ~Paper()
    {
        delete[] noise;
    }

    double sample(double i, double j)
    {
        i -= floor(i);
        j -= floor(j);
        i *= 255.;
        j *= 255.;
        
        int x1 = (int)floor(j);
        int y1 = (int)floor(i);
        int x2 = (int)clamp(x1+1, 0, 255);
        int y2 = (int)clamp(y1+1, 0, 255);

        double a = j - x1;
        double b = i - y1;

        double yy1 = noise[y1 * 256 + x1] * (1-a) + noise[y1 * 256 + x2] * a;
        double yy2 = noise[y2 * 256 + x1] * (1-a) + noise[y2 * 256 + x2] * a;

        return yy1 * (1-b) + yy2 * b;
    }

    virtual void render(Canvas &canvas, const char *filename) override
    {
        Image& image = canvas.image;
#pragma omp parallel for num_threads(8) schedule(dynamic, 4) shared(image)
        for(unsigned i=0;i<image.height; ++i)
        {
            double dy = i/(double)image.height;
            for(unsigned j=0;j<image.width;++j)
            {
                double dx = j /(double)image.width;

                Color col = image(j, i);
                double p = pow(16.0*dx*dy*(1.0-dx)*(1.0-dy), 0.2);
                p *= (0.9 + 0.1 * sample((double)i/256., (double)j/256.));

                col.r = pow(col.r * p, 0.4545);
                col.g = pow(col.g * p, 0.4545);
                col.b = pow(col.b * p, 0.4545);
                image(j, i) = col;
            }
        }

        image.save(filename);
    }

    virtual void initialize() override
    {
        black = Color({0.1, 0.1, 0.1, 1.});
        white = Color({1., 1., 1., 1.});
        marker = Color({0., 0.18, 0.83, 1.});
        water = Color({0.29, 0.29, 1., 1.});
        gray = Color({0.5, 0.5, 0.5, 1.});

        background = white;

        noise = new double[256*256]{};
        
        for(int i=0;i<256*256;++i)
            noise[i] = rand()/(double)RAND_MAX;

        marker_size = 2;
    }

    virtual void plot_arrow(double u, double v,
                            double off_x, double off_y,
                            double len_x, double len_y,
                            double magni,
                            double angle, 
                            Canvas &canvas)
    {
        double unorm = u / magni * len_x;
        double vnorm = v / magni * len_y;

        double x1 = off_x - unorm;
        double y1 = off_y - vnorm;
        double x2 = off_x + unorm;
        double y2 = off_y + vnorm;
        
        Color cc = get_color_palette(magni);

        canvas.drawLine(x1, y1, x2, y2, cc);
        
        angle = angle/180. * 3.14159265358979323846;

        {
            double ax = ((x1-x2) * cos(angle) - (y1-y2) * sin(angle))/2.;
            double ay = ((x1-x2) * sin(angle) + (y1-y2) * cos(angle))/2.;
            canvas.drawLine(ax+x2, ay+y2, x2, y2, cc);
        }

        {
            double ax = ((x1-x2) * cos(-angle) - (y1-y2) * sin(-angle))/2.;
            double ay = ((x1-x2) * sin(-angle) + (y1-y2) * cos(-angle))/2.;
            canvas.drawLine(ax+x2, ay+y2, x2, y2, cc);
        }    
    }

    virtual void plot_boundary(double x1, double y1, 
                            double x2, double y2,
                            Canvas &canvas) override
    {
        if(draw_grid)
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, white, gray);
        else
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, gray, gray);
    }

    virtual void plot_water(double x1, double y1,
                            double x2, double y2,
                            Canvas &canvas) override
    {
        if(draw_grid)
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, white, water);
        else
            canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, water, water);
    }

    virtual void plot_air(double x1, double y1,
                            double x2, double y2,
                            Canvas &canvas) override
    {
        canvas.fillRect((int)x1, (int)y1, (int)x2, (int)y2, white, white);
    }

    virtual void plot_background(Canvas &canvas) override
    {
        canvas.fill(background);
    }

    virtual void plot_marker(double x1, double y1, Canvas &canvas) override
    {
        canvas.drawPoint(x1, y1, marker_size, marker);
    }


    virtual Color get_color_palette(double t)
    {
        return black;
    }
};



// ===== Renderer =====
template<typename Themes>
class Renderer
{
public:
    Themes themes;
    
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

    std::vector<double> u;
    std::vector<double> v;

    bool draw_markers;
    bool draw_field;
    bool draw_grid;

    Renderer(unsigned w, unsigned h, 
            int nx, int ny, double dx, double dy,
            std::vector<vec2> &particle_list, 
            std::vector<double> &u, std::vector<double> &v) 
        : canvas(w, h), grid(nx, ny, dx, dy, particle_list),
          nx(nx), ny(ny), dx(dx), dy(dy), particle_size(2),
          margin_x(0.9), margin_y(0.9),
          u(u), v(v)
    { 
        grid.update(); 
        themes.initialize();
        themes.plot_background(canvas);
    }

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
            themes.plot_boundary(pivot_x+ddx*i, pivot_y, pivot_x+ddx*(i+1), pivot_y+ddy, canvas);
            themes.plot_boundary(pivot_x+ddx*i, pivot_y+ddy*(ny+1), pivot_x+ddx*(i+1), pivot_y+ddy*(ny+2), canvas);
        }

        for(int i=0;i<ny+2;++i)
        {
            themes.plot_boundary(pivot_x, pivot_y+ddy*i, pivot_x+ddx, pivot_y+ddy*(i+1), canvas);
            themes.plot_boundary(pivot_x+ddx*(nx+1), pivot_y+ddy*i, pivot_x+ddx*(nx+2), pivot_y+ddy*(i+1), canvas);
        }

#pragma omp parallel for num_threads(8) schedule(dynamic, 4)
        for(int i=0;i<ny;++i)
        {
            for(int j=0;j<nx;++j)
            {

                int off_x = pivot_x + ddx * (j+1);
                int off_y = pivot_y + ddy * (i+1);

                if(grid(j, i) == 1)  //fluid
                    themes.plot_water(off_x, off_y, off_x+ddx, off_y+ddy, canvas);
                else //air
                    themes.plot_air(off_x, off_y, off_x+ddx, off_y+ddy, canvas);

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
            themes.plot_marker(ix, iy, canvas);
        }
    }

    void drawField()
    {
        double ddx = (double)len_x/(double)(nx+2);
        double ddy = (double)len_y/(double)(ny+2);

        if(u.size() == 0 || v.size() == 0)
            return;

        double max_magni = 0;
        for(int i=0;i<ny;++i)
        {
            for(int j=0;j<nx;++j)
            {
                double fu = (u[i*(nx+1) + j+1] + u[i*(nx+1) + j])/2.;
                double fv = (v[(i+1)*nx+j] + v[i*nx+j])/2.;
                double m = sqrt(fu*fu+fv*fv);

                max_magni = max_magni < m ? m:max_magni;
            }
        }

#pragma omp parallel for num_threads(8) schedule(dynamic, 4) shared(max_magni)
        for(int i=0;i<ny;++i)
        {
            for(int j=0;j<nx;++j)
            {
                double ufield = (u[i*(nx+1) + j+1] + u[i*(nx+1) + j])/2.;
                double vfield = (v[(i+1)*nx+j] + v[i*nx+j])/2.;
                //double m = sqrt(fu*fu+fv*fv);

                int off_x = pivot_x + ddx * (j+1);
                int off_y = pivot_y + ddy * (i+1);


                themes.plot_arrow(ufield, vfield,
                                off_x + ddx/2., off_y + ddy/2.,
                                ddx/5., ddy/5.,
                                max_magni,
                                30, 
                                canvas);    

            }
        }
    }

    void render()
    {
        themes.draw_grid = draw_grid;

        analyzeGrid();
        drawGrid();
        if(draw_markers)
            drawParticles();
        if(draw_field)
            drawField();
    }

    void save(std::string filename){ themes.render(canvas, filename.c_str()); }
};



// ./exe input output width height [field] [themes]
int main(int argc, char **argv)
{

    cmdline::parser options;

    options.add<std::string>("input", 'i', "Input snapshot file", true);
    options.add<std::string>("output", 'o', "Output PNG file", true);
    options.add<unsigned>("width", '\0', "Output width", false, 1920);
    options.add<unsigned>("height", '\0', "Output height", false, 1080);
    options.add("velocity", 'v', "Draw velocity field");
    options.add("marker", 'm', "Draw markers");
    options.add("grid", 'g', "Draw grid");
    options.add<std::string>("theme", 't', "Theme (None/Paper/Nightmare)", false, "None",
                        cmdline::oneof<std::string>("None", "Paper", "Nightmare"));

    //cxxopts::Options options(argv[0], "Homework 6 - Fluid Renderer");

    //options
    //    .add_options()
    //    ("i,input", "Input file, a snapshot file", cxxopts::value<std::string>())
    //    ("o,output", "Output file, in PNG format", cxxopts::value<std::string>())
    //    ("width", "Output width", cxxopts::value<unsigned>()->default_value("1920"))
    //    ("height", "Output height", cxxopts::value<unsigned>()->default_value("1080"))
    //    ("v,velocity", "Draw velocity field")
    //    ("m,marker", "Draw markers")
    //    ("g,grid", "Draw grid")
    //    ("t,theme", "Theme (None/Paper/Nightmare)", cxxopts::value<std::string>()->default_value("None"))
    //    ("h,help", "Print help");
    //assert(argc == 5 || argc == 6 || argc == 7);

    options.parse_check(argc, argv);
    //auto result = options.parse(argc, argv);

    //if(result.count("help"))
    //{
    //    std::cout << options.help({"", "Group"}) << std::endl;
    //    exit(0);
   // }

    //if(!result.count("i") || !result.count("o"))
   // {
   //     std::cout << "[ERROR] Please specifiy the input and output file" << std::endl;
   //     std::cout << std::endl;
   //     std::cout << options.help() << std::endl;
   //     exit(0);
   // }

    //std::cout << "Input: " << result["i"].as<std::string>() << std::endl;
    //std::cout << "Output: " << result["o"].as<std::string>() << std::endl;
    //std::cout << "Width: " << result["width"].as<unsigned>() << std::endl;
    //std::cout << "Height: " << result["height"].as<unsigned>() << std::endl;
    //std::cout << "Draw velocity field: " << std::boolalpha << (result.count("f") > 0) << std::endl;
    //std::cout << "Draw markers: " << std::boolalpha << (result.count("m") > 0) << std::endl;
    //std::cout << "Draw grid flag: " << std::boolalpha << (result.count("g") > 0) << std::endl;
    //std::cout << "Theme: " << result["t"].as<std::string>() << std::endl;


    std::string input_file = options.get<std::string>("input");// result["i"].as<std::string>();
    std::string output_file = options.get<std::string>("output");//result["o"].as<std::string>();
    unsigned width = options.get<unsigned>("width");//result["width"].as<unsigned>();
    unsigned height = options.get<unsigned>("height");//result["height"].as<unsigned>();
    bool draw_markers = options.exist("marker"); //result["m"].as<bool>();
    bool draw_field = options.exist("velocity"); //result["v"].as<bool>();
    bool draw_grid = options.exist("grid"); //result["g"].as<bool>();
    std::string theme = options.get<std::string>("theme"); //result["t"].as<std::string>();


    std::ifstream fin(input_file, std::ios::binary);
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

    std::vector<double> field_u;
    std::vector<double> field_v;

    if(draw_field)
    {
        field_u.reserve( (nx+1)*ny );
        field_v.reserve( nx*(ny+1) );
        double tmp = 0;
        for(int i=0;i<(nx+1)*ny;++i)
        {
            fin.read((char*)&tmp, sizeof(double));
            if(fin.eof())
                break;
            field_u.push_back(tmp);
        }
        for(int i=0;i<nx*(ny+1);++i)
        {
            fin.read((char*)&tmp, sizeof(double));
            if(fin.eof())
                break;
            field_v.push_back(tmp);
        }
    }

    fin.close();

    if(theme == "Nightmare")
    {
        Renderer<Nightmare> renderer(width, height, nx, ny, dx, dy, particles, field_u, field_v);
        renderer.draw_markers = draw_markers;
        renderer.draw_field = draw_field;
        renderer.draw_grid = draw_grid;
        renderer.render();
        renderer.save(output_file);
    }else if(theme == "Paper"){
        Renderer<Paper> renderer(width, height, nx, ny, dx, dy, particles, field_u, field_v);
        renderer.draw_markers = draw_markers;
        renderer.draw_field = draw_field;
        renderer.draw_grid = draw_grid;
        renderer.render();
        renderer.save(output_file);
    }else
    {
        Renderer<Normal> renderer(width, height, nx, ny, dx, dy, particles, field_u, field_v);
        renderer.draw_markers = draw_markers;
        renderer.draw_field = draw_field;
        renderer.draw_grid = draw_grid;
        renderer.render();
        renderer.save(output_file);
    }

    return 0;
}
