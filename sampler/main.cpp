#include <iostream>
#include <vector>
#include <fstream>

#include <cassert>
#include <cmath>

#include "../lodepng/lodepng.h"

struct vec2
{
    double x, y;
};

double fract(double x)
{
    long long n = (long long)x;
    return x - (double)n;
}

double floor(double x)
{
    return (long long )x;
}

void hash22(double x, double y, double &ox, double &oy)
{
    double px = x*17.65465 + y*168.91949+54.156548;
    double py = x*184.64651 + y*517.54651+px;
    ox = fabs(fract(sin(px+151.1564)*65.1598465))*0.8+0.1;
    oy = fabs(fract(sin(py+248.4815)*65.7298465))*0.8+0.1;
}


void decode(const char*filename, unsigned char *image, unsigned int width, unsigned int height)
{
    unsigned int error = lodepng_decode32_file(&image, &width, &height, filename);
    if(error)
        printf("error %u: %s\n", error, lodepng_error_text(error));
}

// ./exe nx ny dx dy rate output j1 i1 j2 i2 ...
int main(int argc, char **argv)
{
    assert(argc > 11);
    std::ofstream fout(argv[6], std::ios::binary);

    int nx = atoi(argv[1]);
    int ny = atoi(argv[2]);
    double dx = atof(argv[3]);
    double dy = atof(argv[4]);
    double rate = atof(argv[5]);
    
    std::vector<vec2> point_list;
    
    for(int i=7;i < argc; i+=4)
    {
        int j1 = atoi(argv[i]);
        int i1 = atoi(argv[i+1]);
        int j2 = atoi(argv[i+2]);
        int i2 = atoi(argv[i+3]);

        double rate_x = (double)abs(j1-j2)*rate;
        double rate_y = (double)abs(i1-i2)*rate;

        for(double rx = 0.; rx < rate_x;rx += 1.)
        {
            for(double ry = 0.; ry < rate_y; ry += 1.)
            {
                double ox, oy;
                hash22(rx, ry, ox, oy);
                double px = j1 + (rx + ox)/rate_x * (j2-j1);
                double py = i1 + (ry + oy)/rate_y * (i2-i1);

                point_list.push_back(vec2({px*dx, py*dy}));
            }
        }
    }

    int num = point_list.size();
    std::cout << "done! " << num << "samples" << std::endl;
    
    fout.write((char*)&nx, sizeof(int));
    fout.write((char*)&ny, sizeof(int));
    fout.write((char*)&dx, sizeof(double));
    fout.write((char*)&dy, sizeof(double));
    fout.write((char*)&num, sizeof(int));

    for(auto &p:point_list)
    {
        fout.write((char*)&p.x, sizeof(double));
        fout.write((char*)&p.y, sizeof(double));
    }

    fout.close();

    return 0;
}
