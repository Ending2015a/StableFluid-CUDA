#include <iostream>
#include <vector>
#include <fstream>

#include <cassert>
#include <cmath>

#include "../lodepng/lodepng.h"

#define max(x, y) ((x)>(y)?(x):(y))

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


void decode(const char* filename, unsigned char **image, unsigned int *width, unsigned int *height)
{
    unsigned int error = lodepng_decode32_file(image, width, height, filename);
    if(error)
        printf("error %u: %s\n", error, lodepng_error_text(error));
}

// ./exe input dx dy rate output
int main(int argc, char **argv)
{
    assert(argc == 6);
    std::ofstream fout(argv[5], std::ios::binary);

    unsigned char* raw_image;
    unsigned int nx;
    unsigned int ny;

    decode(argv[1], &raw_image, &nx, &ny);

    double dx = atof(argv[2]);
    double dy = atof(argv[3]);
    double rate = atof(argv[4]);
    
    std::vector<vec2> point_list;

    unsigned char **image;
    image = new unsigned char*[ny]{};

    for(int i=0;i<ny;++i)
    {
        image[i] = raw_image + i * nx * 4;
    }

    for(int i=0;i<ny;++i)
    {
        for(int j=0;j<nx;++j)
        {
            int mx = max((int)(255-image[i][j*4]), (int)(255-image[i][j*4+1]));
            mx = max(mx, (int)(255-image[i][j*4+2]));
            
            if(mx > 20)  // sample !
            {
                int j1 = j;
                int i1 = ny - i;

                for(double rx = 0.; rx < rate; rx += 1.)
                {
                    for(double ry = 0.; ry < rate; ry += 1.)
                    {
                        double ox, oy;
                        hash22(rx/rate+(double)j1, ry/rate+(double)i1, ox, oy);
                        double px = j1 + (rx+ox)/rate;
                        double py = i1 + (ry+oy)/rate;

                        point_list.push_back(vec2({px*dx, py*dy}));
                    }
                }
            }
        }
    }
    
    int num = point_list.size();
    std::cout << "done! [" << num << "] samples" << std::endl;
    
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

    free(raw_image);
    delete[] image;

    return 0;
}
