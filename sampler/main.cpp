#include <iostream>
#include <vector>
#include <fstream>
#include <string>

#include <cassert>
#include <cmath>

#include <lodepng.h>
#include <cmdline.h>

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


void decode(std::string filename, unsigned char **image, unsigned int *width, unsigned int *height)
{
    unsigned int error = lodepng_decode32_file(image, width, height, filename.c_str());
    if(error)
        printf("error %u: %s\n", error, lodepng_error_text(error));
}

// ./exe input dx dy rate output
int main(int argc, char **argv)
{
    cmdline::parser options;

    options.add<std::string>("input", 'i', "Input file, a binary image in PNG format", true);
    options.add<std::string>("output", 'o', "Output snapshot file", true);
    options.add<double>("dx", 'x', "Delta x", false, 0.1);
    options.add<double>("dy", 'y', "Delta y", false, 0.1);
    options.add<unsigned>("sr", 'r', "Sampling rate", false, 3, cmdline::range(1, 5));

    //cxxopts::Options options(argv[0], "Homework 6 - Fluid Sampler");

    //options
    //    .add_options()
    //    ("i,input", "Input file, a binary image in PNG format", cxxopts::value<std::string>())
    //    ("o,output", "Output file, snapshop", cxxopts::value<std::string>())
    //    ("x,dx", "Delta x", cxxopts::value<double>()->default_value("0.1"))
    //    ("y,dy", "Delta y", cxxopts::value<double>()->default_value("0.1"))
    //    ("r,sr", "Sampling rate", cxxopts::value<unsigned>()->default_value("3"))
    //    ("h,help", "Print help");

    options.parse_check(argc, argv);
    
    //if(result.count("help"))
   // {
   //     std::cout << options.help() << std::endl;
   //     exit(0);
   // }

   // if(!result.count("i") || !result.count("o"))
   // {
   //     std::cout << "[ERROR] Please specifiy the input and output file" << std::endl;
   //     std::cout << std::endl;
   //     std::cout << options.help() << std::endl;
   //     exit(0);
    //}

    std::string input_file = options.get<std::string>("input");
    std::string output_file = options.get<std::string>("output");
    double dx = options.get<double>("dx");
    double dy = options.get<double>("dy");
    unsigned u_rate = options.get<unsigned>("sr");

    //if(u_rate > 5)
   // {
   //     std::cout << "[WARN] The sampling rate " << u_rate << " is too high (must <= 5)" << std::endl;
   //     exit(0);
   // }


    std::ofstream fout(output_file, std::ios::binary);

    unsigned char* raw_image;
    unsigned int nx;
    unsigned int ny;

    decode(input_file, &raw_image, &nx, &ny);
    double rate = (double)u_rate;
    
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
                int i1 = ny - i-1;

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
