#include <arrayfire.h>
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <unistd.h>
#include <chrono>

using namespace af;
using namespace std;

//static const int NSAMPLES = 1e3;
static const int NSAMPLES = 20;
static int WIDTH=1200, HEIGHT=600;

float point_kernel[] = {1,0,0,0,1,
                        0,1,0,1,0,
                        0,0,1,0,0,
                        0,1,0,1,0,
                        1,0,0,0,1};

int main(int argc, char *argv[]) {
    try {
        //// Select a device and display arrayfire info
        int device = argc > 1 ? atoi(argv[1]) : 0;
        int N_order = (argc > 2) ? min(atoi(argv[2]), NSAMPLES) : 4; //set poly-order and make sure system is overdetermined
        af::setDevice(device);
        af::info();
        //initialize arrayfire window
        af::Window myWindow(WIDTH, HEIGHT, "Least Squares Example: ArrayFire");

        while(!myWindow.close()) {
	    //create a set of random particles
	    array particles = randu(NSAMPLES, 2); // <x,y>
	    //particles.col(1) *= particles.col(0)+0.3; // generate something resembling a line
	    particles.col(0) *= WIDTH;
	    particles.col(1) *= HEIGHT;

            //add points in image matrix
            array image = constant(0, HEIGHT, WIDTH, f32);
            array ids(NSAMPLES, u32);
            ids = (particles.col(0).as(u32) * (HEIGHT)) + particles.col(1).as(u32);
            image(ids) = 1;

            array mask_point = array(5, 5, point_kernel);
            image = dilate(image, mask_point);


            //we will fit a line first, so form the regressors(X) matrix
            array X = constant(1, NSAMPLES, 2, f32);
            X.col(1) *= particles.col(0);

            //solve for param vector algebraically
            array b = matmul(matmul(inverse(matmul(X.T(), X)), X.T()), particles.col(1));

            //add line to image
            float *params = b.host<float>();
            array yplot = range(dim4(WIDTH)) * params[1] + params[0];
            delete [] params;
            ids = (range(dim4(WIDTH)).as(u32) * (HEIGHT)) + yplot.as(u32);
            ids = ids(yplot > 0 && yplot < HEIGHT);
            if(ids.dims(0) > 0)
                image(ids) = 1;

            //form regressors(X) matrix for N-order polynomial fit
            X = constant(1, NSAMPLES, N_order, f32);
            X.col(1) *= particles.col(0);
            for(int i=2; i<N_order; ++i)
                X.col(i) = X.col(i-1) * particles.col(0);
            b = matmul(matmul(inverse(matmul(X.T(), X)), X.T()), particles.col(1));

            //solve for param vector using cholesky decomposition
            array U;
            cholesky(U, matmul(X.T(), X), true);
            array z = solve(U.T(), matmul(X.T(), particles.col(1)), AF_MAT_LOWER);
            array bn = solve(U, z, AF_MAT_UPPER);

            array xplot = range(dim4(WIDTH, N_order));
            xplot.col(0) = 1;
            for(int i=2; i<N_order; ++i){
                xplot.col(i) *= xplot.col(i-1);
            }

            //add polynomial to image
            yplot = matmul(xplot, bn);
            ids = (range(dim4(WIDTH)).as(u32) * (HEIGHT)) + yplot.as(u32);
            ids = ids(yplot > 0 && yplot < HEIGHT);
            if(ids.dims(0) > 0)
                image(ids) = 1;

            //make points and lines thicker
            array mask = constant(1, 2, 2);
            image = dilate(image, mask);

            //plot scattered points and fitted function(s)
            myWindow.image(image);

            N_order++;
            if(N_order>6) {N_order=2;}
            usleep(100000);
        }
    } catch (af::exception& e) {
        fprintf(stderr, "%s\n", e.what());
        throw;
    }

    return 0;
}
