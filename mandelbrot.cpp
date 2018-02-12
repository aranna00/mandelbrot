// mandelbrot.cpp
// compile with: g++ -std=c++11 mandelbrot.cpp -o mandelbrot
// view output with: eog mandelbrot.ppm

#include <fstream>
#include <complex> // if you make use of complex number facilities in C++
#include <iostream>;
#include <thread>;
using namespace std;

template<class T>
struct RGB {
    T r, g, b;
};


template<class T>
class Matrix {
public:
    Matrix(const size_t rows, const size_t cols) : _rows(rows), _cols(cols) {
        _matrix = new T *[rows];
        for (size_t i = 0; i < rows; ++i) {
            _matrix[i] = new T[cols];
        }
    }

    Matrix(const Matrix &m) : _rows(m._rows), _cols(m._cols) {
        _matrix = new T *[m._rows];
        for (size_t i = 0; i < m._rows; ++i) {
            _matrix[i] = new T[m._cols];
            for (size_t j = 0; j < m._cols; ++j) {
                _matrix[i][j] = m._matrix[i][j];
            }
        }
    }

    ~Matrix() {
        for (size_t i = 0; i < _rows; ++i) {
            delete[] _matrix[i];
        }
        delete[] _matrix;
    }

    T *operator[](const size_t nIndex) {
        return _matrix[nIndex];
    }

    size_t width() const { return _cols; }

    size_t height() const { return _rows; }

protected:
    size_t _rows, _cols;
    T **_matrix;
};

// Portable PixMap image
class PPMImage : public Matrix<RGB<unsigned char> > {
public:
    PPMImage(const size_t height, const size_t width) : Matrix(height, width) {}

    void save(const std::string &filename) {
        std::ofstream out(filename, std::ios_base::binary);
        out << "P6" << std::endl << _cols << " " << _rows << std::endl << 255 << std::endl;
        for (size_t y = 0; y < _rows; y++) {
            for (size_t x = 0; x < _cols; x++) {
                out << _matrix[y][x].r << _matrix[y][x].g << _matrix[y][x].b;
            }
            if(y%100==0) {
                std::cout << "Processing image " << static_cast<int>(static_cast<double>(y)/static_cast<double>(_rows)*100) << "% \r" << std::flush;
            }
        }
        std::cout << "Image processed" << std::endl;
    }
};

int main() {
    long PhysMemory = 16;
    cout << "found " << thread::hardware_concurrency() << " threads" << endl;
    const unsigned width = 2560*10;
    const unsigned height = 1080*10;
    double max_iterations = 2500;

    PPMImage image(height, width);

//    while (true){
//        cout << "using too much memory, scaling down picture";
//        delete(image);
//        PPMImage image(height*.99,width*.99);
//    }

    double zoom = pow(2, 64);
    double minReal = (-.32175*zoom-1.17)/zoom*width/height; // move the left border
    double maxReal = (-.31675*zoom+1.21)/zoom*width/height; // move the right border
    double minImag = (0.063*zoom-1.16)/zoom; // moves the bottom side
    double maxImag = (minImag+(maxReal-minReal)*height/width); // moves the top side
    double realFactor = (maxReal-minReal)/(width-1);
    double imagFactor = (maxImag-minImag)/(height-1);


    for (int y = 0; y < height; ++y) {
        double c_im = maxImag-y*imagFactor;
        for (int x = 0; x < width; ++x) {
            double c_re = minReal+x*realFactor;
            double Z_re = c_re, Z_im = c_im; // Set Z = c
            unsigned current_ite;
            for(current_ite=0; current_ite<max_iterations; ++current_ite)
            {
                double Z_re2 = Z_re*Z_re, Z_im2 = Z_im*Z_im;
                if(Z_re2 + Z_im2 > 4)
                {
                    break;
                }
                Z_im = 2*Z_re*Z_im + c_im;
                Z_re = Z_re2 - Z_im2 + c_re;
            }
            if(current_ite<max_iterations/2) {
                image[y][x].g = current_ite/(max_iterations/2-1)*255;
            } else if (current_ite<max_iterations-1){
                image[y][x].g = 255;
                image[y][x].r = image[y][x].b = (current_ite - (max_iterations / 2)) / (max_iterations / 2) * 255;
            }
        }
        if(y%(height/100)==0) {
            std::cout << "Calculating " << static_cast<int>(static_cast<double>(y)/static_cast<double>(height)*100) << "% \r" << std::flush;
        }
    }
    std::cout << "Calculating 100%" << std::endl;

    /*
    image[y][x].r = image[y][x].g = image[y][x].b = 255; // white pixel
    
    image[y][x].r = image[y][x].g = image[y][x][b] = 0; // black pixel
   
    // red pixel
    image[y][x].r = 255;
    image[y][x].g = 0;
    image[y][x].b = 0;
    */

    image.save("R:\\Development\\mandelbrot_zoom.ppm");
}

