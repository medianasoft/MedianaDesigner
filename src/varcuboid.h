#pragma once

#include <stdexcept>

class VariableCuboid {
private:
    // size
    unsigned zSize;
    unsigned ySize;
    unsigned xSize;
    // container
    double *body;

public:
    // constructor
    inline VariableCuboid(const unsigned newZSize, const unsigned newYSize, const unsigned newXSize):
        zSize (newZSize),
        ySize (newYSize),
        xSize (newXSize)
    {
        if (zSize == 0 || ySize == 0 || xSize == 0)
            throw std::range_error("Cuboid constructor has 0 size");    // # nocov
        body = new double[zSize * ySize * xSize];
    }
    // destructor
    inline ~VariableCuboid() {
        delete [] body;
    }
    // operator () setter
    inline double& operator()(unsigned z, unsigned y, unsigned x) {
        if (z >= zSize || y >= ySize || x >= xSize)
            throw std::range_error("Cuboid subscript out of bounds");   // # nocov
        return body[z*ySize*xSize + y*xSize + x];
    }
    // operator () getter
    inline double operator()(unsigned z, unsigned y, unsigned x) const {
        if (z >= zSize || y >= ySize || x >= xSize)
            throw std::range_error("Cuboid subscript out of bounds");   // # nocov
        return body[z*ySize*xSize + y*xSize + x];
    }
};
