/**
@file PngWriter.h
@brief Simple PNG writer utility using libpng

Replaces the ImageMagick dependency for PPM to PNG conversion.
*/

#ifndef SRC_THAMESLIB_PNGWRITER_H_
#define SRC_THAMESLIB_PNGWRITER_H_

#include <png.h>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>

/**
@class PngWriter
@brief Utility class for writing PNG files directly without ImageMagick
*/
class PngWriter {
public:
    /**
    @brief Write RGB data to a PNG file

    @param filename Output PNG filename
    @param width Image width in pixels
    @param height Image height in pixels
    @param rgbData Vector of RGB values (size should be width * height * 3)
    @return true on success, false on failure
    */
    static bool writeRGB(const std::string& filename, int width, int height,
                         const std::vector<unsigned char>& rgbData) {
        FILE* fp = fopen(filename.c_str(), "wb");
        if (!fp) {
            return false;
        }

        png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING,
                                                   nullptr, nullptr, nullptr);
        if (!png) {
            fclose(fp);
            return false;
        }

        png_infop info = png_create_info_struct(png);
        if (!info) {
            png_destroy_write_struct(&png, nullptr);
            fclose(fp);
            return false;
        }

        if (setjmp(png_jmpbuf(png))) {
            png_destroy_write_struct(&png, &info);
            fclose(fp);
            return false;
        }

        png_init_io(png, fp);

        // Set image attributes
        png_set_IHDR(png, info, width, height, 8,
                     PNG_COLOR_TYPE_RGB,
                     PNG_INTERLACE_NONE,
                     PNG_COMPRESSION_TYPE_DEFAULT,
                     PNG_FILTER_TYPE_DEFAULT);

        png_write_info(png, info);

        // Write image data row by row
        std::vector<png_bytep> row_pointers(height);
        for (int y = 0; y < height; y++) {
            row_pointers[y] = const_cast<png_bytep>(&rgbData[y * width * 3]);
        }
        png_write_image(png, row_pointers.data());

        png_write_end(png, nullptr);
        png_destroy_write_struct(&png, &info);
        fclose(fp);

        return true;
    }

    /**
    @brief Convert a PPM file to PNG and optionally delete the PPM

    @param ppmFilename Input PPM filename
    @param pngFilename Output PNG filename
    @param deletePpm If true, delete the PPM file after successful conversion
    @return true on success, false on failure
    */
    static bool convertPpmToPng(const std::string& ppmFilename,
                                 const std::string& pngFilename,
                                 bool deletePpm = true) {
        std::ifstream ppmFile(ppmFilename);
        if (!ppmFile.is_open()) {
            return false;
        }

        std::string magic;
        int width, height, maxVal;

        ppmFile >> magic;
        if (magic != "P3") {
            return false;  // Only support ASCII PPM (P3)
        }

        ppmFile >> width >> height >> maxVal;
        if (width <= 0 || height <= 0 || maxVal <= 0) {
            return false;
        }

        std::vector<unsigned char> rgbData(width * height * 3);
        double scale = 255.0 / maxVal;

        for (int i = 0; i < width * height * 3; i++) {
            int val;
            ppmFile >> val;
            rgbData[i] = static_cast<unsigned char>(val * scale);
        }

        ppmFile.close();

        bool success = writeRGB(pngFilename, width, height, rgbData);

        if (success && deletePpm) {
            std::remove(ppmFilename.c_str());
        }

        return success;
    }
};

#endif  // SRC_THAMESLIB_PNGWRITER_H_
