#ifndef AX_MAP_GENERATOR_HPP
#define AX_MAP_GENERATOR_HPP

#include <iostream>
#include <math.h>
#include <memory.h>

#include <string>
#include <functional>

#include <ax_core.hpp>

namespace ax {
  using REAL = float;

  // MAp Generator
  namespace mag {
    REAL Cell_SAN(REAL type) {
      return type;
    }

    REAL Cell_AT() {
      return -1;
    }

    class Map {
    private:
      FILE *out_OLD_;
      FILE *out_BIN_;
      FILE *out_MAP_;
      FILE *out_MAP_VTK;
      int X_, Y_;
      float *map_;
      float *filter_;

      Map() = delete;

      Map(const Map &) = delete;

      Map &operator=(const Map &) = delete;

      inline int xy(int x, int y) {
        return x + X_ * y;
      }

      inline std::pair<int, int> idx(int i) {
        return {i % X_, i / X_};
      }

      std::string format(REAL type) {
        constexpr static char fmt_pos[] = " %.6f";
        constexpr static char fmt_neg[] = "%.6f";
        char buf[10] = {'\0'};
        sprintf(buf, (type < 0) ? (fmt_neg) : (fmt_pos), fmodf(type, 10));
        std::string t(buf);
        return t;
      }

      /*/!
      constexpr const char* fmt(REAL type) {
        constexpr static const char fmt_pos[] = " %.6f";
        constexpr static const char fmt_neg[] = "%.6f";
        return (type < 0) ? (fmt_neg) : (fmt_pos);
      }
      //*/

    public:

      Map(int X, int Y, const char *name) : X_(X), Y_(Y) {
        std::string fname(name);
        out_OLD_ = fopen(fname.c_str(), "w");
        out_BIN_ = fopen((fname + ".bin").c_str(), "wb");
        out_MAP_ = fopen((fname + ".map").c_str(), "w");
        // added
        out_MAP_VTK = fopen((fname +"_viz.vtk").c_str(), "w");

        map_ = new REAL[X * Y];
        filter_ = new REAL[X * Y];
        for (int i = 0; i < X * Y; ++i)
          map_[i] = filter_[i] = 0;
      }

      ~Map() {
        fclose(out_OLD_);
        fclose(out_BIN_);
        fclose(out_MAP_);

        delete[] map_;
        delete[] filter_;
      }

      void Create(std::function<REAL(int, int)> f) {
        for (int i = 0; i < X_ * Y_; i++) {
          auto coords = idx(i);
          map_[i] = f(coords.first, coords.second);
        }
      }

      template<int L>
      inline void Filter(const float matrix[L][L]) {
        static_assert(L % 2 == 1, "Filter matrix must be odd size");
        constexpr int R = L / 2;

        for (int i = 0; i < X_ * Y_; ++i)
          filter_[i] = 0;

        for (int y0 = R; y0 < Y_ - R; y0++) {
          for (int x0 = R; x0 < X_ - R; x0++) {
            for (int y = y0 - L / 2; y <= y0 + L / 2; y++) {
              for (int x = x0 - L / 2; x <= x0 + L / 2; x++) {
                int mi = x - x0 + L / 2;
                int mj = y - y0 + L / 2;
                filter_[xy(x0, y0)] += map_[xy(x, y)] * matrix[mi][mj];
              }
            }
          }
        }
        memcpy(map_, filter_, sizeof(REAL) * X_ * Y_);
      }

      void Write() {
        constexpr static char fmt_pos[] = " %.6f\t";
        constexpr static char fmt_neg[] = "%.6f\t";

        fwrite(map_, sizeof(REAL), static_cast<size_t>(X_ * Y_), out_BIN_);


        // VTK output
        fprintf(out_MAP_VTK, "# vtk DataFile Version 3.0\nSolution\nASCII\nDATASET RECTILINEAR_GRID\nDIMENSIONS 201 201 1\n");
        fprintf(out_MAP_VTK, "X_COORDINATES 201 double\n");

        //auto coords = idx(i);
        for (int i = 0; i < 201; i++)
          fprintf(out_MAP_VTK, "%d ", i);
        fprintf(out_MAP_VTK, "\n");
        
        fprintf(out_MAP_VTK, "Y_COORDINATES 201 double\n");
        
        for (int i = 0; i < 201; i++)
          fprintf(out_MAP_VTK, "%d ", i);
        fprintf(out_MAP_VTK, "\n");

        fprintf(out_MAP_VTK,"Z_COORDINATES 1 double\n0");
        fprintf(out_MAP_VTK,"CELL_DATA %d\n", X_ * Y_);
        
        
        // printing all the data in the file
        fprintf(out_MAP_VTK, "SCALARS type double\nLOOKUP_TABLE default\n");
        for (int idx = 0; idx < X_ * Y_; idx++) 
        {
          REAL type = map_[idx];
          fprintf(out_MAP_VTK, "%e ", type);

          if (( (idx + 1)%X_ ) == 0)
            fprintf(out_MAP_VTK, "\n");
        }


        //fprintf(out_MAP_VTK, "%d,%d\n", X_, Y_);


        fprintf(out_MAP_, "%d\t%d\n", X_, Y_);
        for (int i = 0; i < X_ * Y_; ++i) {
          REAL type = map_[i];
          auto coords = idx(i);
          fprintf(out_MAP_, "%d\t%d\t%.6f\n", coords.first, coords.second, type);
          
          //// VTK
          //fprintf(out_MAP_, "%d,%d,%.6f\n", coords.first, coords.second, type);

          fprintf(out_OLD_, (type < 0) ? (fmt_neg) : (fmt_pos), type);
        }

        // VTK output

        fflush(out_OLD_);
        fflush(out_BIN_);
        fflush(out_MAP_);
      }
    };
  }
}

#endif //AX_MAP_GENERATOR_HPP