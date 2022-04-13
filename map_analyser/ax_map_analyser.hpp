#ifndef AX_MAP_ANALYSER_HPP
#define AX_MAP_ANALYSER_HPP

#include <iostream>
#include <math.h>
#include <memory.h>

#include <vector>
#include <algorithm>
#include <functional>
#include <list>

#include <ax_core.hpp>

namespace ax {
  namespace man {
    std::function<float(float)> interpol(float a, float b) {
      return [&](float t) -> float { return a + (b - a) * t; };
    }

    std::function<float(float)> interpol(float y0, float y1, float y2) {
      return [&](float t) -> float {
        return y0 + 2 * (y1 - y0) * t + (y2 - 2 * y1 + y0) * t * (2 * t - 1);
      };
    }

    using REAL = float;

    class Interpolant {
    public:
      virtual REAL operator()(REAL t) = 0;

      virtual REAL Reverse(REAL ft) = 0;
    };

    class Interpolant1 : public Interpolant {
    private:
      REAL k = 1;
      REAL b = 0;

    public:
      Interpolant1() = default;

      Interpolant1(REAL f1, REAL f0) :
        k(f0 - f1), b(f0) { }

      void Set(REAL f1, REAL f0) {
        *this = std::move(Interpolant1(f1, f0));
      }

      // t = [-1,0]
      virtual REAL operator()(REAL t) override {
        return k * t + b;
      }

      virtual REAL Reverse(REAL ft) override {
        return (ft - b) / k;
      }
    };

    class Interpolant2 : public Interpolant {
    private:
      REAL a = 1;
      REAL b = 0;
      REAL c = 0;

    public:
      Interpolant2() = default;

      Interpolant2(REAL f2, REAL f1, REAL f0) :
        a((f2 + f0) / 2 - f1), b((f2 + 3 * f0) / 2 - 2 * f1), c(f0) { }

      void Set(REAL f2, REAL f1, REAL f0) {
        *this = std::move(Interpolant2(f2, f1, f0));
      }

      // t = [-2,0]
      virtual REAL operator()(REAL t) override {
        return (a * t + b) * t + c;
      }

      virtual REAL Reverse(REAL ft) override {
        if (a == 0)
          return (ft - c) / b;
        REAL D = b * b - 4 * a * (c - ft);
        if (D < 0)
          return -99999;
        auto roots = sorted_pair((-b - sqrtf(D)) / (2 * a), (-b + sqrtf(D)) / (2 * a));
        if (roots.second <= 0)
          return roots.second;
        return roots.first;
      }
    };

    class Interpolant3 : public Interpolant {
    private:
      REAL a = 1;
      REAL b = 0;
      REAL c = 0;
      REAL d = 0;
      Interpolant1 l21;

    public:
      Interpolant3() = default;

      Interpolant3(REAL f3, REAL f2, REAL f1, REAL f0) :
        a(2 * (f2 - f1) + (f2 - f3) / 1 + (f1 - f2) / 1),
        b(9 * (f2 - f1) + 4 * (f2 - f3) / 1 + 5 * (f1 - f2) / 1),
        c(12 * (f2 - f1) + 5 * (f2 - f3) / 1 + 8 * (f1 - f2) / 1),
        d(5 * f2 - 4 * f1 + 2 * (f2 - f3) / 1 + 4 * (f1 - f2) / 1) {

        l21.Set(f2, f1);
      }

      void Set(REAL f3, REAL f2, REAL f1, REAL f0) {
        *this = std::move(Interpolant3(f3, f2, f1, f0));
      }

      // t = [-2,-1]
      virtual REAL operator()(REAL t) override {
        if (t >= -2 && t <= -1)
          return ((a * t + b) * t + c) * t + d;
        return l21(t + 1);
      }

      virtual REAL Reverse(REAL ft) override {
        return l21.Reverse(ft) - 1;
      }
    };

    enum class CELL_STATE {
      SPARK,
      RISING,
      PEAK,
      FALLING,
      UNDERSHOOT
    };

    enum class PROCESS_STATE {
      IGNITION,
      RELAXING
    };

    const bool NUMBERS[11][5][4] = {
      { // 0
        {0, 1, 1, 0},
        {1, 0, 0, 1},
        {1, 0, 0, 1},
        {1, 0, 0, 1},
        {0, 1, 1, 0}
      },
      { // 1
        {0, 0, 0, 1},
        {0, 0, 1, 1},
        {0, 1, 0, 1},
        {0, 0, 0, 1},
        {0, 0, 0, 1}
      },
      { // 2
        {0, 1, 1, 0},
        {1, 0, 0, 1},
        {0, 0, 1, 0},
        {0, 1, 0, 0},
        {1, 1, 1, 1}
      },
      { // 3
        {0, 1, 1, 0},
        {1, 0, 0, 1},
        {0, 0, 1, 0},
        {1, 0, 0, 1},
        {0, 1, 1, 0}
      },
      { // 4
        {1, 0, 0, 0},
        {1, 0, 1, 0},
        {0, 1, 1, 1},
        {0, 0, 1, 0},
        {0, 0, 1, 0}
      },
      { // 5
        {1, 1, 1, 0},
        {1, 0, 0, 0},
        {1, 1, 1, 0},
        {0, 0, 1, 0},
        {1, 1, 0, 0}
      },
      { // 6
        {0, 1, 1, 0},
        {1, 0, 0, 0},
        {1, 1, 1, 0},
        {1, 0, 0, 1},
        {0, 1, 1, 0}
      },
      { // 7
        {1, 1, 1, 1},
        {0, 0, 0, 1},
        {0, 0, 1, 0},
        {0, 1, 0, 0},
        {1, 0, 0, 0}
      },
      { // 8
        {0, 1, 1, 0},
        {1, 0, 0, 1},
        {0, 1, 1, 0},
        {1, 0, 0, 1},
        {0, 1, 1, 0}
      },
      { // 9
        {0, 1, 1, 0},
        {1, 0, 0, 1},
        {0, 1, 1, 1},
        {0, 0, 0, 1},
        {0, 1, 1, 0}
      },
      { // _
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0},
        {0, 0, 0, 0}
      },
    };

    struct CellInfo {
      static constexpr float inertness = 0.5f;
      Moving_average E_min_max[2] = {inertness, inertness};
      Moving_average derE_min_max[2] = {inertness, inertness};

      float E_norm = 1;
      float derE_norm = 1;
    };

//#define INTER2
//#define INTER3

    template<int buffer_len = 4>
    class Map {
    private:
      FILE *_fmap;
      FILE *_fout;
      int _lenX = 0;
      int _lenY = 0;
      float *_Nframes;
      float *_outframe;
      float *_filterframe;
      CellInfo *types_;
      int _N = 1;
      int _current_frame = 0;
      const float _dt = 0.0001;

      int *leading_area_;
      Interpolant1 *interpols1_;
#ifdef INTER2
      Interpolant2 *interpols2_;
#endif
#ifdef INTER3
      Interpolant3 *interpols3_;
#endif

      Map() = delete;

      template <int prev = 0>
      float &E(int x, int y) {
        return _Nframes[x + _lenX * y + ((_current_frame + prev) % _N) * _lenX * _lenY];
        //return operator()(x, y, prev);
      }

      float flowX(int x, int y) {
        return derivative(x + 1, y) - derivative(x - 1, y);
      }

      float flowY(int x, int y) {
        return derivative(x, y + 1) - derivative(x, y - 1);
      }

      float flowAtoB(Point &a, Point &b) {
        return (E(b.x, b.y) - E<-1>(b.x, b.y)) - (E(a.x, a.y) - E<-1>(a.x, a.y));
      }

      inline int xy(int x, int y) {
        return x + _lenX * y;
      }

      // assuming adaptive phase space calibration with Moving average
      inline void UpdateType(int x, int y) {
        static constexpr int local_extr_weight = 10; // 10%
        int i = xy(x, y);
        CellInfo &t = types_[i];

        // ---- E --- //
        auto EE = E<-1>(x, y);
        if (E<-2>(x, y) > EE && EE < E(x, y)) // (-1) local min
          t.E_min_max[0].update<local_extr_weight>(EE);
        else // repeat previous min
          t.E_min_max[0].update(min(E(x, y), t.E_min_max[0].get()));

        if (E<-2>(x, y) < EE && EE > E(x, y)) // (-1) local max
          t.E_min_max[1].update<local_extr_weight>(EE);
        else // repeat previous max
          t.E_min_max[1].update(max(E(x, y), t.E_min_max[1].get()));

        t.E_norm = t.E_min_max[1].get() - t.E_min_max[0].get();

        // ---- derivative E --- //
        auto derEE = derivative<-1>(x, y);
        if (derivative<-2>(x, y) > derEE && derEE < derivative<0>(x, y)) // (-1) local min
          t.derE_min_max[0].update<local_extr_weight>(derEE);
        else // repeat previous min
          t.derE_min_max[0].update(min(derivative(x, y), t.derE_min_max[0].get()));

        if (derivative<-2>(x, y) < derEE && derEE > derivative<0>(x, y)) // (-1) local max
          t.derE_min_max[1].update<local_extr_weight>(derEE);
        else // repeat previous min
          t.derE_min_max[1].update(max(derivative(x, y), t.derE_min_max[1].get()));

        t.derE_norm = t.derE_min_max[1].get() - t.derE_min_max[0].get();

        interpols1_[i].Set(E<-1>(x, y), E(x, y));
#ifdef INTER2
        interpols2_[i].Set(E(x, y, -2), E(x, y, -1), E(x, y, 0));
#endif
#ifdef INTER3
        interpols3_[i].Set(E(x, y, -3), E(x, y, -2), E(x, y, -1), E(x, y, 0));
#endif
      }

    public:
      constexpr static int SPARK_TRESHOLD = -30;

      Map(const char *filename, int X, int Y) {
        static_assert(buffer_len > 3, "Buffer for MapANalyser is too small (need > 3)");
        std::string name(filename);
        _fmap = fopen(filename, "rb");
        _fout = fopen((name + ".bina").c_str(), "wb");
        _lenX = X;
        _lenY = Y;
        _N = buffer_len;
        _outframe = new float[X * Y];
        _filterframe = new float[X * Y];
        types_ = new CellInfo[X * Y];
        _Nframes = new float[buffer_len * X * Y];

        leading_area_ = new int[X * Y];
        interpols1_ = new Interpolant1[X * Y]();
#ifdef INTER2
        interpols2_ = new Interpolant2[X * Y]();
#endif
#ifdef INTER3
        interpols3_ = new Interpolant3[X * Y]();
#endif

        for (_current_frame = 0; _current_frame < _N; _current_frame++) {
          if (fread(&_Nframes[_current_frame * _lenX * _lenY], sizeof(float), static_cast<size_t>(X * Y), _fmap) == 0)
            return;
        }
        _current_frame--;

        for (int y = 0; y < _lenY; y++) {
          for (int x = 0; x < _lenX; x++) {
            UpdateType(x, y);
            leading_area_[xy(x, y)] = 0;
          }
        }
      }


      double GetTimestep() {
        return _dt;
      }

      ~Map() {
        fclose(_fmap);
        fclose(_fout);
        delete[] _Nframes;
        delete[] _outframe;
        delete[] _filterframe;
        delete[] interpols1_;
#ifdef INTER2
        delete[] interpols2_;
#endif
#ifdef INTER3
        delete[] interpols3_;
#endif
      }

      template <int prev = 0>
      inline float derivative(int x, int y) {
        return (E<prev>(x, y) - E<prev - 1>(x, y)) / _dt;
      }

      void ClearOutFrame() {
        for (int i = 0; i < _lenX * _lenY; i++)
          _outframe[i] = 0;
      }

      void DumpOutFrame(float *to) {
        memcpy(to, _outframe, sizeof(float) * _lenX * _lenX);
      }

      void RestoreOutFrame(float *from) {
        memcpy(_outframe, from, sizeof(float) * _lenX * _lenX);
      }

      template<int L>
      inline void _FilterFrame(const float matrix[L][L], float *aim) {
        static_assert(L % 2 == 1, "Filter matrix must be odd size");
        constexpr int R = L / 2;
        for (int i = 0; i < _lenX * _lenY; i++) {
          _filterframe[i] = 0;
        }
        for (int y0 = R; y0 < _lenY - R; y0++) {
          for (int x0 = R; x0 < _lenX - R; x0++) {
            for (int y = y0 - L / 2; y <= y0 + L / 2; y++) {
              for (int x = x0 - L / 2; x <= x0 + L / 2; x++) {
                int mi = x - x0 + L / 2;
                int mj = y - y0 + L / 2;
                _filterframe[xy(x0, y0)] += aim[xy(x, y)]* matrix[mi][mj];
              }
            }
          }
        }
        memcpy(aim, _filterframe, sizeof(float) * _lenX * _lenY);
      }

      template<int L>
      void FilterOutFrame(const float matrix[L][L], int N = 1) {
        for (int i = 0; i < N; ++i)
          _FilterFrame<L>(matrix, _outframe);
      }

      inline void Write_Line_Out_Frame(std::pair<int, int> r0,
                                std::pair<REAL, REAL> dir,
                                std::pair<REAL, REAL> range, int L) {
        if (L < 2)
          return;
        int &x0 = r0.first;
        int &y0 = r0.second;
        REAL &lx = dir.first;
        REAL &ly = dir.second;

        if ((x0 - L) < 0 ||
            (x0 + L) >= _lenX ||
            (y0 - L) < 0 ||
            (y0 + L) >= _lenY)
          return;

        int x1 = static_cast<int>(L * lx) + x0;
        int y1 = static_cast<int>(L * ly) + y0;
        REAL x = (3 * x0 - x1) / 2;
        REAL y = (3 * y0 - y1) / 2;
        REAL color = range.first;
        REAL d_color = (range.second - range.first) / L;

        for (int i = 0; i < L; ++i) {
          int xi = static_cast<int>(x += lx);
          int yi = static_cast<int>(y += ly);
          _outframe[xy(xi, yi)] = (color += d_color);
        }
        /*/! Arrow
        _outframe[xy(x1 - 1, y1)] =
        _outframe[xy(x1 + 1, y1)] =
        _outframe[xy(x1, y1 - 1)] =
        _outframe[xy(x1, y1 + 1)] = range.second;
        //*/
      }

      template <int mult = 100>
      inline void Draw_Digit_Out_Frame (std::pair<int, int> r0, int num) {
        for (int j = 0; j < 5; ++j) {
          for (int i = 0; i < 4; ++i) {
            if (NUMBERS[num][j][i])
              Eout(r0.first + i - 2, r0.second - j - 2) = NUMBERS[num][j][i] * mult;
          }
        }
      }

      template <int mult = 100>
      inline void Draw_Number_Out_Frame (std::pair<int, int> r0, int num) {
        static const int width = 5;

        if (num == 0) {
          Draw_Digit_Out_Frame<mult>(r0, num);
          return;
        }

        int n = static_cast<int>(log10f(num) + 1.f);
        int offset = n * width / 2;
        r0.first += offset;
        for (int i = 0; i < n; ++i) {
          int d = num % 10;
          Draw_Digit_Out_Frame<mult>(r0, d);
          r0.first -= width;
          num /= 10;
        }
      }

      void NormalizeOutFrame(REAL norm, bool integrate = false) {
        REAL integral = 0;
        REAL min = _outframe[0];
        REAL max = _outframe[0];
        for (int i = 0; i < _lenX * _lenY; ++i) {
          auto x = _outframe[i];
          integral += x;
          if (x > max)
            max = x;
          else if (x < min)
            min = x;
        }
        REAL range = max - min;
        for (int i = 0; i < _lenX * _lenY; ++i) {
          if (!integrate)
            _outframe[i] = norm * (_outframe[i] - min) / range;
          else
            _outframe[i] = _outframe[i] / integral;
        }
      }

      void WriteOutFrame(float *from = nullptr) {
        float *source = (from == nullptr) ? (_outframe) : (from);
        fwrite(source, sizeof(float), static_cast<size_t>(_lenX * _lenY), _fout);
        fflush(_fout);
      }

      inline float &Eout(int x, int y) {
        return _outframe[x + _lenX * y];
      }

      int getX() const {
        return _lenX;
      }

      int getY() const {
        return _lenY;
      }

      int getFrame() const {
        return _current_frame + 1;
      }

      float getTime() const {
        return _current_frame * _dt;
      }

      float *GetCurrentFrame() const {
        return &_Nframes[_current_frame % _N];
      }

      /*
      inline float &operator()(int x, int y, int prev = 0) {
        return _Nframes[x + _lenX * y + ((_current_frame + prev) % _N) * _lenX * _lenY];
      }
       */

      bool NextFrame() {

        // FOR DEBUG
        //int num = (_current_frame + 1) % _N;
        //printf("num = %d, _current_frame = %d\n", num, _current_frame);

        bool res = fread(&_Nframes[((++_current_frame) % _N) * _lenX * _lenY],
                         sizeof(float),
                         static_cast<size_t>(_lenX * _lenY), _fmap) == static_cast<size_t>(_lenX * _lenY);   // reading next frames from file to buffer
        // if success
        if (res) {       
          for (int y = 0; y < _lenY; y++) {
            for (int x = 0; x < _lenX; x++) {
              UpdateType(x, y);
              leading_area_[xy(x, y)] = 0;
            }
          }
        }
        return res;
      }

      template <int prev = 0>
      inline float getPhase(int x, int y) {
        CellInfo &type = types_[xy(x, y)];
        float E0 = (type.E_min_max[1].get() + type.E_min_max[0].get()) / 2;
        float derE0 = (type.derE_min_max[1].get() + type.derE_min_max[0].get()) / 2;
        return atan2f((derivative<prev>(x, y) - derE0) / type.derE_norm, (E<prev>(x, y) - E0) / type.E_norm);
      }

      bool isFront(int x, int y) {
        return getPhase(x, y) > (M_PI * 0.70f);
      }

      bool isRefrac(int x, int y) {
        return getPhase(x, y) < (-M_PI * 0.70f);
      }

      bool isHill(int x, int y) {
        return !isFront(x, y) && !isRefrac(x, y);
      }

      float Jump(int x, int y) {
        return E<0>(x, y) - E<-1>(x, y);
      }

      bool isSpark(int x, int y) {
        return E<-1>(x, y) < SPARK_TRESHOLD && E<0>(x, y) > SPARK_TRESHOLD;
      }

      CELL_STATE getState(int x, int y) {
        float a = E(x, y, -2);
        float b = E(x, y, -1);
        float c = E(x, y, 0);
        if (isSpark(x, y))
          return CELL_STATE::SPARK;
        if (b > SPARK_TRESHOLD && c > SPARK_TRESHOLD) {
          if (b < c)
            return CELL_STATE::RISING;
          if (a < b)
            return CELL_STATE::PEAK;
          return CELL_STATE::FALLING;
        }
        return CELL_STATE::UNDERSHOOT;
      }

      bool isEdgeState(int x, int y, int radius, CELL_STATE state) {
        int i, j;
        j = y - radius;
        for (i = x - radius + 1; i < x + radius; i++)
          if (getState(i, j) != state)
            return false;
        j = y + radius;
        for (i = x - radius + 1; i < x + radius; i++)
          if (getState(i, j) != state)
            return false;
        i = x - radius;
        for (j = y - radius + 1; j < y + radius; j++)
          if (getState(i, j) != state)
            return false;
        i = x + radius;
        for (j = y - radius + 1; j < y + radius; j++)
          if (getState(i, j) != state)
            return false;
        return true;
      }

      bool isEdgeState(int x, int y, int radius_from, int radius_to, CELL_STATE state) {
        for (int r = radius_from; r <= radius_to; r++)
          if (!isEdgeState(x, y, r, state))
            return false;
        return true;
      }

      bool isAreaState(int x, int y, int radius, CELL_STATE state) {
        return isEdgeState(x, y, 1, radius, state);
      }

      template<int R, int JT, int ST = SPARK_TRESHOLD>
      bool _isLeadingCentre(int x0, int y0) {
        const float JUMP_THRESHOLD = JT;
        //const int R = 3;
        if (!isSpark(x0, y0))
          return false;
        SphereArea<1, R> ring(x0, y0);
        float E0 = E(x0, y0);
        float J0 = Jump(x0, y0);
        Interpolant1 i(E(x0, y0, -1), E(x0, y0)); // 2
        float ts = i.Reverse(ST);
        if (J0 < JUMP_THRESHOLD)
          return false;
        for (auto &&p : ring.points()) {
          int x = p.x;
          int y = p.y;
          auto i = Interpolant2(E(x, y, -2), E(x, y, -1), E(x, y)); // 3
          if (i(ts) > ST)
            return false;
          if (E(x, y) > E0)
            return false;
          //if(Jump(x, y) > J0)
          //  return false;
        }
        return true;
      }

      bool isLeadingCentre(int x0, int y0) {
        bool a = _isLeadingCentre<4, 50, SPARK_TRESHOLD>(x0, y0);
        bool b = _isLeadingCentre<3, 53, -20>(x0, y0);
        bool c = _isLeadingCentre<2, 60, -10>(x0, y0);
        // ACH fix
        bool aA = _isLeadingCentre<4, 30, SPARK_TRESHOLD>(x0, y0);
        bool bA = _isLeadingCentre<3, 34, -20>(x0, y0);
        bool cA = _isLeadingCentre<2, 40, -10>(x0, y0);
        return (a && b && c) || (aA && bA && cA);
      }

      template<int R = 4, int ST = SPARK_TRESHOLD, int MINDELTA = 1>
      bool _isSmallLeadingCentre(int x0, int y0) {
        if (E(x0, y0, 0) < ST || E(x0, y0, -1) > ST)
          return false;
        if (E(x0, y0, -1) < ST || E(x0, y0, 0) > ST)
          return true;

        if (leading_area_[xy(x0, y0)])
          return false;


        
        SphereArea<1, R> area(x0, y0);
        float ts = interpols1_[xy(x0, y0)].Reverse(ST); // time of spark (x0,y0)
        for (auto &&p : area.points()) {
          int x = p.x;
          int y = p.y;

          if (leading_area_[xy(x, y)])
            return false;
          if (interpols1_[xy(x, y)](ts) > ST) {
            return false;
          }
        }

        
        SphereArea<R, R> ring(x0, y0);
        for (auto &&p : ring.points()) {
          if (ST - interpols1_[xy(p.x, p.y)](ts) < MINDELTA) {
            return false;
          }
        }
        
        /*
        // detect
        SphereArea<1, R + 2> big_ring(x0, y0);
        for (auto &&p : big_ring.points()) {
          leading_area_[xy(p.x, p.y)] = true;
        }*/
        return true;
      }

      template<int R = 9, int ST = SPARK_TRESHOLD, int MINDELTA = 2>
      bool _isBigLeadingCentre(int x0, int y0) {
        if (E(x0, y0, 0) < ST || E(x0, y0, -1) > ST)
          return false;
        if (leading_area_[xy(x0, y0)])
          return false;

        float ts = interpols1_[xy(x0, y0)].Reverse(ST); // time of spark (x0,y0)

        SphereArea<1, R> area(x0, y0);
        for (auto &&p : area.points()) {
          int x = p.x;
          int y = p.y;


          if (leading_area_[xy(x, y)])
            return false;
          if (interpols1_[xy(x, y)](ts) > ST) {
            return false;
          }
        }

        /*
        SphereArea<R, R> ring(x0, y0);
        for (auto &&p : ring.points()) {
          if (ST - interpols1_[xy(p.x, p.y)](ts) < MINDELTA) {
            return false;
          }
        }
        */

        // detect
        SphereArea<1, R> big_ring(x0, y0);
        for (auto &&p : big_ring.points()) {
          leading_area_[xy(p.x, p.y)] = 1;
        }
        return true;
      }



      template<int R, int ST = SPARK_TRESHOLD>
      bool _isSimpleLeadingCentre(int x0, int y0) 
      {
        
        int checkkk = 0;
        if ((x0 == 135) && (y0 == 141) && (E<0>(x0, y0) > ST) && (E<-1>(x0, y0) < ST))
        {
          checkkk = 1;
          printf("\nPoint of check: (%d, %d), E(-1) = %.2f, E(0) = %.2f; belongToLeadingArea = %d\n", x0, y0, E<-1>(x0, y0), E<0>(x0, y0), leading_area_[xy(x0, y0)]);
        }

        if (E<0>(x0, y0) < ST || E<-1>(x0, y0) > ST)
          return false;
        if (leading_area_[xy(x0, y0)] == 1)
          return false;
          //*/ // UNCOMMENT

        SphereArea<0, R> area(x0, y0);
        //printf("LC point: (%d, %d)\n", x0, y0);
        
        //DEBUG
        float ts = interpols1_[xy(x0, y0)].Reverse(ST);// time of spark (x0,y0)
        double E_tmp;
        //std::cout << "ts = " << ts << "\n";
        int numCells = 0;
        // checking points around, in cirle with radius R
        for (auto &&p : area.points_raw())
        {
          area.Correct(p);
          int x = p.x;
          int y = p.y;

          if (x == x0 && y == y0)
              continue;
            /*
            if(E<0>(x, y) > E<0>(x0, y0))
              return false;
            */

            // only if interpols used
            E_tmp = interpols1_[xy(x, y)](ts);
            //printf("check = %d\n", checkkk);
            numCells++;
           
            /*if ((x0 == 135) && (y0 == 141) && (E<0>(x0, y0) > ST) && (E<-1>(x0, y0) < ST))
            {
                printf("frame# = %d, Points around: (%d, %d), E_interpol = %.2f, E(-1) = %.2f, E(0) = %.2f, belongToLeadingArea = %d, numCells=%d\n", (*this).getFrame(), x, y, E_tmp, E<-1>(x, y), E<0>(x, y), leading_area_[xy(x, y)], numCells);

            }*/
            if (E_tmp > ST)
            {
                //if (checkkk==1)
                //  printf("More than -30mv: (%d,%d)\n", x, y);
                return false;
            }
            if (leading_area_[xy(x, y)] == 1)
            {
              //if (checkkk == 1)
              //  printf("xLead, yLead: (%d,%d)\n", x, y);
              return false;
            }
            //printf("cellscheck: (%d, %d)\n", x, y);
            //if (((*this).getFrame() != 2910))
            //{ 
            //}
        }
        //~area();

        
        // mark new area of LC, containing R cells
        SphereArea<0, 0> big_ring(x0, y0);
        for (auto &&p : big_ring.points_raw()) {
          
          big_ring.Correct(p);

          if (checkkk == 1)
           // printf("lead area coords = (%d, %d)\n", p.x, p.y);
          

          leading_area_[xy(p.x, p.y)] = 1;
        }
        //~big_ring();
        //printf("Cell (%d, %d), is LC\n", x0, y0);
        return true;
      }


      bool isLeadingCentreNew(int x0, int y0) {
        //if (_isSmallLeadingCentre(x0, y0)) // small
        //  return true;
        //return _isBigLeadingCentre(x0, y0); // big
        return _isSimpleLeadingCentre<3>(x0, y0);
      }

      template<int R>
      int _isSingularity(int x0, int y0, int prev = 0) {
        SphereArea<R, R> ring(x0, y0);
        auto array = ring.points_sorted();
        int circ = 0;
        int sign = 0;
        Point *a, *b;
        a = &array.back();
        for (auto &&i : array) {
          b = &i;
          // reverse phase
          //if (fmodf(getPhase(a->x, a->y, prev) - getPhase(a->x, a->y, prev - 1), 2 * M_PI) > 0)
          //  return false;
          // crossing barrier (-Pi) -> (Pi)
          float delta = getPhase(a->x, a->y, prev) - getPhase(b->x, b->y, prev);
          bool act = fabsf(delta) > (1.5f * M_PI);
          circ += act;

          if (!act) { // define orientation (rotor sign)
            sign += delta / fabsf(delta);
          }

          a = b;
        }
        //return circ % 2 == 1;
        if (circ == 1)
          return sign / abs(sign);
        return 0;
      }

      template<int R = 4>
      int isSingularity(int x0, int y0) {
        return (_isSingularity<R>(x0, y0, 0) + _isSingularity<R - 1>(x0, y0, -1)) / 2;
      }

      template <int prev = 0>
      std::pair<REAL, REAL> getGradient(int x0, int y0, REAL *value) {
        constexpr static int order = 2;
        if ((x0 - order) < 0 ||
            (x0 + order) >= _lenX ||
            (y0 - order) < 0 ||
            (y0 + order) >= _lenY) {
          if (value)
            *value = 0;
          return {0, 0};
        }

        /*
        REAL Gx = (-E(x0 - 3, y0, prev) + 9 * E(x0 - 2, y0, prev) - 45 * E(x0 - 1, y0, prev)
                   + 45 * E(x0 + 1, y0, prev) - 9 * E(x0 + 2, y0, prev) + E(x0 + 3, y0, prev)) / 60.f;
        REAL Gy = (-E(x0, y0 - 3, prev) + 9 * E(x0, y0 - 2, prev) - 45 * E(x0, y0 - 1, prev)
                   + 45 * E(x0, y0 + 1, prev) - 9 * E(x0, y0 + 2, prev) + E(x0, y0 + 3, prev)) / 60.f;
        */
        /*
        REAL Gx = (-E(x0 + 2, y0, prev) + 8 * E(x0 + 1, y0, prev)
                   - 8 * E(x0 - 1, y0, prev) + E(x0 - 2, y0, prev)) / 12;
        REAL Gy = (-E(x0, y0 + 2, prev) + 8 * E(x0, y0 + 1, prev)
                   - 8 * E(x0, y0 - 1, prev) + E(x0, y0 - 2, prev)) / 12;
        */

        float Ex[5] =
          {
            (2 * E<prev>(x0 - 2, y0) + E<prev>(x0 - 2, y0 + 1) + E<prev>(x0 - 2, y0 - 1)) / 4,
            (3 * E<prev>(x0 - 1, y0) + E<prev>(x0 - 1, y0 + 1) + E<prev>(x0 - 1, y0 - 1)) / 5,
            0,
            (3 * E<prev>(x0 + 1, y0) + E<prev>(x0 + 1, y0 + 1) + E<prev>(x0 + 1, y0 - 1)) / 5,
            (2 * E<prev>(x0 + 2, y0) + E<prev>(x0 + 2, y0 + 1) + E<prev>(x0 + 2, y0 - 1)) / 4
          };
        float Ey[5] =
          {
            (2 * E<prev>(x0, y0 - 2) + E<prev>(x0 + 1, y0 - 2) + E<prev>(x0 - 1, y0 - 2)) / 4,
            (3 * E<prev>(x0, y0 - 1) + E<prev>(x0 + 1, y0 - 1) + E<prev>(x0 - 1, y0 - 1)) / 5,
            0,
            (3 * E<prev>(x0, y0 + 1) + E<prev>(x0 + 1, y0 + 1) + E<prev>(x0 - 1, y0 + 1)) / 5,
            (2 * E<prev>(x0, y0 + 2) + E<prev>(x0 + 1, y0 + 2) + E<prev>(x0 - 1, y0 + 2)) / 4,
          };

        REAL Gx = (-Ex[4] + 8 * Ex[3] - 8 * Ex[1] + Ex[0]) / 12;
        REAL Gy = (-Ey[4] + 8 * Ey[3] - 8 * Ey[1] + Ey[0]) / 12;

        REAL mod = static_cast<REAL>(sqrt(sq(Gx) + sq(Gy)));
        if (value)
          *value = mod;
        if (mod == 0)
          return {0, 0};
        return {Gx / mod, Gy / mod};
      }

      template<int THRESHOLD = SPARK_TRESHOLD, int R = 20>
      REAL getWaveSpeed(int x0, int y0) {
        constexpr static int steps_limit = 20;

        constexpr static int prev = -2;
        auto isSpark_ = [&](int x, int y) {
          return (E<prev - 1>(x, y) < THRESHOLD) && (E<prev>(x, y) > THRESHOLD);
        };

        if (!isSpark_(x0, y0))
          return -1;

        REAL gr_mod;
        auto dir = getGradient<-2>(x0, y0, &gr_mod);
        if (gr_mod == 0)
          return -1;

        dir.first = -dir.first;
        dir.second = -dir.second;

        int steps = abs(static_cast<int>(R / max(dir.first, dir.second)));
        if (steps > steps_limit)
          return -1;

        REAL x(x0), y(y0);
        for (int i = 0; i < steps; ++i) {
          x += dir.first;
          y += dir.second;

          int ix = static_cast<int>(x);
          int iy = static_cast<int>(y);

          if (ix == x0 && iy == y0)
            continue;

          if (E(ix, iy) < THRESHOLD) {
            return static_cast<REAL>(sqrt(sq(x - x0) + sq(y - y0)) / (2 * _dt));
          }
        }
        return -1;
      }
    };
  }
}

#endif //AX_MAP_ANALYSER_HPP
