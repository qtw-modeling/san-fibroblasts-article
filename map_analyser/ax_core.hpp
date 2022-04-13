#ifndef AX_CORE_HPP
#define AX_CORE_HPP

#include <vector>
#include <algorithm>
#include <functional>
#include <list>
//#include <cstdlib>
//#include <ctime>
#include <time.h>
#include <iostream>
#include <stdlib.h>

namespace ax {
  template<typename T>
  inline std::pair<T, T> sorted_pair(T a, T b) {
    return (a > b) ? std::pair<T, T>(b, a) : std::pair<T, T>(a, b);
  }


  template<typename T>
  inline T max(T x, T y) {
    return (x < y) ? y : x;
  }

  template<typename T>
  inline T min(T x, T y) {
    return (x > y) ? y : x;
  }

  template<typename T>
  inline T abs(T x) {
    return (x < 0) ? -x : x;
  }

  template<typename T>
  inline T sq(T x) {
    return x * x;
  }

  template <typename T = double>
  inline T Get_Normal(T mean, T sigma) {
    static T z[2];
    static bool idx = 1;

    if (idx) {
      T x, y, s;
      do {
        x = static_cast<T>(2) * rand() / RAND_MAX - static_cast<T>(1);
        y = static_cast<T>(2) * rand() / RAND_MAX - static_cast<T>(1);
        s = sq(x) + sq(y);
      } while (s > 1 || s == 0);
      s = static_cast<T>(sqrt((-2 * log(s)) / s));
      z[0] = x * s;
      z[1] = y * s;
    }
    idx = !idx;
    return mean + z[idx] * sigma;
  }

  // mean - dispersion - standard deviation
  template <typename T>
  inline std::tuple<T, T, T> Get_M_D_SD(const std::list<T>& list) {
    int n = static_cast<int>(list.size());
    if (n == 1)
      return std::make_tuple(list.front(), static_cast<T>(0), static_cast<T>(0));
    T mean = 0;
    T dispersion = 0;
    T sd = 0;
    for (auto &&p : list) {
      mean += p;
      dispersion += p * p;
    }
    dispersion = (dispersion - mean * mean / n) / (n - 1);
    sd = static_cast<T>(sqrt(dispersion));
    mean /= n;
    return std::make_tuple(mean, dispersion, sd);
  }

  class Moving_average {
  private:
    float hystory;
    float hystory_importance;
    bool initialized;

    Moving_average() = delete;

  public:
    Moving_average(float h) :
      hystory(0), hystory_importance(h), initialized(false) { }

    Moving_average(const Moving_average &x) = default;

    // weight of new element in percents
    template <int weight = -1>
    inline void update(float value) {
      if (initialized) {
        float hi = (weight == -1) ? (hystory_importance) : (1.f - static_cast<float>(weight) / 100.f);
        hystory = value + hi * (hystory - value);
      } else {
        hystory = value;
        initialized = true;
      }
    }

    inline float get() {
      return hystory;
    }
  };

  class Point {
    template<int r, int R> friend class SphereArea;

  public:
    int x = {-1};
    int y = {-1};

    Point() = default;

    Point(int _x, int _y) : x(_x), y(_y) {
      default_constructor_ = false;
    }

    void set(int _x, int _y) {
      x = _x;
      y = _y;
      default_constructor_ = false;
    }

    operator std::pair<int, int>() {
      return {x, y};
    }

  private:
    bool default_constructor_ = {true};
  };

  inline float r_points(const Point &a, const Point &b) {
    return sqrtf(sq(a.x - b.x) + sq(a.y - b.y));
  }

  template<int r, int R>
  class SphereArea {
  private:
    constexpr static int range_ = 2 * R + 1;
    static bool m[range_][range_];
    static int points_;
    static bool init_;
    static std::array<Point, range_ * range_> points_list_raw_;
    static std::array<Point, range_ * range_> points_list_sorted_;
    int x0_;
    int y0_;

    SphereArea() = delete;
    SphereArea& operator=(SphereArea const&) = delete;

  public:
    SphereArea(int x0, int y0) {
      if (!init_) {
        static_assert((r >= 0 && r <= R), "Wrong 'r' or/and 'R'");
        int c0 = R;
        for (int i = 0; i < range_; i++) {
          for (int j = 0; j < range_; j++) {
            m[i][j] = ((sq(i - c0) + sq(j - c0)) >= sq(r)) && ((sq(i - c0) + sq(j - c0)) < sq(R + 1));
            points_ += m[i][j];
            if (m[i][j])
              points_list_raw_[points_ - 1].set(i, j);
          }
        }
        printf("\nCheck matrix: points = %d\n", points_);
        for (int i = 0; i < range_; i++) {
          for (int j = 0; j < range_; j++) {
            printf("%s", (m[i][j]) ? ("#") : ("."));
          }
          printf("\n");
        }
        printf("\n");
        points_list_sorted_ = points_list_raw_;
        auto end_it = points_list_sorted_.end();
        for (int i = 0; i < (range_ * range_ - points_); ++i)
          end_it--;
        std::sort(points_list_sorted_.begin(),
                  end_it,
                  [](const Point &a, const Point &b) -> bool {

                    return atan2f(a.y - R, a.x - R) < atan2f(b.y - R, b.x - R);
                  });
        init_ = true;

        //! Print points
        printf("{");
        for (auto &&p : points_list_raw_) {
          printf("(%d|%d)", p.x, p.y);
        }
        printf("}\n");
        printf("{");
        for (auto &&p : points_list_sorted_) {
          printf("(%d|%d)", p.x, p.y);
        }
        printf("}\n");
        //*/
      }
      x0_ = x0;
      y0_ = y0;
    }

    bool contains(int x, int y) const {
      if (abs(x - x0_) > R || abs(y - y0_) > R)
        return false;
      return m[x - x0_ + R][y - y0_ + R];
    }

    const static int N() {
      return points_;
    }

    const static size_t range() {
      return static_cast<size_t>(range_);
    }

    inline std::vector<Point> points_sorted() const {
      return std::vector<Point>(points_list_sorted_.data(), points_list_sorted_.data() + points_);
    }

    std::vector<Point> points_raw() const {
      return std::vector<Point>(points_list_raw_.data(), points_list_raw_.data() + points_);
    }

    inline void Correct(Point& p) {
      p.x += x0_ - R;
      p.y += y0_ - R;
    }
  };

  template<int r, int R>
  bool SphereArea<r, R>::init_ = {false};

  template<int r, int R>
  int SphereArea<r, R>::points_ = {0};

  template<int r, int R>
  bool SphereArea<r, R>::m[SphereArea<r, R>::range_][SphereArea<r, R>::range_] = {false};

  //
  template<int r, int R>
  std::array<Point, (ax::SphereArea<r, R>::range_ * ax::SphereArea<r, R>::range_)>
    SphereArea<r, R>::points_list_raw_;

  template<int r, int R>
  std::array<Point, (ax::SphereArea<r, R>::range_ * ax::SphereArea<r, R>::range_)>
    SphereArea<r, R>::points_list_sorted_;
   //*/

  const float FILTER_mean3x3[3][3] = {
    {1.f / 9, 1.f / 9, 1.f / 9},
    {1.f / 9, 1.f / 9, 1.f / 9},
    {1.f / 9, 1.f / 9, 1.f / 9}
  };


  const float FILTER_mean5x5[5][5] = {
    {1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25},
    {1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25},
    {1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25},
    {1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25},
    {1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25, 1.f / 25}
  };

  const float FILTER_sharp[3][3] = {
    {-1.f, -1.f, -1.f},
    {-1.f, 9.f,  -1.f},
    {-1.f, -1.f, -1.f}
  };
}

#endif //AX_CORE_HPP
