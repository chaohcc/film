// Copyright (c) Microsoft Corporation.
// Licensed under the MIT license.
// from ALEX
// the github address of ALEX:

// Zipf generator, inspired by
// https://github.com/brianfrankcooper/YCSB/blob/master/core/src/main/java/site/ycsb/generator/ScrambledZipfianGenerator.java
// https://github.com/brianfrankcooper/YCSB/blob/master/core/src/main/java/site/ycsb/generator/ZipfianGenerator.java

#ifndef ZIPF_H
#define ZIPF_H

#include <random>

class ScrambledZipfianGenerator {
 public:
  static constexpr double ZETAN = 26.46902820178302;
  static constexpr double ZIPFIAN_CONSTANT = 0.99;

  int num_keys_;
  double alpha_;
  double eta_;
  std::mt19937_64 gen_;
  std::uniform_real_distribution<double> dis_;

  explicit ScrambledZipfianGenerator(int num_keys)
      : num_keys_(num_keys), gen_(std::random_device{}()), dis_(0, 1) {
    double zeta2theta = zeta(2);
    alpha_ = 1. / (1. - ZIPFIAN_CONSTANT);
    eta_ = (1 - std::pow(2. / num_keys_, 1 - ZIPFIAN_CONSTANT)) /
           (1 - zeta2theta / ZETAN);
  }

  int nextValue() {
    double u = dis_(gen_);
    double uz = u * ZETAN;

    int ret;
    if (uz < 1.0) {
      ret = 0;
    } else if (uz < 1.0 + std::pow(0.5, ZIPFIAN_CONSTANT)) {
      ret = 1;
    } else {
      ret = (int)(num_keys_ * std::pow(eta_ * u - eta_ + 1, alpha_));
    }

    ret = fnv1a(ret) % num_keys_;
    return ret;
  }

  double zeta(long n) {
    double sum = 0.0;
    for (long i = 0; i < n; i++) {
      sum += 1 / std::pow(i + 1, ZIPFIAN_CONSTANT);
    }
    return sum;
  }

  // FNV hash from https://create.stephan-brumme.com/fnv-hash/
  static const uint32_t PRIME = 0x01000193;  //   16777619
  static const uint32_t SEED = 0x811C9DC5;   // 2166136261
  /// hash a single byte
  inline uint32_t fnv1a(unsigned char oneByte, uint32_t hash = SEED) {
    return (oneByte ^ hash) * PRIME;
  }
  /// hash a 32 bit integer (four bytes)
  inline uint32_t fnv1a(int fourBytes, uint32_t hash = SEED) {
    const unsigned char* ptr = (const unsigned char*)&fourBytes;
    hash = fnv1a(*ptr++, hash);
    hash = fnv1a(*ptr++, hash);
    hash = fnv1a(*ptr++, hash);
    return fnv1a(*ptr, hash);
  }
};


//
// Created by CCMa on 2022/7/19.
//

#ifndef ZIPFGENRATION_SERIESZIPF_H
#define ZIPFGENRATION_SERIESZIPF_H
//
// Created by CCMa on 2022/7/19.
//

#include <algorithm>
#include <cmath>
#include <random>

/** Zipf-like random distribution.
 *
 * "Rejection-inversion to generate variates from monotone discrete
 * distributions", Wolfgang Hörmann and Gerhard Derflinger
 * ACM TOMACS 6.3 (1996): 169-184
 */
template<class IntType = unsigned long, class RealType = double>
class zipf_distribution
{
public:
    typedef RealType input_type;
    typedef IntType result_type;

    static_assert(std::numeric_limits<IntType>::is_integer, "");
    static_assert(!std::numeric_limits<RealType>::is_integer, "");

    zipf_distribution(const IntType n=std::numeric_limits<IntType>::max(),
                      const RealType q=1.0)
            : n(n)
            , q(q)
            , H_x1(H(1.5) - 1.0)
            , H_n(H(n + 0.5))
            , dist(H_x1, H_n)
    {}

    IntType operator()(std::mt19937& rng)
    {
        while (true) {
            const RealType u = dist(rng);
            const RealType x = H_inv(u);
            const IntType  k = clamp<IntType>(std::round(x), 1, n);
            if (u >= H(k + 0.5) - h(k)) {
                return k;
            }
        }
    }

private:
    /** Clamp x to [min, max]. */
    template<typename T>
    static constexpr T clamp(const T x, const T min, const T max)
    {
        return std::max(min, std::min(max, x));
    }

    /** exp(x) - 1 / x */
    static double
    expxm1bx(const double x)
    {
        return (std::abs(x) > epsilon)
               ? std::expm1(x) / x
               : (1.0 + x/2.0 * (1.0 + x/3.0 * (1.0 + x/4.0)));
    }

    /** H(x) = log(x) if q == 1, (x^(1-q) - 1)/(1 - q) otherwise.
     * H(x) is an integral of h(x).
     *
     * Note the numerator is one less than in the paper order to work with all
     * positive q.
     */
    const RealType H(const RealType x)
    {
        const RealType log_x = std::log(x);
        return expxm1bx((1.0 - q) * log_x) * log_x;
    }

    /** log(1 + x) / x */
    static RealType
    log1pxbx(const RealType x)
    {
        return (std::abs(x) > epsilon)
               ? std::log1p(x) / x
               : 1.0 - x * ((1/2.0) - x * ((1/3.0) - x * (1/4.0)));
    }

    /** The inverse function of H(x) */
    const RealType H_inv(const RealType x)
    {
        const RealType t = std::max(-1.0, x * (1.0 - q));
        return std::exp(log1pxbx(t) * x);
    }

    /** That hat function h(x) = 1 / (x ^ q) */
    const RealType h(const RealType x)
    {
        return std::exp(-q * std::log(x));
    }

    static constexpr RealType epsilon = 1e-8;

    IntType                                  n;     ///< Number of elements
    RealType                                 q;     ///< Exponent
    RealType                                 H_x1;  ///< H(x_1)
    RealType                                 H_n;   ///< H(n)
    std::uniform_real_distribution<RealType> dist;  ///< [H(x_1), H(n)]
};

#endif //ZIPFGENRATION_SERIESZIPF_H


#endif