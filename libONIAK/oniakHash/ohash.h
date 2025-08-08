#pragma once

#define wyhash_final_version_3

#include "../global.h"
#include <wyhash/wyhash.h>
#include <wyhash/wyhash32.h>


namespace ONIAK {

enum class OHashBackend {
  WYHash,
  WYHash64,
  Unknown
};

template <OHashBackend Backend>
class OHash {
public:
  template <typename T>
  OHash(T t) {
    static_assert(AlwaysFalse<T>(), "Backend not supported!");
  }
};

// WYHash Backend
template <>
class OHash<OHashBackend::WYHash> {
public:
  OHash() : mod_(0xffffffff), mask_(0xffffffff), seed_(3000750715), seed2_(2241463753), mode64_(false) {}
  // Valid masks must be of the form 2^n - 1 (not checked in this)
  OHash(uint32_t mask, uint32_t mod, uint32_t seed, uint32_t seed2 = 2241463753, bool mode = false) : 
      mod_(mod), mask_(mask), seed_(seed), seed2_(seed2), mode64_(mode) {}

  template <typename T>
  uint32_t operator()(const T& x) const {
    uint32_t result = wyhash32(&x, sizeof(T), seed_);
    return result & mask_;
  }

  // mask is applied only to the higher 32 bits.
  template <typename T>
  uint64_t hash64(const T& x) const {
    uint64_t result = wyhash32(&x, sizeof(T), seed_) & mask_;
    if (mode64_) {
      result <<= 32;
      result |= wyhash32(&x, sizeof(T), seed2_);
    }
    return result;
  }

  // fast alternative when mod_ is much smaller than 2^32
  template <typename T>
  uint32_t hash_in_range(const T& x) const {
    int64_t result = wyhash32(&x, sizeof(T), seed_);
    int64_t product = result * mod_;
    return product >> 32;
  }
  // TODO: If mod_ needs to be larger
  // refer to the solution here
  // https://lemire.me/blog/2019/06/06/nearly-divisionless-random-integer-generation-on-various-systems/

  uint32_t operator()(int x, int y) const {
    int xy[2];
    xy[0] = x;
    xy[1] = y;
    uint32_t result = wyhash32(xy, sizeof(int) * 2, seed_);
    return result & mask_;
  }

  bool operator==(const OHash& other) const = default;
  uint32_t& seed() { return seed_; }
  uint64_t& mod() { return mod_; }
  uint32_t& mask() { return mask_; }

private:
  uint64_t mod_;
  uint32_t mask_;
  uint32_t seed_, seed2_;
  bool mode64_;
};

// WYHash64 Backend
template <>
class OHash<OHashBackend::WYHash64> {
public:
  OHash() : mod_(0xffffffff), mask_(0xffffffff), seed_(3000750715) {
    make_secret(seed_, secret_);
  }
  // Valid masks must be of the form 2^n - 1 (not checked in this)
  OHash(uint32_t seed, uint32_t mask) : mod_(mask), mask_(mask), seed_(seed) {
    make_secret(seed_, secret_);
  }

  uint32_t operator()(int x) const {
    uint32_t result = wyhash(&x, sizeof(int), seed_, secret_);
    return result & mask_;
  }

  // Caution: division is very slow.
  // Only use in time-insensitive applications.
  uint32_t hash_in_range(int x) const {
    uint32_t result = wyhash(&x, sizeof(int), seed_, secret_);
    return result % mod_;
  }

  uint32_t operator()(int x, int y) const {
    int xy[2];
    xy[0] = x;
    xy[1] = y;
    uint32_t result = wyhash(xy, sizeof(int) * 2, seed_, secret_);
    return result & mask_;
  }

  const uint32_t& seed() const { return seed_; }
  const uint32_t& mod() const { return mod_; }
  const uint32_t& mask() const { return mask_; }
  // secret is not shown to you
  const uint64_t* secret() { return secret_; }

private:
  uint32_t mod_;
  uint32_t mask_;
  uint32_t seed_;
  uint64_t secret_[4];
};

using WYHash = OHash<OHashBackend::WYHash>;

}

