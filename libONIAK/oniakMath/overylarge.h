#pragma once 

#include <bitset>
#include <format>
#include <sstream>

namespace ONIAK {
  template <size_t N>
  class VeryLargeInt : public std::bitset<N> {
    public:
    static constexpr size_t size = N;
    static constexpr VeryLargeInt mask = -1ull;
    using std::bitset<N>::bitset; // inherit constructors
    VeryLargeInt(std::bitset<N> bs) : std::bitset<N>(bs) {}
    VeryLargeInt(): std::bitset<N>() {}
    
    std::string to_hex_string() const {
      std::ostringstream oss;
      oss << std::hex << "0x";
      int rem_width = size;
      while (rem_width > 0) {
        rem_width -= 64;
        int move_bits = std::max(rem_width, 0);
        int pad_bits = std::max(-rem_width, 0);
        VeryLargeInt val = (*this >> move_bits) & (mask >> pad_bits);
        oss << val.to_ullong();
      }
      return oss.str();
    }

    std::strong_ordering operator<=>(const VeryLargeInt& other) const {
      int rem_width = size;
      while (rem_width > 0) {
        rem_width -= 64;
        int move_bits = std::max(rem_width, 0);
        int pad_bits = std::max(-rem_width, 0);
        VeryLargeInt val = (*this >> move_bits) & (mask >> pad_bits);
        VeryLargeInt other_val = (other >> move_bits) & (mask >> pad_bits);
        auto comp_result = val.to_ullong() <=> other_val.to_ullong();
        if (comp_result != std::strong_ordering::equal) {
          return comp_result;
        }
      }
      return std::strong_ordering::equal;
    }

    // only keeps the lower 64 bits
    explicit operator uint64_t() const {
      return (*this & mask).to_ullong();
    }
  };

  template <typename T>
  concept isVeryLargeInt = std::is_same_v<T, VeryLargeInt<T::size>>;

}

template <size_t N>
struct std::formatter<ONIAK::VeryLargeInt<N>> : std::formatter<std::string> {
auto format(const ONIAK::VeryLargeInt<N>& p, format_context& ctx) const {
  return formatter<string>::format(
    p.to_hex_string(), ctx);
}
};

template <size_t N>
struct std::hash<ONIAK::VeryLargeInt<N>> {
  size_t operator()(const ONIAK::VeryLargeInt<N>& vli) const {
    return std::hash<std::bitset<N>>()(vli);  // use the hash method of bitset
  }
};