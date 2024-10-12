#include <iostream>
#include <random>
#include <vector>
#include "bch_wrapper.hpp"

using namespace std;
using namespace Doro;

int main() {
  int capability = 5, repeat = 10;
  BCHWrapper bch(10, capability);
  int message_size = bch.message_size();
  cout << "Message size: " << message_size << endl;

  uniform_int_distribution<uint8_t> dist(0, 1);
  uniform_int_distribution<size_t> corrupt_dist(0, message_size - 1);
  mt19937 mt(4239384);
  vector<uint8_t> data(message_size);
  for (int i = 0; i < message_size; ++i) {
    data[i] = dist(mt);
  }

  vector<uint8_t> parity_bits = bch.encode(data, true);
  for (int rpt = 0; rpt < repeat; ++rpt) {
    vector<uint8_t> corrupted_data(data);
    for (int i = 0; i < capability; ++i) {
      int idx = corrupt_dist(mt);
      corrupted_data[idx] = !corrupted_data[idx];
    }

    vector<uint8_t> decoded_data = bch.decode(corrupted_data, parity_bits, true);
    cout << "Correction is successful: " << ((data == decoded_data) ? "True" : "False" )<< endl;
  }
  return 0;
}