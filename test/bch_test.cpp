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
    cout << "Correction is successful: " << ((data == decoded_data) ? "True" : "False") << endl;
  }

  // test encoding using multiple blocks.
  message_size = message_size * 10 + 250;
  data.resize(message_size);
  for (int i = 0; i < message_size; ++i) {
    data[i] = dist(mt);
  }
  vector<uint8_t> corrupted_data(data);
  vector<size_t> cidx_pos;
  for (int j = 0; j < 11; ++j) {
    for (int i = 0; i < capability; ++i) {
      int idx = corrupt_dist(mt);
      int cidx = j * bch.message_size() + idx;
      if (cidx < message_size){
        corrupted_data[cidx] = !corrupted_data[cidx];
        cidx_pos.push_back(cidx);
      }
    }
  }
  parity_bits = bch.encode(data, true);
  vector<uint8_t> decoded_data = bch.decode(corrupted_data, parity_bits, true);
  auto decoded_pos = bch.get_error_positions(corrupted_data, parity_bits, true);
  cout << "Encoded data size: " << bch.encoded_parity_size(data) << endl;
  cout << "Correction is successful: " << ((data == decoded_data) ? "True" : "False") << endl;
  std::sort(cidx_pos.begin(), cidx_pos.end());
  std::sort(decoded_pos.begin(), decoded_pos.end());
  cout << "Error positions are the same: " << ((cidx_pos == decoded_pos) ? "True" : "False") << endl;
  for(size_t i = 0; i < cidx_pos.size(); ++i) {
    if (cidx_pos[i] != decoded_pos[i]) {
      cout << "Error at " << i << "\t" << cidx_pos[i] << "\t" << decoded_pos[i] << endl;
      break;
    }
  }

  return 0;
}