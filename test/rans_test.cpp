#include "rans_wrapper.hpp"
#include "probability.hpp"
#include <cstdint>
#include <iostream>
#include <random>

#include <oniakDataStructure/ohist.h>
#include <oniakTimer/otime.h>

using namespace std;
using namespace ONIAK;
using namespace Doro;

constexpr int message_size = 10000;

template <typename FrequencyType>
void test_routine(const std::vector<uint8_t>& message, const FrequencyType& frequencies) {
  StopWatch sw;
  RansWrapper rans_wrapper(frequencies);
  RansCode code = rans_wrapper.encode(message);
  auto decoded_message = rans_wrapper.decode(code);

  cout << "Time: " << sw.peek() << " s" << endl;
  cout << "Message size: " << message_size << " bytes" << endl;
  cout << "Compressed message size: " << code.size() << " bytes" << endl;
  cout << "Decompressed message size: " << decoded_message.size() << " bytes" << endl;
  cout << "The decoded message is the same as the original message: " << (message == decoded_message) << endl;
  for (int i = 0; i < message_size; i++) {
    if (message[i] != decoded_message[i]) {
      cout << "Mismatch at index " << i << ": " << (int)message[i] << " != " << (int)decoded_message[i] << endl;
      break;
    }
  }
  double message_entropy = entropy(frequencies);
  cout << "Message entropy: " << message_entropy * message.size() / 8 << " bytes" << endl;
}

int main() {
  //std::random_device rd;
  std::mt19937 gen(2966300);
  std::uniform_int_distribution<uint8_t> dis(0, 100);
  std::poisson_distribution<uint8_t> poisson_dis(10), poisson_dis2(5);
  Poisson poisson(10);
  Skellam skellam(5, 10);
  // std::geometric_distribution<uint8_t> poisson_dis(0.1);

  std::vector<uint8_t> message1(message_size), message2(message_size), message3(message_size);
  for (int i = 0; i < message_size; i++) {
    message1[i] = dis(gen);
    message2[i] = poisson_dis(gen);
    message3[i] = poisson_dis2(gen) - poisson_dis(gen);
  }
  message1[320] = 144;  // unexpected symbol

  auto frequencies = frequency_count(message1);
  test_routine(message1, frequencies);
  
  auto frequencies2 = poisson.pmf_map<uint8_t>();
  test_routine(message2, frequencies2);

  frequencies2 = skellam.pmf_map<uint8_t>();
  test_routine(message3, frequencies2);
  return 0;
}