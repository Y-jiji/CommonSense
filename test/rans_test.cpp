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

int main() {
  int message_size = 10000;

  //std::random_device rd;
  std::mt19937 gen(2966300);
  std::uniform_int_distribution<uint8_t> dis(0, 100);
  // std::poisson_distribution<uint8_t> poisson_dis(10);
  std::geometric_distribution<uint8_t> poisson_dis(0.1);

  std::vector<uint8_t> message1(message_size), message2(message_size);
  for (int i = 0; i < message_size; i++) {
    message1[i] = dis(gen);
    message2[i] = poisson_dis(gen);
  }

  StopWatch sw;
  auto frequencies = frequency_count(message1);
  RansWrapper<uint8_t> rans_wrapper(frequencies);
  RansCode code = rans_wrapper.encode(message1);
  auto decoded_message = rans_wrapper.decode(code);

  cout << "Time: " << sw.peek() << " s" << endl;
  cout << "Message size: " << message_size << " bytes" << endl;
  cout << "Compressed message size: " << code.size() << " bytes" << endl;
  cout << "Decompressed message size: " << decoded_message.size() << " bytes" << endl;
  cout << "The decoded message is the same as the original message: " << (message1 == decoded_message) << endl;
  for (int i = 0; i < message_size; i++) {
    if (message1[i] != decoded_message[i]) {
      cout << "Mismatch at index " << i << ": " << (int)message1[i] << " != " << (int)decoded_message[i] << endl;
      break;
    }
  }
  double message_entropy = entropy(frequencies);
  cout << "Message entropy: " << message_entropy / 8 << " bytes" << endl;

  sw.reset_and_start();
  frequencies = frequency_count(message2);
  RansWrapper<uint8_t> rans_wrapper2(frequencies);
  code = rans_wrapper2.encode(message2);
  decoded_message = rans_wrapper2.decode(code);

  cout << "Time: " << sw.peek() << " s" << endl;
  cout << "Message size: " << message_size << " bytes" << endl;
  cout << "Compressed message size: " << code.size() << " bytes" << endl;
  cout << "The decoded message is the same as the original message: " << (message2 == decoded_message) << endl;
  cout << "Decompressed message size: " << decoded_message.size() << " bytes" << endl;
  for (int i = 0; i < message_size; i++) {
    if (message2[i] != decoded_message[i]) {
      cout << "Mismatch at index " << i << ": " << (int)message2[i] << " != " << (int)decoded_message[i] << endl;
      break;
    }
  }
  message_entropy = entropy(frequencies);
  cout << "Message entropy: " << message_entropy / 8 << " bytes" << endl;
  return 0;
}