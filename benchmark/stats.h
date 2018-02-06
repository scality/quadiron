/* -*- mode: c++ -*- */
#ifndef __NTL_NETL_BENCH_STATS_H__
#define __NTL_NETL_BENCH_STATS_H__

#include <cmath>
#include <iostream>
#include <string>

class Stats_t {
public:
  Stats_t(const std::string& str, size_t work_load) {
    this->name = str;
    this->work_load = work_load;

    this->nb = 0;
    this->sum = 0;
    this->sum_2 = 0;
  }

  void begin() {
    nb = 0;
    sum = 0;
    sum_2 = 0;
  }

  void add(uint64_t val) {
    nb++;
    sum += val;
    sum_2 += (val * val);
  }

  void end() {
    avg = (double)sum / (double)nb;
    std_dev = sqrt((double)sum_2 / (double)nb - avg * avg);
  }

  void show() {
    std::cout << name << ":\tLatency(us) " << avg << " +/- " << std_dev;
    std::cout << "\t\tThroughput " << work_load/avg << " (MB/s)" << std::endl;
  }

private:
  uint64_t nb;
  uint64_t sum;
  uint64_t sum_2;
  double avg;
  double std_dev;
  size_t work_load;
  std::string name;
};

#endif
