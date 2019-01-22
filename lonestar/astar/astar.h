/** Single source shortest paths -*- C++ -*-
 * @file
 * @section License
 *
 * Galois, a framework to exploit amorphous data-parallelism in irregular
 * programs.
 *
 * Copyright (C) 2011, The University of Texas at Austin. All rights reserved.
 * UNIVERSITY EXPRESSLY DISCLAIMS ANY AND ALL WARRANTIES CONCERNING THIS
 * SOFTWARE AND DOCUMENTATION, INCLUDING ANY WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR ANY PARTICULAR PURPOSE, NON-INFRINGEMENT AND WARRANTIES OF
 * PERFORMANCE, AND ANY WARRANTY THAT MIGHT OTHERWISE ARISE FROM COURSE OF
 * DEALING OR USAGE OF TRADE.  NO WARRANTY IS EITHER EXPRESS OR IMPLIED WITH
 * RESPECT TO THE USE OF THE SOFTWARE OR DOCUMENTATION. Under no circumstances
 * shall University be liable for incidental, special, indirect, direct or
 * consequential damages or loss of profits, interruption of business, or
 * related expenses which may arise from use of Software or Documentation,
 * including but not limited to those resulting from defects in Software and/or
 * Documentation, or loss or inaccuracy of data of any kind.
 *
 * @section Description
 *
 * Single source shortest paths.
 *
 * @author Andrew Lenharth <andrewl@lenharth.org>
 */
#ifndef SSSP_H
#define SSSP_H

#include <limits>
#include <string>
#include <sstream>
#include <stdint.h>

static const unsigned int DIST_INFINITY =
  std::numeric_limits<unsigned int>::max()/2 - 1;

template<typename GrNode>
struct UpdateRequestCommon {
  GrNode n;
  uint64_t gScore;
  uint64_t fScore;
  uint64_t parent;

  UpdateRequestCommon(const GrNode& N, unsigned int gScore)
    :n(N), gScore(gScore)
  {}

  UpdateRequestCommon(const GrNode& N, uint64_t gScore, uint64_t fScore, uint64_t parent)
    :n(N), gScore(gScore), fScore(fScore), parent(parent)
  {}

  UpdateRequestCommon()
    :n(), gScore(0)
  {}

  bool operator>(const UpdateRequestCommon& rhs) const {
    if (gScore > rhs.gScore) return true;
    if (gScore < rhs.gScore) return false;
    return n > rhs.n;
  }

  bool operator<(const UpdateRequestCommon& rhs) const {
    if (gScore < rhs.gScore) return true;
    if (gScore > rhs.gScore) return false;
    return n < rhs.n;
  }

  bool operator!=(const UpdateRequestCommon& other) const {
    if (gScore != other.gScore) return true;
    return n != other.n;
  }

  uintptr_t getID() const {
    return reinterpret_cast<uintptr_t>(n);
  }
};

struct SNode {
  unsigned int id;
  unsigned int dist;
  double lat;
  double lon;
  uint64_t parent;

  SNode(int _id = -1) : id(_id), dist(DIST_INFINITY), lat(0), lon(0), parent(0) {}
  std::string toString() {
    std::ostringstream s;
    s << '[' << id << "] dist: " << dist;
    return s.str();
  }
};

template <typename T, typename BucketFunc, size_t MAX_BUCKETS = 543210ul>
class SerialBucketWL {

  using Bucket      = std::deque<T>;
  using BucketsCont = std::vector<Bucket>;

  size_t m_minBucket;
  BucketFunc m_func;
  BucketsCont m_buckets;
  Bucket m_lastBucket;

  static_assert(MAX_BUCKETS > 0, "MAX_BUCKETS must be > 0");

public:
  explicit SerialBucketWL(const BucketFunc& f) : m_minBucket(0ul), m_func(f) {
    // reserve enough so that resize never reallocates memory
    // otherwise, minBucket may return an invalid reference
    m_buckets.reserve(MAX_BUCKETS);
  }

  void push(const T& item) {
    size_t b = m_func(item);
    assert(b >= m_minBucket && "can't push below m_minBucket");

    if (b < m_buckets.size()) {
      m_buckets[b].push_back(item);
      return;
    } else {
      if (b >= MAX_BUCKETS) {
        std::cerr << "Increase MAX_BUCKETS limit" << std::endl;
        m_lastBucket.push_back(item);
      } else {
        m_buckets.resize(b + 1);
        m_buckets[b].push_back(item);
      }
    }
  }

  void goToNextBucket(void) {
    while (m_minBucket < m_buckets.size() && m_buckets[m_minBucket].empty()) {
      ++m_minBucket;
    }
  }

  Bucket& minBucket(void) {
    if (m_minBucket < m_buckets.size()) {
      return m_buckets[m_minBucket];
    } else {
      return m_lastBucket;
    }
  }

  bool empty(void) const { return emptyImpl(m_minBucket); }

  bool allEmpty(void) const { return emptyImpl(0ul); }

private:
  bool emptyImpl(size_t start) const {
    for (size_t i = start; i < m_buckets.size(); ++i) {
      if (!m_buckets[i].empty()) {
        return false;
      }
    }

    return m_lastBucket.empty();
  }
};

#endif
