#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <vector>

namespace {
struct comp {
  template <typename T> bool operator()(const T &l, const T &r) const {
    if (l.first == r.first) {
      return l.second > r.second;
    }

    return l.first < r.first;
  }
};

struct RandContext {
  explicit RandContext() : mt(dev()) {}
  uint64_t generateRandNumWithinRange(uint64_t min, uint64_t max) {
    /*
    std::chi_squared_distribution<double> dist(10);

    uint64_t num;
    do {
      num = (uint64_t)(max*dist(mt)/10) + min;
    }
    while (num > max);
    return num;
    */

    std::uniform_int_distribution<std::mt19937_64::result_type> dist(min, max);
    // std::fisher_f_distribution<std::mt19937_64::result_type> dist(10, 10);
    uint64_t num = dist(mt);
    return num;
  }

  std::set<uint64_t> generateDocumentIds(uint64_t numItems) {
    std::set<uint64_t> randSet;

    for (uint64_t i = 0; i < numItems; ++i) {
      randSet.insert(std::pow(i, 5));
    }

    /*
    std::string str;
    std::ifstream file("/home/cpjulia/Downloads/1000000.txt");
    while (std::getline(file, str, ',')) {
      randSet.insert(stoull(str));
    }
     */
    return randSet;
  }

  std::set<uint64_t> generateRandSet(uint64_t numItems) {
    std::set<uint64_t> randSet;
    while (randSet.size() != numItems) {
      randSet.insert(generateRandNumWithinRange(
          //      0, 100));
          0, (std::numeric_limits<uint64_t>::max() - 1)));
      if (randSet.size() % 100000 == 0) {
        std::cout << randSet.size() << " \n";
      }
    }
    return randSet;
  }
  std::random_device dev;
  std::mt19937_64 mt;
};

std::set<uint64_t> setDefaultSet() {

  return {11486, 23846, 58131, 54165, 58072, 80333, 95314, 18655, 84940, 84371};
}

} // namespace

class PartitionGenerator {
public:
  using interval = std::pair<uint64_t, uint64_t>;

  PartitionGenerator(uint64_t budget, std::set<uint64_t> docIds)
      : _budget{budget}, _docIds{std::move(docIds)} {}

  std::set<uint64_t> const &getDocIds() { return _docIds; }

  double stddev(std::vector<size_t> const &v) {
    assert(!v.empty());
    size_t sum = std::accumulate(v.begin(), v.end(), 0);
    double mean = static_cast<double>(sum) / v.size();

    std::vector<double> diff(v.size());
    std::transform(v.begin(), v.end(), diff.begin(),
                   std::bind2nd(std::minus<double>(), mean));
    double sqSum =
        std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    return std::sqrt(sqSum / v.size());
  }

  std::vector<size_t>
  getCounts(std::vector<std::pair<uint64_t, uint64_t>> const &partitions) {
    std::vector<size_t> counts;
    counts.reserve(partitions.size());

    for (auto const &p : partitions) {
      auto &counter = counts.emplace_back(0);

      auto it = _docIds.lower_bound(p.first);
      //      std::cout << "p.first " << p.first << " p.second " << p.second <<
      //      " _docIds.lower_bound " << *it << '\n';
      while (it != _docIds.end() && (*it) <= p.second) {
        //        std::cout << "value is within boundaries " << *it << '\n';
        ++counter;
        ++it;
      }
    }

    return counts;
  }

  void
  getStatistics(std::vector<std::pair<uint64_t, uint64_t>> const &partitions) {
    auto counts = getCounts(partitions);

    std::cout << "result (" << partitions.size()
              << " partition(s)):" << std::endl;
    size_t j = 0;
    uint64_t docsInPartitions = 0;
    for (auto const &p : partitions) {
      docsInPartitions += counts[j];
      std::cout << " - [" << p.first << " - " << p.second
                << "): size: " << counts[j++] << std::endl;
    }
    std::cout << "size " << _docIds.size() << '\n';
    std::cout << "diff " << _docIds.size() - docsInPartitions << '\n';

    std::cout << "stddev: " << stddev(counts) << std::endl;
  }

  void detectGapsIncrementally(uint64_t min, uint64_t max) {
    std::deque<std::pair<uint64_t, uint64_t>> ranges{{min, max}};
    uint64_t currBudget = _budget;
    while (currBudget > 0 && !ranges.empty()) {
      currBudget--;
      auto currRange = ranges.front();
      ranges.pop_front();
      uint64_t middleNum = (currRange.first + currRange.second) / 2;
      //  std::cout << "middle number " << middleNum << '\n';
      auto middleRangeIt = _docIds.equal_range(middleNum);
      //   std::cout << *(middleRangeIt.first) << " " << *(middleRangeIt.second)
      //   << '\n';
      if (*(middleRangeIt.first) == *(middleRangeIt.second)) {
        middleRangeIt.first = std::prev(middleRangeIt.first);
      }
      //   std::cout << *(middleRangeIt.first) << " " << *(middleRangeIt.second)
      //   << '\n';
      _gaps.emplace_back(*(middleRangeIt.first), *(middleRangeIt.second));
      if (*middleRangeIt.first - currRange.first > 1) {
        ranges.emplace_back(currRange.first, *middleRangeIt.first);
      }
      if (currRange.second - *middleRangeIt.second > 1) {
        ranges.emplace_back(*middleRangeIt.second, currRange.second);
      }
      //   std::cout << "current range " << ranges.back().first << " "
      //             << ranges.back().second << '\n';
    }
  }

  void computePotentiallyContinuousSegments() {
    std::sort(_gaps.begin(), _gaps.end(),
              [](std::pair<uint64_t, uint64_t> pair1,
                 std::pair<uint64_t, uint64_t> pair2) {
                return pair1.first < pair2.first;
              });
    _potentiallyContinuousSegments.emplace_back(*(_docIds.begin()),
                                                _gaps[0].first);
    for (uint64_t i = 0; i < _gaps.size() - 1; ++i) {
      _potentiallyContinuousSegments.emplace_back(_gaps[i].second,
                                                  _gaps[i + 1].first);
    }
    //   std::cout << "_docIds begin " << *(_docIds.begin()) << " end " <<
    //   *(_docIds.rbegin()) << " " << _gaps.back().second << '\n';
    _potentiallyContinuousSegments.emplace_back(_gaps.back().second,
                                                *(_docIds.rbegin()));
  }

  void mapDownscale() {
    uint64_t lastValidNum = *(_docIds.begin());
    for (auto const &[l, u] : _potentiallyContinuousSegments) {
      uint64_t lMapped = lastValidNum;
      uint64_t uMapped = lMapped + u - l;
      lastValidNum = uMapped + 1;
      //    std::cout << l << " " << u << " maps to " << lMapped << " " <<
      //    uMapped
      //              << '\n';
      _originalToDownscaled.insert({{l, u}, {lMapped, uMapped}});
      _downscaledToOriginal.insert({{lMapped, uMapped}, {l, u}});
    }
  }

  void divideDownscaled(size_t nThreads) {
    uint64_t continuousSegmentsSize =
        std::accumulate(_potentiallyContinuousSegments.begin(),
                        _potentiallyContinuousSegments.end(), uint64_t{},
                        [](uint64_t acc, auto pair1) {
                          uint64_t subtraction = pair1.second - pair1.first;
                          //     std::cout << acc << " " << subtraction + 1 <<
                          //     '\n';
                          return acc + subtraction + 1;
                        });
    //   std::cout << "SIZE " << continuousSegmentsSize << '\n';
    uint64_t const rangeSize =
        //   (continuousSegmentsSize + _potentiallyContinuousSegments[0].first)
        //   /
        continuousSegmentsSize / nThreads;
    //  std::cout << "range size " << rangeSize << '\n';
    uint64_t lastUpper = _potentiallyContinuousSegments[0].first;
    for (uint64_t i = 0; i < nThreads; ++i) {
      uint64_t lower = lastUpper;
      uint64_t upper = lower + rangeSize - 1;
      //   std:: cout << "lower " << lower << '\n' << " upper " << upper <<
      //   '\n';

      //      std::cout << "thread " << i << " range " << lower << " " << upper
      //      <<
      //      '\n';
      _downscaledPartitions.emplace_back(lower, upper);
      lastUpper = lastUpper + rangeSize;
      //        std::cout << "lastUpper " << lastUpper << '\n';
    }
    _downscaledPartitions.back().second =
        _originalToDownscaled.rbegin()->second.second;
  }

  std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>>
  findRangeThatContains(uint64_t num) {
    auto it =
        std::find_if(_downscaledToOriginal.begin(), _downscaledToOriginal.end(),
                     [&num](auto p) {
                       //             std::cout << num << " p " << p.first.first
                       //             << " "
                       //                       << p.first.second << '\n';
                       return num >= p.first.first && num <= p.first.second;
                     });
    if (it == _downscaledToOriginal.end()) {
      std::cout << "NAO DEU \n";
      exit(0);
    }
    return *it;
  }

  uint64_t downscaledToOriginal(uint64_t num) {
    //   std::cout << "num " << num << '\n';
    auto [downscaled, original] = findRangeThatContains(num);
    return (num - downscaled.first + original.first);
  }

  void upscale() {
    for (auto const &[l, u] : _downscaledPartitions) {
      //   std::cout << "l " << l << ", u " << u << '\n';
      auto lOriginalPos = downscaledToOriginal(l);
      auto uOriginalPos = downscaledToOriginal(u);

      _partitions.emplace_back(lOriginalPos, uOriginalPos);
    }
  }

  std::vector<std::pair<uint64_t, uint64_t>>
  partitionDocIdsForThreads(size_t nThreads) {
    std::cout << "Generating partitions for " << _docIds.size()
              << " documents, " << nThreads << " threads, "
              << "budget = " << _budget << '\n';
    assert(_budget > 0);
    assert(nThreads >= 1);
    assert(nThreads <= 16);

    detectGapsIncrementally(*(_docIds.begin()), *(_docIds.rbegin()));
    computePotentiallyContinuousSegments();

    /*
    std::cout << "GAPS \n";
    for (const auto &[k, v] : _gaps) {
      std::cout << k << ": " << v << '\n';
    }
    std::cout << "---------------------------\n";
    std::cout << "SEGMENTS \n";
    for (const auto &[k, v] : _potentiallyContinuousSegments) {
      std::cout << k << ": " << v << '\n';
    }
     */

    mapDownscale();

    /*
    std::cout << "MAP \n";
    for (auto const &[k, v] : _originalToDownscaled) {
      std::cout << k.first << " " << k.second << ", " << v.first << " "
                << v.second << '\n';
    }
    */

    divideDownscaled(nThreads);
    upscale();
    return _partitions;
  }

  std::vector<std::pair<uint64_t, uint64_t>>
  partitionsDocIdsNaive(size_t nThreads) {
    assert(nThreads >= 1);
    assert(nThreads <= 16);

    std::vector<std::pair<uint64_t, uint64_t>> result;

    if (_docIds.empty()) {
      result.push_back({0, 1});
    } else {
      uint64_t first = *_docIds.begin();
      uint64_t last = *_docIds.rbegin();
      double step = static_cast<double>(last - first + 1) / nThreads;
      std::cout << "step " << step << " " << static_cast<uint64_t>(step);
      for (size_t i = 0; i < nThreads; ++i) {
        last = std::max<uint64_t>(first + 1, std::ceil(i + 1) * step);
        if ((i == nThreads - 1) && (last != *_docIds.rbegin())) {
          last = *_docIds.rbegin();
        }
        result.push_back({first, last});
        first = last;
      }
    }

    return result;
  }

private:
  uint64_t const _budget;
  std::set<uint64_t> _docIds;
  std::vector<std::pair<uint64_t, uint64_t>> _gaps;
  std::vector<std::pair<uint64_t, uint64_t>> _potentiallyContinuousSegments;
  std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>, comp>
      _originalToDownscaled;
  std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>, comp>
      _downscaledToOriginal;
  std::vector<std::pair<uint64_t, uint64_t>> _downscaledPartitions;
  std::vector<std::pair<uint64_t, uint64_t>> _partitions;
};

int main() {
  RandContext randContext;
  uint64_t numDocs = 300;
  //std::set<uint64_t> docs = randContext.generateRandSet(numDocs);
  std::set<uint64_t> docs = randContext.generateDocumentIds(numDocs);
  for (uint8_t i = 1; i < 12; ++i) {
    //   PartitionGenerator partitionGen(5*i, ::setDefaultSet());
    PartitionGenerator partitionGen(std::pow(2, i) - 1, docs);
    std::vector<std::pair<uint64_t, uint64_t>> partitions =
        //    partitionGen.partitionDocIdsForThreads(randContext.generateRandNumWithinRange(3));
        partitionGen.partitionDocIdsForThreads(5);
    std::cout << "Results for approach 1:" << std::endl;
    partitionGen.getStatistics(partitions);
  }
  PartitionGenerator partitionGen(1, docs);
  std::vector<std::pair<uint64_t, uint64_t>> partitionsNaive =
      partitionGen.partitionsDocIdsNaive(5);
  std::cout << "Results for approach 2:" << std::endl;
  partitionGen.getStatistics(partitionsNaive);
  return 0;
}
