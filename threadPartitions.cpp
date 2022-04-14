#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <set>
#include <vector>

struct comp {
    template<typename T>
    bool operator()(const T &l, const T &r) const {
        if (l.first == r.first) {
            return l.second > r.second;
        }
 
        return l.first < r.first;
    }
};

class PartitionGenerator {
   public:
    struct RandContext {
        explicit RandContext(uint64_t seed) : mt(seed) {}

        std::mt19937_64 mt;
    };
    void setDocIds(uint64_t numDocs) {
        for (uint64_t i = 0; i < numDocs; ++i) {
            _docIds.insert(generateRandNumWithinRange(
                0, std::numeric_limits<uint64_t>::max() - 1));
        }
    }
    void setDefaultDocIds() {
        _docIds = {2,    3,    1450, 1550,  1650,  1670,
                   1680, 1690, 2000, 10000, 11000, 12000};
        // _docIds = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1000};
    }
    uint64_t generateRandNumWithinRange(uint64_t min, uint64_t max) {
        uint64_t num = min + (_randContext.mt() % (max + 1 - min));
        return num;
    }
    PartitionGenerator(uint64_t budget)
        : _budget(budget), _seed(std::random_device()()), _randContext{_seed} {}

    std::vector<std::pair<uint64_t, uint64_t>>
    partitionDocIdsForThreadsVeryNaive(uint64_t nThreads, uint64_t budget) {
        std::vector<std::pair<uint64_t, uint64_t>> partitionsRanges;
        uint32_t distPerThread = _docIds.size() / nThreads;
        uint64_t const lowestBound = *(_docIds.begin());
        uint64_t const highestBound = *(_docIds.rbegin());

        uint64_t const rangeSize = (lowestBound + highestBound) / nThreads;
        uint64_t ub = lowestBound - 1;
        for (uint32_t i = 0; i < nThreads; ++i) {
            uint64_t lb = ub + 1;
            ub = std::min(highestBound, lb + rangeSize);
            partitionsRanges.emplace_back(lb, ub);
        }
        return partitionsRanges;
    }

    void detectGapsRecursively(uint64_t min, uint64_t max) {
        static uint64_t currBudget = _budget;
        currBudget--;
        if (currBudget != 0 && (max - min > 1)) {
            uint64_t middleNum = (min + max) / 2;
            auto middleRangeIt = _docIds.equal_range(middleNum);
            if (*(middleRangeIt.first) == *(middleRangeIt.second)) {
                middleRangeIt.first = std::prev(middleRangeIt.first);
            }
            _gaps.emplace_back(*(middleRangeIt.first), *(middleRangeIt.second));
            detectGapsRecursively(min, *middleRangeIt.first);
            detectGapsRecursively(*middleRangeIt.second, max);
        }
    }

    void detectGapsIncrementally(uint64_t min, uint64_t max) {
        std::deque<std::pair<uint64_t, uint64_t>> ranges{{min, max}};
        uint64_t currBudget = _budget;
        while (currBudget > 0 && !ranges.empty()) {
            currBudget--;
            auto currRange = ranges.front();
            ranges.pop_front();
            uint64_t middleNum = (currRange.first + currRange.second) / 2;
            auto middleRangeIt = _docIds.equal_range(middleNum);
            if (*(middleRangeIt.first) == *(middleRangeIt.second)) {
                middleRangeIt.first = std::prev(middleRangeIt.first);
            }
            _gaps.emplace_back(*(middleRangeIt.first), *(middleRangeIt.second));
            if (*middleRangeIt.first - currRange.first > 1) {
                ranges.emplace_back(currRange.first, *middleRangeIt.first);
            }
            if (currRange.second - *middleRangeIt.second > 1) {
                ranges.emplace_back(*middleRangeIt.second, currRange.second);
            }
        }
    }

    void computePotentiallyContinuousSegments() {
        std::sort(_gaps.begin(), _gaps.end(), [](std::pair<uint64_t, uint64_t> pair1, std::pair<uint64_t, uint64_t> pair2) {
            return pair1.first < pair2.first;
        }); 
        _potentiallyContinuousSegments.emplace_back(*(_docIds.begin()), _gaps[0].first);
        for (uint64_t i = 0; i < _gaps.size() - 1; ++i) {
           _potentiallyContinuousSegments.emplace_back(_gaps[i].second, _gaps[i + 1].first); 
        }
        _potentiallyContinuousSegments.emplace_back(*(_docIds.rbegin()), _gaps.back().second);

    }

    void mapDownscale() {
       
        uint64_t lastValidNum = *(_docIds.begin());
        for (auto const& [l, u] : _potentiallyContinuousSegments) {
            uint64_t lMapped = lastValidNum;
            uint64_t uMapped = lMapped + u - l;
            lastValidNum = uMapped + 1;
            std::cout << l << " " << u << " maps to " << lMapped << " " << uMapped << '\n';
            _originalToDownscaled.insert({{l, u}, {lMapped, uMapped}});
            _downscaledToOriginal.insert({{lMapped, uMapped}, {l, u}});
        }
    }

    void divideDownscaled(uint64_t nThreads) {
         uint64_t continuousSegmentsSize = std::accumulate(_potentiallyContinuousSegments.begin(), _potentiallyContinuousSegments.end(), 0, [](auto acc, auto pair1) {
            uint64_t subtraction = pair1.second - pair1.first;
            std::cout << acc << " " << subtraction << " " << ((subtraction  > 1) ? 2 : 1) << '\n';
            return acc + subtraction + ((subtraction  > 1) ? 2 : 1);
        });
        std::cout << "SIZE " << continuousSegmentsSize << '\n';
        uint64_t const rangeSize = (continuousSegmentsSize + _potentiallyContinuousSegments[0].first) / nThreads;
        std::cout << "range size " << rangeSize << '\n';
        int64_t lastUpper = _potentiallyContinuousSegments[0].first - 1;
        for (uint64_t i = 0; i < nThreads; ++i) {
            uint64_t lower = lastUpper + 1;
            uint64_t upper = std::min(lower + rangeSize, continuousSegmentsSize);
            std::cout << "thread " << i << " range " << lower << " " << upper << '\n';
            _downscaledPartitions.emplace_back(lower, upper);
            lastUpper = lastUpper + 1 + rangeSize; 
        }
    }

    std::pair<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>> findRangeThatContains(uint64_t num) {
        auto it = std::find_if(_downscaledToOriginal.begin(), _downscaledToOriginal.end(), [&num](auto p) {
            std::cout << num << " p " << p.first.first << " " << p.first.second << '\n';
            return num >= p.first.first && num <= p.first.second;
        });
        if (it == _downscaledToOriginal.end()) {
                std::cout << "NAO DEU \n"; 
                exit(0);
        }
        return *it;
    }

    uint64_t downscaledToOriginal(uint64_t num) {
        for (auto const&[downscaled, original] : _downscaledToOriginal) {
            std::cout << "map " << downscaled.first << " " << downscaled.second << " " << original.first << " " << original.second << '\n';
        }
        auto [downscaled, original] = findRangeThatContains(num);
        return(num - downscaled.first + original.first);

    }
 
    void upscale() {
        for (auto const& [l, u] : _downscaledPartitions) {
            auto lOriginalPos = downscaledToOriginal(l);
            auto uOriginalPos = downscaledToOriginal(u);
            _partitions.emplace_back(lOriginalPos, uOriginalPos);
        }
    }

    std::vector<std::pair<uint64_t, uint64_t>> partitionDocIdsForThreads(
        uint64_t nThreads) {
        detectGapsIncrementally(*(_docIds.begin()), *(_docIds.rbegin()));
        computePotentiallyContinuousSegments();
        std::cout << "GAPS \n";
        for (const auto &[k, v] : _gaps) {
            std::cout << k << ": " << v << '\n';
        }
        std::cout << "---------------------------\n";
        std::cout << "SEGMENTS \n";
        for (const auto &[k, v] : _potentiallyContinuousSegments) {
            std::cout << k << ": " << v << '\n';
        }
        mapDownscale();
        std::cout << "MAP \n";
        for (auto const& [k, v] : _originalToDownscaled) {
            std::cout << k.first << " " << k.second << ", " << v.first << " " << v.second << '\n';
        }
        divideDownscaled(nThreads);
        upscale();
        return _partitions;
    }

   private:
    uint64_t const _budget;
    uint64_t _seed;
    RandContext _randContext;
    std::set<uint64_t> _docIds;
    std::vector<std::pair<uint64_t, uint64_t>> _gaps;
    std::vector<std::pair<uint64_t, uint64_t>> _potentiallyContinuousSegments;
    std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>, comp> _originalToDownscaled;
    std::map<std::pair<uint64_t, uint64_t>, std::pair<uint64_t, uint64_t>, comp> _downscaledToOriginal;
    std::vector<std::pair<uint64_t, uint64_t>> _downscaledPartitions;
    std::vector<std::pair<uint64_t, uint64_t>> _partitions;
};

int main() {
    PartitionGenerator partitionGen(5);
    partitionGen.setDefaultDocIds();
    std::vector<std::pair<uint64_t, uint64_t>> partitions =
        partitionGen.partitionDocIdsForThreads(3);
    std::cout << "--------------------------------- \n";
    for (const auto &[k, v] : partitions) {
        std::cout << k << ": " << v << '\n';
    }
    return 0;
}
