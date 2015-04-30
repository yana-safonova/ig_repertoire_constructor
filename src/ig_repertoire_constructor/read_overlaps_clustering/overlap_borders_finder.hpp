#pragma once

#include <memory>
#include <vector>
#include <utility>

namespace ig_repertoire_constructor {

class OverlapBordersFinder {
public:
    virtual std::pair <unsigned, unsigned> FindOptimalOverlapBorders(std::shared_ptr <ReadArchive> read_archive_ptr,
                                                                     const AlignedReadCluster & cluster,
                                                                     const std::vector <size_t> & indices) const = 0;

    virtual ~OverlapBordersFinder() {}
};

class OverlapBordersFinderFactory {
public:
    std::shared_ptr <OverlapBordersFinder> Create();
};

class CutTailsOverlapBordersFinder : public OverlapBordersFinder {
public:
    virtual std::pair <unsigned, unsigned> FindOptimalOverlapBorders(std::shared_ptr <ReadArchive> read_archive_ptr,
                                                                     const AlignedReadCluster & cluster,
                                                                     const std::vector <size_t> & indices) const;
};

class OptimalOverlapBordersFinder : public OverlapBordersFinder {
public:
    virtual std::pair <unsigned, unsigned> FindOptimalOverlapBorders(std::shared_ptr <ReadArchive> read_archive_ptr,
                                                                     const AlignedReadCluster & cluster,
                                                                     const std::vector <size_t> & indices) const;
};


std::shared_ptr <OverlapBordersFinder> OverlapBordersFinderFactory::Create() {
    return std::make_shared <OptimalOverlapBordersFinder>();
//    return std::make_shared <CutTailsOverlapBordersFinder>();
}

inline std::pair <unsigned, unsigned> CutTailsOverlapBordersFinder::FindOptimalOverlapBorders(std::shared_ptr <ReadArchive> read_archive_ptr,
                                                                                              const AlignedReadCluster & cluster,
                                                                                              const std::vector <size_t> & indices) const {
    std::vector <unsigned> lefts;
    std::vector <unsigned> rights;
    lefts.reserve(indices.size());
    rights.reserve(indices.size());
    for (auto index : indices) {
        const auto & aligned_read = cluster[index];
        const io::SingleRead & read = read_archive_ptr->operator[](aligned_read.read_number);
        lefts.push_back(aligned_read.shift);
        rights.push_back((unsigned)read.size() + aligned_read.shift);
    }
    std::sort(lefts.begin(), lefts.end(), std::greater <unsigned>());
    std::sort(rights.begin(), rights.end());

    unsigned exact_left_overlap = lefts[0];
    unsigned exact_right_overlap = rights[0];

    const double outliers_percent = 0.01; // TODO tune this constants
    const double threshold = 0.6;

    unsigned inexact_left_overlap = lefts[(int)(outliers_percent * (double)indices.size())];
    unsigned inexact_right_overlap = rights[(int)(outliers_percent * (double)indices.size())];

    if ((1.0 + threshold) * ((int)exact_right_overlap - (int)exact_left_overlap)
            < (1.0 - outliers_percent) * ((int)inexact_right_overlap - (int)inexact_left_overlap)) {
        INFO("size = " << indices.size() << ", exact = [" << exact_left_overlap << ", " << exact_right_overlap
                        << "], inexact = [" << inexact_left_overlap << ", " << inexact_right_overlap << "]");
        return std::make_pair(inexact_left_overlap, inexact_right_overlap);
    } else {
        return std::make_pair(exact_left_overlap, exact_right_overlap);
    }
}

inline std::pair <unsigned, unsigned> OptimalOverlapBordersFinder::FindOptimalOverlapBorders(std::shared_ptr <ReadArchive> read_archive_ptr,
                                                                                             const AlignedReadCluster & cluster,
                                                                                             const std::vector <size_t> & indices) const {
    size_t n = indices.size();
    std::vector <unsigned> lefts;
    std::vector <unsigned> rights;
    std::vector <unsigned> positions;
    lefts.reserve(n);
    rights.reserve(n);
    positions.reserve(2 * n);
    for (auto index : indices) {
        const auto & aligned_read = cluster[index];
        const io::SingleRead & read = read_archive_ptr->operator[](aligned_read.read_number);
        lefts.push_back(aligned_read.shift);
        rights.push_back((unsigned)read.size() + aligned_read.shift);
        positions.push_back(lefts.back());
        positions.push_back(rights.back());
    }

    std::sort(positions.begin(), positions.end());
    positions.erase(std::unique(positions.begin(), positions.end()), positions.end());

    size_t m = positions.size();

    std::vector <std::vector <int> > count(m, vector <int>(m, 0));
    for (size_t i = 0; i < n; ++i) {
        size_t L = std::lower_bound(positions.begin(), positions.end(), lefts[i]) - positions.begin();
        size_t R = std::lower_bound(positions.begin(), positions.end(), rights[i]) - positions.begin();
        count[L][R] += 1;
    }
    std::vector <std::vector <int> > dp(m, vector <int>(m, 0));
    for (size_t len = m - 1; len > 0; --len) {
        for (size_t L = 0, R = len; R < m; ++L, ++R) {
            dp[L][R] = count[L][R];
            if (L > 0) {
                dp[L][R] += dp[L - 1][R];
            }
            if (R < m - 1) {
                dp[L][R] += dp[L][R + 1];
                if (L > 0) {
                    dp[L][R] -= dp[L - 1][R + 1];
                }
            }
        }
    }
    std::pair <unsigned, unsigned> best_borders(-1, -1);
    double best_score = -1;
    double max_overlap_length = positions.back() - positions[0];
    for (size_t L = 0; L < m; ++L) {
        std::pair <unsigned, unsigned> current_borders;
        current_borders.first = positions[L];
        for (size_t R = L + 1; R < m; ++R) {
            current_borders.second = positions[R];
            unsigned current_overlap_length = current_borders.second - current_borders.first;
            if (current_overlap_length  < 300) {
                continue;
            }
            double current_score = double(dp[L][R]) * max_overlap_length + current_overlap_length;
            if (current_score > best_score) {
                best_score = current_score;
                best_borders = current_borders;
            }
        }
    }
    return best_borders;
}


}
