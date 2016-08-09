#include "logger/logger.hpp"
#include "info_based_d_hits_calculator.hpp"

using namespace vdj_labeler;

alignment_utils::ImmuneGeneAlignmentPositions InfoBasedDHitsCalculator::CreateDAlignmentPositions(
        alignment_utils::AlignmentPositions d_alignment_positions,
        const germline_utils::ImmuneGene* gene_ptr,
        const core::Read* read_ptr) const
{
    d_alignment_positions.subject_pos.second = seqan::length(gene_ptr->seq()) - 1;
    assert(read_ptr != nullptr);
    return alignment_utils::ImmuneGeneAlignmentPositions(std::move(d_alignment_positions), gene_ptr, read_ptr);
}

std::vector<alignment_utils::ImmuneGeneReadAlignment> InfoBasedDHitsCalculator::CreateDGeneAlignments(
        const core::Read* read_ptr,
        const alignment_utils::AlignmentPositions d_positions) const
{
    std::vector<alignment_utils::ImmuneGeneReadAlignment> d_gene_hits;
    for (const auto& d_gene: d_gene_database_) {
        alignment_utils::ImmuneGeneAlignmentPositions d_alignment_pos =
            CreateDAlignmentPositions(d_positions, &d_gene, read_ptr);
        auto d_alignment = d_gene_aligner_.ComputeAlignment(d_alignment_pos);
        if (quality_checker_.AlignmentIsGood(d_alignment)) {
            d_gene_hits.emplace_back(std::move(d_alignment));
        }
    }
    return d_gene_hits;
}

std::vector<size_t> InfoBasedDHitsCalculator::CreatePreviousDGenesPositions(
        const std::vector<alignment_utils::ImmuneGeneReadAlignment>& d_gene_hits) const
{
    std::vector<size_t> prev(d_gene_hits.size());
    if (prev.size() == 0) { return prev; }

    prev[0] = INDICATE_START_ANSWER;
    for (size_t i = 1; i < d_gene_hits.size(); ++i) {
        size_t &j = prev[i];
        j = i - 1;

        while (j != INDICATE_START_ANSWER and
               d_gene_hits[j].EndQueryPosition() >= d_gene_hits[i].StartQueryPosition())
        { --j; }

        if (j == INDICATE_START_ANSWER) { continue; }

        size_t k = j;
        double current_max_score = d_gene_hits[j].Score();
        size_t gap = d_gene_hits[i].StartQueryPosition() - d_gene_hits[j].EndQueryPosition();
        while (k != INDICATE_START_ANSWER and
               d_gene_hits[j].EndQueryPosition() - d_gene_hits[k].EndQueryPosition() == gap)
        {
            double k_score = d_gene_hits[k].Score();
            if (k_score > current_max_score) {
                j = k;
                current_max_score = k_score;
            }
            --k;
        }
    }
    return prev;
}

std::vector<double> InfoBasedDHitsCalculator::CalcOptimalScore(
    const std::vector<alignment_utils::ImmuneGeneReadAlignment>& d_gene_hits,
    const std::vector<size_t>& prev_d_gene_pos) const
{
    std::vector<double> opt_score(d_gene_hits.size());
    if (opt_score.size() == 0) {
        return opt_score;
    }
    opt_score[0] = d_gene_hits[0].Score();
    for (size_t i = 1; i < d_gene_hits.size(); ++i) {
        double score_w_curr_gene = d_gene_hits[i].Score();
        if (prev_d_gene_pos[i] != INDICATE_START_ANSWER) {
            score_w_curr_gene += opt_score[prev_d_gene_pos[i]];
            if (d_gene_hits[i].AlignmentLength() <= 6) {
                score_w_curr_gene -= 3;
            }
        }
        opt_score[i] = std::max<double>(opt_score[i - 1], score_w_curr_gene);
    }
    return opt_score;
}

DGeneHits InfoBasedDHitsCalculator::CalcAnswer(const std::vector<alignment_utils::ImmuneGeneReadAlignment>& d_gene_hits,
                                               const std::vector<size_t>& prev_d_gene_pos,
                                               const std::vector<double>& opt_score,
                                               const core::Read* read_ptr) const
{
    DGeneHits d_hits(read_ptr);
    if (opt_score.size() == 0) { return d_hits; }

    std::vector<alignment_utils::ImmuneGeneReadAlignment> answer;
    size_t i = opt_score.size() - 1;
    while(i != INDICATE_START_ANSWER) {
        if (opt_score[i] == opt_score[i - 1]) {
            --i;
        } else {
            answer.push_back(d_gene_hits[i]);
            i = prev_d_gene_pos[i];
        }
    }
    std::reverse(answer.begin(), answer.end());
    d_hits.AddHit(DGeneHit(read_ptr, answer));
    return d_hits;
}

DGeneHits InfoBasedDHitsCalculator::ComputeDHits(const core::Read* read_ptr,
                                                 const std::vector<vj_finder::VGeneHit> &v_hits,
                                                 const std::vector<vj_finder::JGeneHit> &j_hits) const
{
    auto d_positions = d_alignment_positions_calculator_.ComputeDAlignmentPositions(v_hits, j_hits);

    if (!d_alignment_position_checker_.DAlignmentPositionsAreGood(d_positions)) {
        TRACE("D positions are too short to generate alignment");
        TRACE(d_positions);
        // add single empty alignment and return hits storage
        return DGeneHits(read_ptr);
    }

    auto d_gene_hits = CreateDGeneAlignments(read_ptr, d_positions);
    // for (const auto & d_gene_hit : d_gene_hits) {
    //     INFO(d_gene_hit.Subject() << " " << d_gene_hit.Score() << " " << d_gene_hit.Alignment());
    // }

    // Sort by right position on read of each gene (in alignment)
    std::sort(d_gene_hits.begin(), d_gene_hits.end(),
              [](alignment_utils::ImmuneGeneReadAlignment a, alignment_utils::ImmuneGeneReadAlignment b) {
                return a.EndQueryPosition() < b.EndQueryPosition();
              });

    // int i = 0;
    // for (const auto &d_gene_hit : d_gene_hits) {
    //     INFO(i++ << " " << d_gene_hit.Subject().name() << " " << d_gene_hit.StartQueryPosition() <<
    //                 " " << d_gene_hit.EndQueryPosition() << " " << d_gene_hit.Score());
    // }
    // Calculate "previous" d genes positions O(n^2) where n = d_gene_hits.size(). Bin Search is possible: O(n log n).
    auto prev_d_gene_pos = CreatePreviousDGenesPositions(d_gene_hits);
    // i = 0;
    // for (const auto &p: prev_d_gene_pos)
    //     INFO(i++ << " " << p);
    auto opt_score = CalcOptimalScore(d_gene_hits, prev_d_gene_pos);
    // i = 0;
    // for (const auto &p: opt_score)
    //     INFO(i++ << " " << p);
    auto answer = CalcAnswer(d_gene_hits, prev_d_gene_pos, opt_score, read_ptr);
    // i = 0;
    // for (const auto& ans: answer) {
    //     INFO("");
    //     INFO(read_ptr->seq);
    //     for (const auto &p: ans)
    //         INFO(i++ << " " << p.StartQueryPosition() << " " << p.EndQueryPosition() << " " << p.Subject().name());
    // }

    return answer;
}
