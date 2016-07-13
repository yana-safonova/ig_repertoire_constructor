#pragma once

#include <seqan/sequence.h>
#include <seqan/basic.h>
#include <logger/logger.hpp>
#include <logger/log_writers.hpp>

namespace pog {

using id_t = seqan::CharString;
using nt_t = seqan::Dna5;
using seq_t = seqan::String<nt_t>;

struct seqan_read {
    seqan_read(id_t id, seq_t sequence)
            : id(id)
            , sequence(sequence) {
        // EMPTY
    }

    id_t id;
    seq_t sequence;
};

struct directed_seqan_read {
    directed_seqan_read(seqan_read* read, bool straight = true)
        : read(read)
        , straight(straight) {}

    directed_seqan_read(id_t id, seq_t sequence, bool straight = true)
        : read(new seqan_read(id, sequence))
        , straight(straight) {}

    ~directed_seqan_read() {
        delete read;
    }

    directed_seqan_read* reverse_complement() {
        return new directed_seqan_read(read, !straight);
    }

    id_t& id() {
        return read->id;
    }

    seq_t& sequence() {
        return read->sequence;
    }

    seqan_read* read;
    bool straight; // false is for reverse complement
};

} // namespace pog
