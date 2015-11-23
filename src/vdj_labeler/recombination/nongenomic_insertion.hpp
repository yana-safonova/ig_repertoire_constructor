#pragma once

class NongenomicInsertion {
    size_t start_position_;
    size_t end_position_;

public:
    NongenomicInsertion(size_t start_position, size_t end_position) :
            start_position_(start_position),
            end_position_(end_position) { }

    size_t StartPosition() const { return start_position_; }

    size_t EndPosition() const { return end_position_; }

    size_t length() const { return end_position_ - start_position_ + 1; }
};