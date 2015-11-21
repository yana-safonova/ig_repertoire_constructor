#pragma once

template<class Recombination, class RecombinationIterator>
class AbstractRecombinationGenerator {

public:
    virtual void ComputeRecombinations() = 0;

    virtual RecombinationIterator begin() = 0;

    virtual RecombinationIterator end() = 0;

    //virtual std::size_t size() const = 0;

    virtual ~AbstractRecombinationGenerator() { }
};