//
// Created by Andrew Bzikadze on 11/9/16.
//

#pragma once

#include <cmath>
#include <vector>

#include "seqan/basic.h"
#include "verify.hpp"

#include "kmer_utils.hpp"

namespace shm_kmer_matrix_estimator {

template <typename T>
class KmerIndexedVector {
private:
    std::vector<T> inner_data_;

public:
    KmerIndexedVector() = default;

    KmerIndexedVector(const KmerIndexedVector &) = default;
    KmerIndexedVector(KmerIndexedVector &&)      = default;
    KmerIndexedVector &operator=(const KmerIndexedVector &) = default;
    KmerIndexedVector &operator=(KmerIndexedVector &&)      = default;

    explicit KmerIndexedVector(size_t n, const T& val = T()) :
        inner_data_(n, val)
    { }

    virtual ~KmerIndexedVector() = default;

    const std::vector<T>& inner_data() { return inner_data_; }

    typedef typename std::vector<T>::iterator iterator;
    typedef typename std::vector<T>::const_iterator const_iterator;
    typedef typename std::vector<T>::reverse_iterator reverse_iterator;
    typedef typename std::vector<T>::const_reverse_iterator const_reverse_iterator;

    iterator       begin ()       { return inner_data_.begin (); }
    const_iterator begin () const { return inner_data_.begin (); }
    const_iterator cbegin() const { return inner_data_.cbegin(); }
    iterator       end   ()       { return inner_data_.end   (); }
    const_iterator end   () const { return inner_data_.end   (); }
    const_iterator cend  () const { return inner_data_.cend  (); }

    reverse_iterator       rbegin ()       { return inner_data_.rbegin (); }
    const_reverse_iterator rbegin () const { return inner_data_.rbegin (); }
    const_reverse_iterator crbegin() const { return inner_data_.crbegin(); }
    reverse_iterator       rend   ()       { return inner_data_.rend   (); }
    const_reverse_iterator rend   () const { return inner_data_.rend   (); }
    const_reverse_iterator crend  () const { return inner_data_.crend  (); }

    size_t size() const noexcept { return inner_data_.size(); }
    size_t max_size() const noexcept { return inner_data_.max_size(); }
    void resize (size_t n) { inner_data_.resize(n); }
    void resize (size_t n, const T& val) { inner_data_.resize(n, val); }
    size_t capacity() const noexcept { return inner_data_.size(); }
    bool empty() const noexcept { return inner_data_.size() == 0; }
    void reserve (size_t n) { inner_data_.reserve(n); }
    void shrink_to_fit() { inner_data_.shrink_to_fit(); }

    T& operator[] (size_t n) { return inner_data_[n]; }
    const T& operator[] (size_t n) const { return inner_data_[n]; }

    // Warning! No checks are made!
    T& operator[] (const std::string& kmer) {
        return (*this)[KmerUtils::GetIndexByKmer(kmer)];
    }

    const T& operator[] (const std::string& kmer) const {
        return (*this)[KmerUtils::GetIndexByKmer(kmer)];
    }

    T& at(size_t n) { return inner_data_.at(n); }
    const T& at(size_t n) const { return inner_data_.at(n); }

    T& at(const std::string& kmer) {
        return this -> at(KmerUtils::GetIndexByKmer(kmer));
    }

    const T& at(const std::string& kmer) const {
        return this -> at(KmerUtils::GetIndexByKmer(kmer));
    }

    T& front() { return inner_data_.front(); }
    const T& front() const { return inner_data_.front(); }
    T& back() { return inner_data_.back(); }
    const T& back() const { return inner_data_.back(); }

    T* data() noexcept { return inner_data_.data(); }
    const T* data() const noexcept { return inner_data_.data(); }

    template <class InputIterator>
    void assign (InputIterator first, InputIterator last) { inner_data_.assign(first, last); }
    void assign (size_t n, const T& val) { inner_data_.assign(n, val); }
    void assign (std::initializer_list<T> il) { inner_data_.assign(il); }

    void push_back (const T& val) { inner_data_.push_back(val); }
    void push_back (T&& val) { inner_data_.push_back(val); }

    void pop_back() { inner_data_.pop_back(); }

    iterator insert (const_iterator position, const T& val) { return inner_data_.insert(position, val); }
    iterator insert (const_iterator position, size_t n, const T& val) { return inner_data_.insert(position, n, val); }

    template <class InputIterator>
    iterator insert (const_iterator position, InputIterator first, InputIterator last) {
        return inner_data_.insert(position, first, last);
    }

    iterator insert (const_iterator position, T&& val) { return inner_data_.insert(position, val); }
    iterator insert (const_iterator position, std::initializer_list<T> il) { return inner_data_.insert(position, il); }

    iterator erase (const_iterator position) { return inner_data_.erase(position); }
    iterator erase (const_iterator first, const_iterator last) { return inner_data_.erase(first, last); }

    void clear() noexcept { return inner_data_.clear(); }

    template <class... Args>
    iterator emplace (const_iterator position, Args&&... args) {
        return inner_data_.emplace(position, std::forward<Args>(args)...);
    }

    template <class... Args>
    void emplace_back (Args&&... args) { return inner_data_.emplace_back(std::forward<Args>(args)...); }


    /*********************************************************************/
    size_t kmer_len() const {
        double kmer_len_ = log(static_cast<double>(size())) /
                           log(static_cast<double>(seqan::ValueSize<seqan::Dna>::VALUE));
        VERIFY_MSG(floor(kmer_len_) == kmer_len_, "kmer_len_ is not integer");
        return static_cast<size_t>(floor(kmer_len_));
    }
};

template <class T, class Alloc>
bool operator== (const std::vector<T,Alloc>& lhs, const std::vector<T,Alloc>& rhs)
{ return lhs.inner_data() == rhs.inner_data(); }

template <class T, class Alloc>
bool operator!= (const std::vector<T,Alloc>& lhs, const std::vector<T,Alloc>& rhs)
{ return lhs.inner_data() != rhs.inner_data(); }

template <class T, class Alloc>
bool operator<  (const std::vector<T,Alloc>& lhs, const std::vector<T,Alloc>& rhs)
{ return lhs.inner_data() <  rhs.inner_data(); }

template <class T, class Alloc>
bool operator<= (const std::vector<T,Alloc>& lhs, const std::vector<T,Alloc>& rhs)
{ return lhs.inner_data() <= rhs.inner_data(); }

template <class T, class Alloc>
bool operator>  (const std::vector<T,Alloc>& lhs, const std::vector<T,Alloc>& rhs)
{ return lhs.inner_data() >  rhs.inner_data(); }

template <class T, class Alloc>
bool operator>= (const std::vector<T,Alloc>& lhs, const std::vector<T,Alloc>& rhs)
{ return lhs.inner_data() >= rhs.inner_data(); }

} // End namespace shm_kmer_matrix_estimator
