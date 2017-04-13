//
// Created by Andrew Bzikadze on 4/10/17.
//

#pragma once

#include <cstddef>
#include <memory>
#include "verify.hpp"
#include "simulation_routines.hpp"

namespace ig_simulator {

template<class KeyType=size_t, class PriorType=size_t, class FreqType=size_t>
class Treap {
private:
    struct TreapNode;
    // TODO change to unique_ptr
    // using TreapNodePtr = std::shared_ptr<TreapNode>;
    using TreapNodePtr = TreapNode*;

    struct TreapNode {
        // @field key -- index in Tree (unique)
        // @field freq -- the "rational" probability in discrete distribution
        // @field sum -- sum of freq in the subtree possesing @this as a root
        // @field prior -- normally random priority for the heap
        KeyType key;
        FreqType freq, sum;
        PriorType prior;

        TreapNodePtr left, right;

        TreapNode(KeyType key, FreqType freq,
                  PriorType prior,
                  TreapNodePtr left = nullptr, TreapNodePtr right = nullptr) :
            key(key), freq(freq), sum(freq),
            prior(prior),
            left(left), right(right)
        { }

        ~TreapNode() {
            delete left;
            delete right;
        }

        static FreqType Sum(const TreapNodePtr &t) {
            if (t != nullptr)
                return t->sum;
            return 0;
        }

        static void Upd(TreapNodePtr &t) {
            if (t != nullptr)
                t->sum = Sum(t->left) + Sum(t->right) + t->freq;
        }
    };

    TreapNodePtr root;
    size_t treap_size;

private:
    static void Merge(TreapNodePtr *pt, TreapNodePtr &l, TreapNodePtr &r) {
        if (l == nullptr)
            *pt = r;
        else if (r == nullptr)
            *pt = l;
        else if (l->prior < r->prior) {
            Merge(&l->right, l->right, r);
            *pt = l;
        } else {
            Merge(&r->left, l, r->left);
            *pt = r;
        }
        TreapNode::Upd(*pt);
    }

    static void Split(TreapNodePtr t, KeyType key, TreapNodePtr *l, TreapNodePtr *r) {
        if (t == nullptr) {
            *l = nullptr;
            *r = nullptr;
        }
        else if (t->key < key) {
            Split(t->right, key, &t->right, r);
            *l = t;
            TreapNode::Upd(*l);
        } else {
            Split(t->left, key, l, &t->left);
            *r = t;
            TreapNode::Upd(*r);
        }
    }

public:
    Treap(): root(nullptr), treap_size(0) { }
    ~Treap() { delete root; }

    void Insert(KeyType key, FreqType freq, PriorType prior = random_index()) {
        TreapNodePtr * pt = &root;
        while (*pt and (*pt)->prior < prior) {
            (*pt)->sum += freq;
            if (key < (*pt)->key)
                pt = &(*pt)->left;
            else
                pt = &(*pt)->right;
        }
        TreapNodePtr l, r;
        Split(*pt, key, &l, &r);
        *pt = TreapNodePtr(new TreapNode(key, freq, prior, l, r));
        TreapNode::Upd(*pt);
        treap_size++;
    }

    void Erase(KeyType key, FreqType freq) {
        TreapNodePtr * pt = &root;
        while ((*pt)->key != key) {
            (*pt)->sum -= freq;
            if (key < (*pt)->key)
                pt = &(*pt)->left;
            else
                pt = &(*pt)->right;
        }
        TreapNodePtr p;
        Merge(&p, (*pt)->left, (*pt)->right);
        (*pt)->left = nullptr;
        (*pt)->right = nullptr;
        // TODO Fix bug with not setting to nullptr pointer of parent of *pt if *pt has no children
        delete *pt;
        *pt = p;
        treap_size--;
    }

    KeyType Find(FreqType sum) const {
        TreapNodePtr t = root;
        FreqType temp;
        while((temp = TreapNode::Sum(t->right) + 1) != sum) {
            if (temp > sum)
                t = t->right;
            else {
                t = t->left;
                sum -= temp;
            }
        }
        return t->key;
    }

    void SetFreq(KeyType key, FreqType old_freq, FreqType new_freq) {
        TreapNodePtr t = root;
        while(t->key != key) {
            // if FreqType is unsigned `new_freq - old_freq` is dangerous
            // std::cout << t->sum << " " << old_freq << "\n";
            VERIFY(t->sum >= old_freq);
            t->sum = t->sum - old_freq + new_freq;
            if (t->key > key)
                t = t->left;
            else
                t = t->right;
        }
        t->freq = new_freq;
        t->sum = t->sum - old_freq + new_freq;
        // TreapNode::Upd(t);
    }

    std::pair<KeyType, FreqType> LowerBound(double sum) const {
        TreapNodePtr t = root;
        FreqType sum_left, sum_right;

        while(true) {
            sum_left = TreapNode::Sum(t->left);
            sum_right = TreapNode::Sum(t->right);

            if (sum_left + sum_right < sum)
                break;

            if (sum_left > sum )
                t = t->left;
            else {
                t = t->right;
                sum -= sum_left;
            }
        }
        return { t->key, t->freq };
    }

    FreqType Sum() const {
        return TreapNode::Sum(root);
    }

size_t Size() const { return treap_size; }
};

} // End namespace ig_simulator
