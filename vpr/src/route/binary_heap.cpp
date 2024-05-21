#include <iostream>
#include "binary_heap.h"
#include "rr_graph_fwd.h"
#include "vtr_log.h"
#include "vtr_assert.h"

#if defined(VPR_ROUTER_USE_CUSTOMIZED_BINARY_HEAP)

static size_t parent(size_t i) { return i >> 1; }
// child indices of a heap
static size_t left(size_t i) { return i << 1; }
static size_t right(size_t i) { return (i << 1) + 1; }

BinaryHeap::BinaryHeap()
    : heap_()
    , heap_size_(0)
    , heap_tail_(0)
    , max_index_(std::numeric_limits<size_t>::max())
    , prune_limit_(std::numeric_limits<size_t>::max()) {
    VTR_LOG("VPR router use customized binary heap\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    size_t target_heap_size = (grid.width() - 1) * (grid.height() - 1);
    if (heap_.empty() || heap_size_ < target_heap_size) {
        if (!heap_.empty()) {
            // coverity[offset_free : Intentional]
            heap_.clear();
        }
        heap_size_ = (grid.width() - 1) * (grid.height() - 1);
        heap_.resize(heap_size_ + 1); /* heap_size_ + 1 because heap stores from [1..heap_size] */
    }
    heap_tail_ = 1;
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    expand_heap_if_full();
    // start with undefined hole
    ++heap_tail_;
    sift_up(heap_tail_ - 1, hptr);

    // If we have pruned, rebuild the heap now.
    if (check_prune_limit()) {
        build_heap();
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_tail_ == 1);
}

t_heap* BinaryHeap::get_heap_head() {
    /* Returns a pointer to the smallest element on the heap, or NULL if the     *
     * heap is empty.  Invalid (index == OPEN) entries on the heap are never     *
     * returned -- they are just skipped over.                                   */

    t_heap* cheapest;
    size_t hole, child;

    do {
        if (heap_tail_ == 1) { /* Empty heap. */
            VTR_LOG_WARN("Empty heap occurred in get_heap_head.\n");
            return (nullptr);
        }

        cheapest = heap_[1];

        hole = 1;
        child = 2;
        --heap_tail_;
        while (child < heap_tail_) {
            if (heap_[child + 1]->cost < heap_[child]->cost)
                ++child; // become right child
            heap_[hole] = heap_[child];
            hole = child;
            child = left(child);
        }
        sift_up(hole, heap_[heap_tail_]);

    } while (!cheapest->index.is_valid()); /* Get another one if invalid entry. */

    return (cheapest);
}

void BinaryHeap::empty_heap() {
    for (size_t i = 1; i < heap_tail_; i++)
        free(heap_[i]);

    heap_tail_ = 1;
}

size_t BinaryHeap::size() const { return heap_tail_ - 1; } // heap[0] is not valid element

// make a heap rooted at index hole by **sifting down** in O(lgn) time
void BinaryHeap::sift_down(size_t hole) {
    t_heap* head{heap_[hole]};
    size_t child{left(hole)};
    while (child < heap_tail_) {
        if (child + 1 < heap_tail_ && heap_[child + 1]->cost < heap_[child]->cost)
            ++child;
        if (heap_[child]->cost < head->cost) {
            heap_[hole] = heap_[child];
            hole = child;
            child = left(child);
        } else
            break;
    }
    heap_[hole] = head;
}

// runs in O(n) time by sifting down; the least work is done on the most elements: 1 swap for bottom layer, 2 swap for 2nd, ... lgn swap for top
// 1*(n/2) + 2*(n/4) + 3*(n/8) + ... + lgn*1 = 2n (sum of i/2^i)
void BinaryHeap::build_heap() {
    // second half of heap are leaves
    for (size_t i = heap_tail_ >> 1; i != 0; --i)
        sift_down(i);
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    if (prune_limit != std::numeric_limits<size_t>::max()) {
        VTR_ASSERT(max_index < prune_limit);
    }
    max_index_ = max_index;
    prune_limit_ = prune_limit;
}

// O(lgn) sifting up to maintain heap property after insertion (should sift down when building heap)
void BinaryHeap::sift_up(size_t leaf, t_heap* const node) {
    while ((leaf > 1) && (node->cost < heap_[parent(leaf)]->cost)) {
        // sift hole up
        heap_[leaf] = heap_[parent(leaf)];
        leaf = parent(leaf);
    }
    heap_[leaf] = node;
}

//expands heap by "realloc"
void BinaryHeap::expand_heap_if_full() {
    if (heap_tail_ > heap_size_) { /* Heap is full */
        heap_size_ *= 2;
        heap_.resize(heap_size_ + 1);
    }
}

// adds an element to the back of heap and expand if necessary, but does not maintain heap property
void BinaryHeap::push_back(t_heap* const hptr) {
    expand_heap_if_full();
    heap_[heap_tail_] = hptr;
    ++heap_tail_;

    check_prune_limit();
}

bool BinaryHeap::is_valid() const {
    if (heap_.empty()) {
        return false;
    }

    for (size_t i = 1; i <= heap_tail_ >> 1; ++i) {
        if (left(i) < heap_tail_ && heap_[left(i)]->cost < heap_[i]->cost) return false;
        if (right(i) < heap_tail_ && heap_[right(i)]->cost < heap_[i]->cost) return false;
    }
    return true;
}

void BinaryHeap::free_all_memory() {
    if (!heap_.empty()) {
        empty_heap();

        // coverity[offset_free : Intentional]
        heap_.clear();
    }

    //  heap_ = nullptr; /* Defensive coding:  crash hard if I use these. */

    storage_.free_all_memory();
}

bool BinaryHeap::check_prune_limit() {
    if (heap_tail_ > prune_limit_) {
        prune_heap();
        return true;
    }

    return false;
}

void BinaryHeap::prune_heap() {
    VTR_ASSERT(max_index_ < prune_limit_);

    std::vector<t_heap*> best_heap_item(max_index_, nullptr);

    // Find the cheapest instance of each index and store it.
    for (size_t i = 1; i < heap_tail_; i++) {
        if (heap_[i] == nullptr) {
            continue;
        }

        if (!heap_[i]->index.is_valid()) {
            free(heap_[i]);
            heap_[i] = nullptr;
            continue;
        }

        auto idx = size_t(heap_[i]->index);

        VTR_ASSERT(idx < max_index_);

        if (best_heap_item[idx] == nullptr || best_heap_item[idx]->cost > heap_[i]->cost) {
            best_heap_item[idx] = heap_[i];
        }
    }

    // Free unused nodes.
    for (size_t i = 1; i < heap_tail_; i++) {
        if (heap_[i] == nullptr) {
            continue;
        }

        auto idx = size_t(heap_[i]->index);

        if (best_heap_item[idx] != heap_[i]) {
            free(heap_[i]);
            heap_[i] = nullptr;
        }
    }

    heap_tail_ = 1;

    for (size_t i = 0; i < max_index_; ++i) {
        if (best_heap_item[i] != nullptr) {
            heap_[heap_tail_++] = best_heap_item[i];
        }
    }
}

#elif defined(VPR_ROUTER_USE_CUSTOMIZED_BINARY_HEAP_FIRST_INDEX_ZERO)

static size_t parent(size_t i) { return (i - 1) / 2; }
// child indices of a heap
static size_t left(size_t i) { return (i * 2) + 1; }
static size_t right(size_t i) { return (i * 2) + 2; }

BinaryHeap::BinaryHeap()
    : heap_()
    , heap_size_(0)
    , heap_tail_(0)
    , max_index_(std::numeric_limits<size_t>::max())
    , prune_limit_(std::numeric_limits<size_t>::max()) {
    VTR_LOG("VPR router use customized binary heap\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    size_t target_heap_size = (grid.width() - 1) * (grid.height() - 1);
    if (heap_.empty() || heap_size_ < target_heap_size) {
        if (!heap_.empty()) {
            // coverity[offset_free : Intentional]
            heap_.clear();
        }
        heap_size_ = (grid.width() - 1) * (grid.height() - 1);
        heap_.resize(heap_size_);
    }
    heap_tail_ = 0;
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    expand_heap_if_full();
    // start with undefined hole
    ++heap_tail_;
    sift_up(heap_tail_ - 1, hptr);

    // If we have pruned, rebuild the heap now.
    if (check_prune_limit()) {
        build_heap();
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_tail_ == 0);
}

t_heap* BinaryHeap::get_heap_head() {
    /* Returns a pointer to the smallest element on the heap, or NULL if the     *
     * heap is empty.  Invalid (index == OPEN) entries on the heap are never     *
     * returned -- they are just skipped over.                                   */

    t_heap* cheapest;
    size_t hole, child;

    do {
        if (heap_tail_ == 0) { /* Empty heap. */
            VTR_LOG_WARN("Empty heap occurred in get_heap_head.\n");
            return (nullptr);
        }

        cheapest = heap_[0];

        hole = 0;
        child = 1;
        --heap_tail_;
        while (child < heap_tail_) {
            if (heap_[child + 1]->cost < heap_[child]->cost)
                ++child; // become right child
            heap_[hole] = heap_[child];
            hole = child;
            child = left(child);
        }
        sift_up(hole, heap_[heap_tail_]);

    } while (!cheapest->index.is_valid()); /* Get another one if invalid entry. */

    return (cheapest);
}

void BinaryHeap::empty_heap() {
    for (size_t i = 0; i < heap_tail_; i++)
        free(heap_[i]);

    heap_tail_ = 0;
}

size_t BinaryHeap::size() const { return heap_tail_; }

// make a heap rooted at index hole by **sifting down** in O(lgn) time
void BinaryHeap::sift_down(size_t hole) {
    t_heap* head{heap_[hole]};
    size_t child{left(hole)};
    while (child < heap_tail_) {
        if (child + 1 < heap_tail_ && heap_[child + 1]->cost < heap_[child]->cost)
            ++child;
        if (heap_[child]->cost < head->cost) {
            heap_[hole] = heap_[child];
            hole = child;
            child = left(child);
        } else
            break;
    }
    heap_[hole] = head;
}

// runs in O(n) time by sifting down; the least work is done on the most elements: 1 swap for bottom layer, 2 swap for 2nd, ... lgn swap for top
// 1*(n/2) + 2*(n/4) + 3*(n/8) + ... + lgn*1 = 2n (sum of i/2^i)
void BinaryHeap::build_heap() {
    // second half of heap are leaves
    for (int i = (int)(heap_tail_ >> 1); i != -1; --i)
        sift_down(i);
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    if (prune_limit != std::numeric_limits<size_t>::max()) {
        VTR_ASSERT(max_index < prune_limit);
    }
    max_index_ = max_index;
    prune_limit_ = prune_limit;
}

// O(lgn) sifting up to maintain heap property after insertion (should sift down when building heap)
void BinaryHeap::sift_up(size_t leaf, t_heap* const node) {
    while ((leaf > 0) && (node->cost < heap_[parent(leaf)]->cost)) {
        // sift hole up
        heap_[leaf] = heap_[parent(leaf)];
        leaf = parent(leaf);
    }

    heap_[leaf] = node;
}

//expands heap by "realloc"
void BinaryHeap::expand_heap_if_full() {
    if (heap_tail_ >= heap_size_) { /* Heap is full */
        heap_size_ *= 2;
        heap_.resize(heap_size_);
    }
}

// adds an element to the back of heap and expand if necessary, but does not maintain heap property
void BinaryHeap::push_back(t_heap* const hptr) {
    expand_heap_if_full();
    heap_[heap_tail_] = hptr;
    ++heap_tail_;

    check_prune_limit();
}

bool BinaryHeap::is_valid() const {
    if (heap_.empty()) {
        return false;
    }

    for (size_t i = 0; i <= heap_tail_ >> 1; ++i) {
        if (left(i) < heap_tail_ && heap_[left(i)]->cost < heap_[i]->cost) return false;
        if (right(i) < heap_tail_ && heap_[right(i)]->cost < heap_[i]->cost) return false;
    }
    return true;
}

void BinaryHeap::free_all_memory() {
    if (!heap_.empty()) {
        empty_heap();

        // coverity[offset_free : Intentional]
        heap_.clear();
    }

    //  heap_ = nullptr; /* Defensive coding:  crash hard if I use these. */

    storage_.free_all_memory();
}

bool BinaryHeap::check_prune_limit() {
    if (heap_tail_ > prune_limit_) {
        prune_heap();
        return true;
    }

    return false;
}

void BinaryHeap::prune_heap() {
    VTR_ASSERT(max_index_ < prune_limit_);

    std::vector<t_heap*> best_heap_item(max_index_, nullptr);

    // Find the cheapest instance of each index and store it.
    for (size_t i = 0; i < heap_tail_; i++) {
        if (heap_[i] == nullptr) {
            continue;
        }

        if (!heap_[i]->index.is_valid()) {
            free(heap_[i]);
            heap_[i] = nullptr;
            continue;
        }

        auto idx = size_t(heap_[i]->index);

        VTR_ASSERT(idx < max_index_);

        if (best_heap_item[idx] == nullptr || best_heap_item[idx]->cost > heap_[i]->cost) {
            best_heap_item[idx] = heap_[i];
        }
    }

    // Free unused nodes.
    for (size_t i = 0; i < heap_tail_; i++) {
        if (heap_[i] == nullptr) {
            continue;
        }

        auto idx = size_t(heap_[i]->index);

        if (best_heap_item[idx] != heap_[i]) {
            free(heap_[i]);
            heap_[i] = nullptr;
        }
    }

    heap_tail_ = 0;

    for (size_t i = 0; i < max_index_; ++i) {
        if (best_heap_item[i] != nullptr) {
            heap_[heap_tail_++] = best_heap_item[i];
        }
    }
}

#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP)

static inline size_t parent(size_t i) { return (i + 2) >> 2; }
static inline size_t first_child(size_t i) { return (i << 2) - 2; } // weird expressions because first index is 1

size_t BinaryHeap::smallest_child(size_t i) {
    // Returns heap_tail_ if i has no children

    size_t child_1 = first_child(i);
    size_t child_2 = child_1 + 1;
    size_t child_3 = child_1 + 2;
    size_t child_4 = child_1 + 3;

    int num_children = (((int)heap_tail_ - (int)child_1) > 4) ? 4 : (int)heap_tail_ - (int)child_1;

    switch (num_children) {
        case 4: {
            size_t minA = (heap_[child_1]->cost < heap_[child_2]->cost) ? child_1 : child_2;
            size_t minB = (heap_[child_3]->cost < heap_[child_4]->cost) ? child_3 : child_4;
            return (heap_[minA]->cost < heap_[minB]->cost) ? minA : minB;
        }
        case 3: {
            size_t minA = (heap_[child_1]->cost < heap_[child_2]->cost) ? child_1 : child_2;
            return (heap_[minA]->cost < heap_[child_3]->cost) ? minA : child_3;
        }
        case 2:
            return (heap_[child_1]->cost < heap_[child_2]->cost) ? child_1 : child_2;
        case 1:
            return child_1;
        default:
            return heap_tail_;
    }
}

BinaryHeap::BinaryHeap()
    : heap_()
    , heap_size_(0)
    , heap_tail_(0)
    , max_index_(std::numeric_limits<size_t>::max())
    , prune_limit_(std::numeric_limits<size_t>::max()) {
    VTR_LOG("VPR router use customized 4-ary heap\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    size_t target_heap_size = (grid.width() - 1) * (grid.height() - 1);
    if (heap_.empty() || heap_size_ < target_heap_size) {
        if (!heap_.empty()) {
            // coverity[offset_free : Intentional]
            heap_.clear();
        }
        heap_size_ = (grid.width() - 1) * (grid.height() - 1);
        heap_.resize(heap_size_ + 1); /* heap_size_ + 1 because heap stores from [1..heap_size] */
    }
    heap_tail_ = 1;
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    expand_heap_if_full();
    // start with undefined hole
    ++heap_tail_;
    sift_up(heap_tail_ - 1, hptr);

    // If we have pruned, rebuild the heap now.
    if (check_prune_limit()) {
        build_heap();
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_tail_ == 1);
}

t_heap* BinaryHeap::get_heap_head() {
    /* Returns a pointer to the smallest element on the heap, or NULL if the     *
     * heap is empty.  Invalid (index == OPEN) entries on the heap are never     *
     * returned -- they are just skipped over.                                   */

    t_heap* cheapest;
    size_t hole, child;

    do {
        if (heap_tail_ == 1) { /* Empty heap. */
            VTR_LOG_WARN("Empty heap occurred in get_heap_head.\n");
            return (nullptr);
        }

        cheapest = heap_[1];

        hole = 1;
        child = smallest_child(hole);
        --heap_tail_;

        while (child < heap_tail_) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        }

        sift_up(hole, heap_[heap_tail_]);
    } while (!cheapest->index.is_valid()); /* Get another one if invalid entry. */

    return (cheapest);
}

void BinaryHeap::empty_heap() {
    for (size_t i = 1; i < heap_tail_; i++)
        free(heap_[i]);

    heap_tail_ = 1;
}

size_t BinaryHeap::size() const { return heap_tail_ - 1; } // heap[0] is not valid element

// make a heap rooted at index hole by **sifting down** in O(lgn) time
void BinaryHeap::sift_down(size_t hole) {
    t_heap* head{heap_[hole]};
    size_t child{smallest_child(hole)};
    while (child < heap_tail_) {
        if (heap_[child]->cost < head->cost) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        } else
            break;
    }
    heap_[hole] = head;
}

// runs in O(n) time by sifting down; the least work is done on the most elements: 1 swap for bottom layer, 2 swap for 2nd, ... lgn swap for top
// 1*(n/4) + 2*(n/16) + 3*(n/64) + ... + lgn*1 = (4/9) * n (sum of i/4^i)
void BinaryHeap::build_heap() {
    // start at highest index branch node
    for (size_t i = parent(heap_tail_); i != 0; --i)
        sift_down(i);
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    if (prune_limit != std::numeric_limits<size_t>::max()) {
        VTR_ASSERT(max_index < prune_limit);
    }
    max_index_ = max_index;
    prune_limit_ = prune_limit;
}

// O(lgn) sifting up to maintain heap property after insertion (should sift down when building heap)
void BinaryHeap::sift_up(size_t leaf, t_heap* const node) {
    while ((leaf > 1) && (node->cost < heap_[parent(leaf)]->cost)) {
        // sift hole up
        heap_[leaf] = heap_[parent(leaf)];
        leaf = parent(leaf);
    }
    heap_[leaf] = node;
}

// expands heap by "realloc"
void BinaryHeap::expand_heap_if_full() {
    if (heap_tail_ > heap_size_) { /* Heap is full */
        heap_size_ *= 2;
        heap_.resize(heap_size_ + 1);
    }
}

// adds an element to the back of heap and expand if necessary, but does not maintain heap property
void BinaryHeap::push_back(t_heap* const hptr) {
    expand_heap_if_full();
    heap_[heap_tail_] = hptr;
    ++heap_tail_;

    check_prune_limit();
}

bool BinaryHeap::is_valid() const {
    if (heap_.empty()) {
        return false;
    }

    for (size_t i = 1; i <= parent(heap_tail_); ++i) {
        size_t leftmost_child = first_child(i);

        for (size_t j = 0; j < 4; ++j) {
            if (leftmost_child + j >= heap_tail_)
                break;
            else if (heap_[leftmost_child + j]->cost < heap_[i]->cost)
                return false;
        }
    }

    return true;
}

void BinaryHeap::free_all_memory() {
    if (!heap_.empty()) {
        empty_heap();

        // coverity[offset_free : Intentional]
        heap_.clear();
    }

    //  heap_ = nullptr; /* Defensive coding:  crash hard if I use these. */

    storage_.free_all_memory();
}

bool BinaryHeap::check_prune_limit() {
    if (heap_tail_ > prune_limit_) {
        prune_heap();
        return true;
    }

    return false;
}

void BinaryHeap::prune_heap() {
    VTR_ASSERT(max_index_ < prune_limit_);

    std::vector<t_heap*> best_heap_item(max_index_, nullptr);

    // Find the cheapest instance of each index and store it.
    for (size_t i = 1; i < heap_tail_; i++) {
        if (heap_[i] == nullptr) {
            continue;
        }

        if (!heap_[i]->index.is_valid()) {
            free(heap_[i]);
            heap_[i] = nullptr;
            continue;
        }

        auto idx = size_t(heap_[i]->index);

        VTR_ASSERT(idx < max_index_);

        if (best_heap_item[idx] == nullptr || best_heap_item[idx]->cost > heap_[i]->cost) {
            best_heap_item[idx] = heap_[i];
        }
    }

    // Free unused nodes.
    for (size_t i = 1; i < heap_tail_; i++) {
        if (heap_[i] == nullptr) {
            continue;
        }

        auto idx = size_t(heap_[i]->index);

        if (best_heap_item[idx] != heap_[i]) {
            free(heap_[i]);
            heap_[i] = nullptr;
        }
    }

    heap_tail_ = 1;

    for (size_t i = 0; i < max_index_; ++i) {
        if (best_heap_item[i] != nullptr) {
            heap_[heap_tail_++] = best_heap_item[i];
        }
    }
}

#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FIRST_INDEX_ZERO)

static size_t parent(size_t i) { return i >> 2; }
static size_t first_child(size_t i) { return (i != 0) ? (i << 2) : 1; }

size_t BinaryHeap::smallest_child(size_t i) {
    // Returns heap_tail_ if i has no children

    size_t child_1 = first_child(i);
    size_t child_2 = child_1 + 1;
    size_t child_3 = child_1 + 2;
    size_t child_4 = child_1 + 3;

    int num_children = (((int)heap_tail_ - (int)child_1) > 4) ? 4 : (int)heap_tail_ - (int)child_1;
    num_children = (i == 0 && num_children > 3) ? 3 : num_children; // the second layer of the heap has only 3 nodes

    switch (num_children) {
        case 4: {
            size_t minA = (heap_[child_1].elem_ptr->cost < heap_[child_2].elem_ptr->cost) ? child_1 : child_2;
            size_t minB = (heap_[child_3].elem_ptr->cost < heap_[child_4].elem_ptr->cost) ? child_3 : child_4;
            return (heap_[minA].elem_ptr->cost < heap_[minB].elem_ptr->cost) ? minA : minB;
        }
        case 3: {
            size_t minA = (heap_[child_1].elem_ptr->cost < heap_[child_2].elem_ptr->cost) ? child_1 : child_2;
            return (heap_[minA].elem_ptr->cost < heap_[child_3].elem_ptr->cost) ? minA : child_3;
        }
        case 2:
            return (heap_[child_1].elem_ptr->cost < heap_[child_2].elem_ptr->cost) ? child_1 : child_2;
        case 1:
            return child_1;
        default:
            return heap_tail_;
    }
}

BinaryHeap::BinaryHeap()
    : heap_()
    , heap_size_(0)
    , heap_tail_(0)
    , max_index_(std::numeric_limits<size_t>::max())
    , prune_limit_(std::numeric_limits<size_t>::max()) {
    VTR_LOG("VPR router use customized 4-ary heap\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    size_t target_heap_size = (grid.width() - 1) * (grid.height() - 1);
    if (heap_.empty() || heap_size_ < target_heap_size) {
        if (!heap_.empty()) {
            // coverity[offset_free : Intentional]
            heap_.clear();
        }
        heap_size_ = (grid.width() - 1) * (grid.height() - 1);
        heap_.resize(heap_size_);
    }
    heap_tail_ = 0;
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    expand_heap_if_full();
    // start with undefined hole
    ++heap_tail_;
    sift_up(heap_tail_ - 1, hptr);

    // If we have pruned, rebuild the heap now.
    if (check_prune_limit()) {
        build_heap();
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_tail_ == 0);
}

t_heap* BinaryHeap::get_heap_head() {
    /* Returns a pointer to the smallest element on the heap, or NULL if the     *
     * heap is empty.  Invalid (index == OPEN) entries on the heap are never     *
     * returned -- they are just skipped over.                                   */

    t_heap* cheapest;
    size_t hole, child;

    do {
        if (heap_tail_ == 0) { /* Empty heap. */
            VTR_LOG_WARN("Empty heap occurred in get_heap_head.\n");
            return (nullptr);
        }

        cheapest = heap_[0].elem_ptr;

        hole = 0;
        child = smallest_child(hole);
        --heap_tail_;

        while (child < heap_tail_) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        }

        sift_up(hole, heap_[heap_tail_].elem_ptr);
    } while (!cheapest->index.is_valid()); /* Get another one if invalid entry. */

    return (cheapest);
}

void BinaryHeap::empty_heap() {
    for (size_t i = 0; i < heap_tail_; i++)
        free(heap_[i].elem_ptr);

    heap_tail_ = 0;
}

size_t BinaryHeap::size() const { return heap_tail_; }

// make a heap rooted at index hole by **sifting down** in O(lgn) time
void BinaryHeap::sift_down(size_t hole) {
    t_heap* head{heap_[hole].elem_ptr};
    size_t child{smallest_child(hole)};
    while (child < heap_tail_) {
        if (heap_[child].elem_ptr->cost < head->cost) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        } else
            break;
    }
    heap_[hole].elem_ptr = head;
}

// runs in O(n) time by sifting down; the least work is done on the most elements: 1 swap for bottom layer, 2 swap for 2nd, ... lgn swap for top
// 1*(n/4) + 2*(n/16) + 3*(n/64) + ... + lgn*1 = (4/9) * n (sum of i/4^i)
void BinaryHeap::build_heap() {
    // start at highest index branch node
    for (int i = (int)parent(heap_tail_); i != -1; --i)
        sift_down(i);
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    if (prune_limit != std::numeric_limits<size_t>::max()) {
        VTR_ASSERT(max_index < prune_limit);
    }
    max_index_ = max_index;
    prune_limit_ = prune_limit;
}

// O(lgn) sifting up to maintain heap property after insertion (should sift down when building heap)
void BinaryHeap::sift_up(size_t leaf, t_heap* const node) {
    while ((leaf > 0) && (node->cost < heap_[parent(leaf)].elem_ptr->cost)) {
        // sift hole up
        heap_[leaf] = heap_[parent(leaf)];
        leaf = parent(leaf);
    }
    heap_[leaf].elem_ptr = node;
}

// expands heap by "realloc"
void BinaryHeap::expand_heap_if_full() {
    if (heap_tail_ >= heap_size_) { /* Heap is full */
        heap_size_ *= 2;
        heap_.resize(heap_size_);
    }
}

// adds an element to the back of heap and expand if necessary, but does not maintain heap property
void BinaryHeap::push_back(t_heap* const hptr) {
    expand_heap_if_full();
    heap_[heap_tail_].elem_ptr = hptr;
    ++heap_tail_;

    check_prune_limit();
}

bool BinaryHeap::is_valid() const {
    if (heap_.empty()) {
        return false;
    }

    for (size_t i = 0; i <= parent(heap_tail_); ++i) {
        size_t leftmost_child = first_child(i);
        size_t num_potential_children = (i == 0) ? 3 : 4; // second layer in heap has 3 nodes

        for (size_t j = 0; j < num_potential_children; ++j) {
            if (leftmost_child + j >= heap_tail_)
                break;
            else if (heap_[leftmost_child + j].elem_ptr->cost < heap_[i].elem_ptr->cost)
                return false;
        }
    }

    return true;
}

void BinaryHeap::free_all_memory() {
    if (!heap_.empty()) {
        empty_heap();

        // coverity[offset_free : Intentional]
        heap_.clear();
    }

    //  heap_ = nullptr; /* Defensive coding:  crash hard if I use these. */

    storage_.free_all_memory();
}

bool BinaryHeap::check_prune_limit() {
    if (heap_tail_ > prune_limit_) {
        prune_heap();
        return true;
    }

    return false;
}

void BinaryHeap::prune_heap() {
    VTR_ASSERT(max_index_ < prune_limit_);

    std::vector<t_heap*> best_heap_item(max_index_, nullptr);

    // Find the cheapest instance of each index and store it.
    for (size_t i = 0; i < heap_tail_; i++) {
        if (heap_[i].elem_ptr == nullptr) {
            continue;
        }

        if (!heap_[i].elem_ptr->index.is_valid()) {
            free(heap_[i].elem_ptr);
            heap_[i].elem_ptr = nullptr;
            continue;
        }

        auto idx = size_t(heap_[i].elem_ptr->index);

        VTR_ASSERT(idx < max_index_);

        if (best_heap_item[idx] == nullptr || best_heap_item[idx]->cost > heap_[i].elem_ptr->cost) {
            best_heap_item[idx] = heap_[i].elem_ptr;
        }
    }

    // Free unused nodes.
    for (size_t i = 0; i < heap_tail_; i++) {
        if (heap_[i].elem_ptr == nullptr) {
            continue;
        }

        auto idx = size_t(heap_[i].elem_ptr->index);

        if (best_heap_item[idx] != heap_[i].elem_ptr) {
            free(heap_[i].elem_ptr);
            heap_[i].elem_ptr = nullptr;
        }
    }

    heap_tail_ = 0;

    for (size_t i = 0; i < max_index_; ++i) {
        if (best_heap_item[i] != nullptr) {
            heap_[heap_tail_++].elem_ptr = best_heap_item[i];
        }
    }
}

#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM)

static size_t parent(size_t i) { return i >> 2; }
static size_t first_child(size_t i) { return (i != 0) ? (i << 2) : 1; }

size_t BinaryHeap::smallest_child(size_t i) {
    // Returns heap_tail_ if i has no children

    size_t child_1 = first_child(i);
    size_t child_2 = child_1 + 1;
    size_t child_3 = child_1 + 2;
    size_t child_4 = child_1 + 3;

    int num_children = (((int)heap_tail_ - (int)child_1) > 4) ? 4 : (int)heap_tail_ - (int)child_1;
    num_children = (i == 0 && num_children > 3) ? 3 : num_children; // the second layer of the heap has only 3 nodes

    switch (num_children) {
        case 4: {
            size_t minA = (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
            size_t minB = (heap_[child_3].cost < heap_[child_4].cost) ? child_3 : child_4;
            return (heap_[minA].cost < heap_[minB].cost) ? minA : minB;
        }
        case 3: {
            size_t minA = (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
            return (heap_[minA].cost < heap_[child_3].cost) ? minA : child_3;
        }
        case 2:
            return (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
        case 1:
            return child_1;
        default:
            return heap_tail_;
    }
}

BinaryHeap::BinaryHeap()
    : heap_()
    , heap_size_(0)
    , heap_tail_(0)
    , max_index_(std::numeric_limits<size_t>::max())
    , prune_limit_(std::numeric_limits<size_t>::max()) {
    VTR_LOG("VPR router use customized 4-ary heap with new heap_elem struct\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    size_t target_heap_size = (grid.width() - 1) * (grid.height() - 1);
    if (heap_.empty() || heap_size_ < target_heap_size) {
        if (!heap_.empty()) {
            // coverity[offset_free : Intentional]
            heap_.clear();
        }
        heap_size_ = (grid.width() - 1) * (grid.height() - 1);
        heap_.resize(heap_size_);
    }
    heap_tail_ = 0;
}

void BinaryHeap::add_to_heap_non_virtual(heap_elem elem) {
    expand_heap_if_full();
    // start with undefined hole
    ++heap_tail_;
    sift_up(heap_tail_ - 1, elem);

    // If we have pruned, rebuild the heap now.
    if (check_prune_limit()) {
        build_heap();
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_tail_ == 0);
}

uint32_t BinaryHeap::get_heap_head_non_virtual() {
    /* Returns a pointer to the smallest element on the heap, or NULL if the     *
     * heap is empty.  Invalid (index == OPEN) entries on the heap are never     *
     * returned -- they are just skipped over.                                   */

    uint32_t cheapest;
    size_t hole, child;

    do {
        if (heap_tail_ == 0) { /* Empty heap. */
            VTR_LOG_WARN("Empty heap occurred in get_heap_head.\n");
            return (0);
        }

        cheapest = heap_[0].ID;

        hole = 0;
        child = smallest_child(hole);
        --heap_tail_;

        while (child < heap_tail_) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        }

        sift_up(hole, heap_[heap_tail_]);
    } while (!cheapest->index.is_valid()); /* Get another one if invalid entry. */

    return (cheapest);
}

void BinaryHeap::empty_heap() {
    //    for (size_t i = 0; i < heap_tail_; i++)
    //        free(heap_[i]);

    heap_tail_ = 0;
}

size_t BinaryHeap::size() const { return heap_tail_; }

// make a heap rooted at index hole by **sifting down** in O(lgn) time
void BinaryHeap::sift_down(size_t hole) {
    heap_elem head{heap_[hole]};
    size_t child{smallest_child(hole)};
    while (child < heap_tail_) {
        if (heap_[child].cost < head.cost) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        } else
            break;
    }
    heap_[hole] = head;
}

// runs in O(n) time by sifting down; the least work is done on the most elements: 1 swap for bottom layer, 2 swap for 2nd, ... lgn swap for top
// 1*(n/4) + 2*(n/16) + 3*(n/64) + ... + lgn*1 = (4/9) * n (sum of i/4^i)
void BinaryHeap::build_heap() {
    // start at highest index branch node
    for (int i = (int)parent(heap_tail_); i != -1; --i)
        sift_down(i);
}

//void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
//    if (prune_limit != std::numeric_limits<size_t>::max()) {
//        VTR_ASSERT(max_index < prune_limit);
//    }
//    max_index_ = max_index;
//    prune_limit_ = prune_limit;
//}

// O(lgn) sifting up to maintain heap property after insertion (should sift down when building heap)
void BinaryHeap::sift_up(size_t leaf, heap_elem const& node) {
    while ((leaf > 0) && (node.cost < heap_[parent(leaf)].cost)) {
        // sift hole up
        heap_[leaf] = heap_[parent(leaf)];
        leaf = parent(leaf);
    }
    heap_[leaf] = node;
}

// expands heap by "realloc"
void BinaryHeap::expand_heap_if_full() {
    if (heap_tail_ >= heap_size_) { /* Heap is full */
        heap_size_ *= 2;
        heap_.resize(heap_size_);
    }
}

// adds an element to the back of heap and expand if necessary, but does not maintain heap property
void BinaryHeap::push_back_non_virtual(heap_elem const& elem) {
    expand_heap_if_full();
    heap_[heap_tail_] = elem;
    ++heap_tail_;

    //    check_prune_limit();
}

bool BinaryHeap::is_valid() const {
    if (heap_.empty()) {
        return false;
    }

    for (size_t i = 0; i <= parent(heap_tail_); ++i) {
        size_t leftmost_child = first_child(i);
        size_t num_potential_children = (i == 0) ? 3 : 4; // second layer in heap has 3 nodes

        for (size_t j = 0; j < num_potential_children; ++j) {
            if (leftmost_child + j >= heap_tail_)
                break;
            else if (heap_[leftmost_child + j].cost < heap_[i].cost)
                return false;
        }
    }

    return true;
}

void BinaryHeap::free_all_memory() {
    if (!heap_.empty()) {
        empty_heap();

        // coverity[offset_free : Intentional]
        heap_.clear();
    }

    //  heap_ = nullptr; /* Defensive coding:  crash hard if I use these. */

    storage_.free_all_memory();
}

//bool BinaryHeap::check_prune_limit() {
//    if (heap_tail_ > prune_limit_) {
//        prune_heap();
//        return true;
//    }
//
//    return false;
//}

//void BinaryHeap::prune_heap() {
//    VTR_ASSERT(max_index_ < prune_limit_);
//
//    std::vector<t_heap*> best_heap_item(max_index_, nullptr);
//
//    // Find the cheapest instance of each index and store it.
//    for (size_t i = 0; i < heap_tail_; i++) {
//        if (heap_[i] == nullptr) {
//            continue;
//        }
//
//        if (!heap_[i]->index.is_valid()) {
//            free(heap_[i]);
//            heap_[i] = nullptr;
//            continue;
//        }
//
//        auto idx = size_t(heap_[i]->index);
//
//        VTR_ASSERT(idx < max_index_);
//
//        if (best_heap_item[idx] == nullptr || best_heap_item[idx]->cost > heap_[i].cost) {
//            best_heap_item[idx] = heap_[i];
//        }
//    }
//
//    // Free unused nodes.
//    for (size_t i = 0; i < heap_tail_; i++) {
//        if (heap_[i] == nullptr) {
//            continue;
//        }
//
//        auto idx = size_t(heap_[i]->index);
//
//        if (best_heap_item[idx] != heap_[i]) {
//            free(heap_[i]);
//            heap_[i] = nullptr;
//        }
//    }
//
//    heap_tail_ = 0;
//
//    for (size_t i = 0; i < max_index_; ++i) {
//        if (best_heap_item[i] != nullptr) {
//            heap_[heap_tail_++] = best_heap_item[i];
//        }
//    }
//}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {}
void BinaryHeap::add_to_heap(t_heap* hptr) {}
void BinaryHeap::push_back(t_heap* const hptr) {}
t_heap* BinaryHeap::get_heap_head() {}

#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM_LARGE)

static size_t parent(size_t i) { return i >> 2; }
static size_t first_child(size_t i) { return (i != 0) ? (i << 2) : 1; }

size_t BinaryHeap::smallest_child(size_t i) {
    // Returns heap_tail_ if i has no children

    size_t child_1 = first_child(i);
    size_t child_2 = child_1 + 1;
    size_t child_3 = child_1 + 2;
    size_t child_4 = child_1 + 3;

    int num_children = (((int)heap_tail_ - (int)child_1) > 4) ? 4 : (int)heap_tail_ - (int)child_1;
    num_children = (i == 0 && num_children > 3) ? 3 : num_children; // the second layer of the heap has only 3 nodes

    switch (num_children) {
        case 4: {
            size_t minA = (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
            size_t minB = (heap_[child_3].cost < heap_[child_4].cost) ? child_3 : child_4;
            return (heap_[minA].cost < heap_[minB].cost) ? minA : minB;
        }
        case 3: {
            size_t minA = (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
            return (heap_[minA].cost < heap_[child_3].cost) ? minA : child_3;
        }
        case 2:
            return (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
        case 1:
            return child_1;
        default:
            return heap_tail_;
    }
}

BinaryHeap::BinaryHeap()
    : heap_()
    , heap_size_(0)
    , heap_tail_(0)
    , max_index_(std::numeric_limits<size_t>::max())
    , prune_limit_(std::numeric_limits<size_t>::max()) {
    VTR_LOG("VPR router use customized 4-ary heap with large heap_elem struct\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    size_t target_heap_size = (grid.width() - 1) * (grid.height() - 1);
    if (heap_.empty() || heap_size_ < target_heap_size) {
        if (!heap_.empty()) {
            // coverity[offset_free : Intentional]
            heap_.clear();
        }
        heap_size_ = (grid.width() - 1) * (grid.height() - 1);
        heap_.resize(heap_size_);
    }
    heap_tail_ = 0;
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    expand_heap_if_full();
    // start with undefined hole
    ++heap_tail_;
    heap_elem new_elem = {hptr, hptr->cost};
    sift_up(heap_tail_ - 1, new_elem);

    // If we have pruned, rebuild the heap now.
    if (check_prune_limit()) {
        build_heap();
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_tail_ == 0);
}

t_heap* BinaryHeap::get_heap_head() {
    /* Returns a pointer to the smallest element on the heap, or NULL if the     *
     * heap is empty.  Invalid (index == OPEN) entries on the heap are never     *
     * returned -- they are just skipped over.                                   */

    t_heap* cheapest;
    size_t hole, child;

    do {
        if (heap_tail_ == 0) { /* Empty heap. */
            VTR_LOG_WARN("Empty heap occurred in get_heap_head.\n");
            return nullptr;
        }

        cheapest = heap_[0].elem_ptr;

        hole = 0;
        child = smallest_child(hole);
        --heap_tail_;

        while (child < heap_tail_) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        }

        sift_up(hole, heap_[heap_tail_]);
    } while (!cheapest->index.is_valid()); /* Get another one if invalid entry. */

    return (cheapest);
}

void BinaryHeap::empty_heap() {
    for (size_t i = 0; i < heap_tail_; i++)
        free(heap_[i].elem_ptr);

    heap_tail_ = 0;
}

size_t BinaryHeap::size() const { return heap_tail_; }

// make a heap rooted at index hole by **sifting down** in O(lgn) time
void BinaryHeap::sift_down(size_t hole) {
    heap_elem head{heap_[hole]};
    size_t child{smallest_child(hole)};
    while (child < heap_tail_) {
        if (heap_[child].cost < head.cost) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        } else
            break;
    }
    heap_[hole] = head;
}

// runs in O(n) time by sifting down; the least work is done on the most elements: 1 swap for bottom layer, 2 swap for 2nd, ... lgn swap for top
// 1*(n/4) + 2*(n/16) + 3*(n/64) + ... + lgn*1 = (4/9) * n (sum of i/4^i)
void BinaryHeap::build_heap() {
    // start at highest index branch node
    for (int i = (int)parent(heap_tail_); i != -1; --i)
        sift_down(i);
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    if (prune_limit != std::numeric_limits<size_t>::max()) {
        VTR_ASSERT(max_index < prune_limit);
    }
    max_index_ = max_index;
    prune_limit_ = prune_limit;
}

// O(lgn) sifting up to maintain heap property after insertion (should sift down when building heap)
void BinaryHeap::sift_up(size_t leaf, heap_elem const& node) {
    while ((leaf > 0) && (node.cost < heap_[parent(leaf)].cost)) {
        // sift hole up
        heap_[leaf] = heap_[parent(leaf)];
        leaf = parent(leaf);
    }
    heap_[leaf] = node;
}

// expands heap by "realloc"
void BinaryHeap::expand_heap_if_full() {
    if (heap_tail_ >= heap_size_) { /* Heap is full */
        heap_size_ *= 2;
        heap_.resize(heap_size_);
    }
}

// adds an element to the back of heap and expand if necessary, but does not maintain heap property
void BinaryHeap::push_back(t_heap* const hptr) {
    expand_heap_if_full();
    heap_elem new_elem = {hptr, hptr->cost};
    heap_[heap_tail_] = new_elem;
    ++heap_tail_;

    check_prune_limit();
}

bool BinaryHeap::is_valid() const {
    if (heap_.empty()) {
        return false;
    }

    for (size_t i = 0; i <= parent(heap_tail_); ++i) {
        size_t leftmost_child = first_child(i);
        size_t num_potential_children = (i == 0) ? 3 : 4; // second layer in heap has 3 nodes

        for (size_t j = 0; j < num_potential_children; ++j) {
            if (leftmost_child + j >= heap_tail_)
                break;
            else if (heap_[leftmost_child + j].cost < heap_[i].cost)
                return false;
        }
    }

    return true;
}

void BinaryHeap::free_all_memory() {
    if (!heap_.empty()) {
        empty_heap();

        // coverity[offset_free : Intentional]
        heap_.clear();
    }

    //  heap_ = nullptr; /* Defensive coding:  crash hard if I use these. */

    storage_.free_all_memory();
}

bool BinaryHeap::check_prune_limit() {
    if (heap_tail_ > prune_limit_) {
        prune_heap();
        return true;
    }

    return false;
}

void BinaryHeap::prune_heap() {
    VTR_ASSERT(max_index_ < prune_limit_);

    std::vector<heap_elem> best_heap_item(max_index_, {nullptr, 0.0});

    // Find the cheapest instance of each index and store it.
    for (size_t i = 0; i < heap_tail_; i++) {
        if (heap_[i].elem_ptr == nullptr) {
            continue;
        }

        if (!heap_[i].elem_ptr->index.is_valid()) {
            free(heap_[i].elem_ptr);
            heap_[i].elem_ptr = nullptr;
            continue;
        }

        auto idx = size_t(heap_[i].elem_ptr->index);

        VTR_ASSERT(idx < max_index_);

        if (best_heap_item[idx].elem_ptr == nullptr || best_heap_item[idx].cost > heap_[i].cost) {
            best_heap_item[idx] = heap_[i];
        }
    }

    // Free unused nodes.
    for (size_t i = 0; i < heap_tail_; i++) {
        if (heap_[i].elem_ptr == nullptr) {
            continue;
        }

        auto idx = size_t(heap_[i].elem_ptr->index);

        if (best_heap_item[idx].elem_ptr != heap_[i].elem_ptr) {
            free(heap_[i].elem_ptr);
            heap_[i].elem_ptr = nullptr;
        }
    }

    heap_tail_ = 0;

    for (size_t i = 0; i < max_index_; ++i) {
        if (best_heap_item[i].elem_ptr != nullptr) {
            heap_[heap_tail_++] = best_heap_item[i];
        }
    }
}

#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ARRAY)

static size_t parent(size_t i) { return i >> 2; }
static size_t first_child(size_t i) { return (i != 0) ? (i << 2) : 1; }
static void free_heap(t_heap** heap_) {
    // std::cout << "\tfree_heap" << std::endl;
    free(heap_);
}

size_t BinaryHeap::smallest_child(size_t i) {
    // std::cout << "\tsmallest_child" << std::endl;
    // Returns heap_tail_ if i has no children

    size_t child_1 = first_child(i);
    size_t child_2 = child_1 + 1;
    size_t child_3 = child_1 + 2;
    size_t child_4 = child_1 + 3;

    int num_children = (((int)heap_tail_ - (int)child_1) > 4) ? 4 : (int)heap_tail_ - (int)child_1;
    num_children = (i == 0 && num_children > 3) ? 3 : num_children; // the second layer of the heap has only 3 nodes

    switch (num_children) {
        case 4: {
            size_t minA = (heap_[child_1]->cost < heap_[child_2]->cost) ? child_1 : child_2;
            size_t minB = (heap_[child_3]->cost < heap_[child_4]->cost) ? child_3 : child_4;
            return (heap_[minA]->cost < heap_[minB]->cost) ? minA : minB;
        }
        case 3: {
            size_t minA = (heap_[child_1]->cost < heap_[child_2]->cost) ? child_1 : child_2;
            return (heap_[minA]->cost < heap_[child_3]->cost) ? minA : child_3;
        }
        case 2:
            return (heap_[child_1]->cost < heap_[child_2]->cost) ? child_1 : child_2;
        case 1:
            return child_1;
        default:
            return heap_tail_;
    }
}

BinaryHeap::BinaryHeap()
    : heap_()
    , heap_size_(0)
    , heap_tail_(0)
    , max_index_(std::numeric_limits<size_t>::max())
    , prune_limit_(std::numeric_limits<size_t>::max()) {
    VTR_LOG("VPR router use customized 4-ary heap using array\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    size_t target_heap_size = (grid.width() - 1) * (grid.height() - 1);
    if (size() != 0) {
        // coverity[offset_free : Intentional]
        empty_heap();
        free_heap(heap_);
        // heap_size_ = (grid.width() - 1) * (grid.height() - 1)
    }

    heap_size_ = 1000000000;
    heap_ = (heap_item_t*)aligned_alloc(256, heap_size_ * sizeof(heap_item_t));
    // heap_.resize(heap_size_);

    heap_tail_ = 0;
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    // std::cout << "\tadd_to_heap" << std::endl;
    expand_heap_if_full();
    // start with undefined hole
    ++heap_tail_;
    sift_up(heap_tail_ - 1, hptr);

    // If we have pruned, rebuild the heap now.
    if (check_prune_limit()) {
        build_heap();
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_tail_ == 0);
}

t_heap* BinaryHeap::get_heap_head() {
    // std::cout << "\tget_heap_head" << std::endl;
    /* Returns a pointer to the smallest element on the heap, or NULL if the     *
     * heap is empty.  Invalid (index == OPEN) entries on the heap are never     *
     * returned -- they are just skipped over.                                   */

    t_heap* cheapest;
    size_t hole, child;

    do {
        if (heap_tail_ == 0) { /* Empty heap. */
            VTR_LOG_WARN("Empty heap occurred in get_heap_head.\n");
            return (nullptr);
        }

        cheapest = heap_[0];

        hole = 0;
        child = smallest_child(hole);
        --heap_tail_;

        while (child < heap_tail_) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        }

        sift_up(hole, heap_[heap_tail_]);
    } while (!cheapest->index.is_valid()); /* Get another one if invalid entry. */

    return (cheapest);
}

void BinaryHeap::empty_heap() {
    // std::cout << "\tget_heaempty_heap" << std::endl;
    for (size_t i = 0; i < heap_tail_; i++)
        free(heap_[i]);

    heap_tail_ = 0;
}

size_t BinaryHeap::size() const { return heap_tail_; }

// make a heap rooted at index hole by **sifting down** in O(lgn) time
void BinaryHeap::sift_down(size_t hole) {
    // std::cout << "\tsift_down" << std::endl;
    t_heap* head{heap_[hole]};
    size_t child{smallest_child(hole)};
    while (child < heap_tail_) {
        if (heap_[child]->cost < head->cost) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        } else
            break;
    }
    heap_[hole] = head;
}

// runs in O(n) time by sifting down; the least work is done on the most elements: 1 swap for bottom layer, 2 swap for 2nd, ... lgn swap for top
// 1*(n/4) + 2*(n/16) + 3*(n/64) + ... + lgn*1 = (4/9) * n (sum of i/4^i)
void BinaryHeap::build_heap() {
    // std::cout << "\tbuild_heap" << std::endl;
    // start at highest index branch node
    for (int i = (int)parent(heap_tail_); i != -1; --i)
        sift_down(i);
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    if (prune_limit != std::numeric_limits<size_t>::max()) {
        VTR_ASSERT(max_index < prune_limit);
    }
    max_index_ = max_index;
    prune_limit_ = prune_limit;
}

// O(lgn) sifting up to maintain heap property after insertion (should sift down when building heap)
void BinaryHeap::sift_up(size_t leaf, t_heap* const node) {
    // std::cout << "\tsift_up" << std::endl;
    // std::cout << "\t\tleaf = " << leaf << std::endl;
    // std::cout << "\t\tnode = " << node << std::endl;
    // std::cout << "\t\tnode->cost = " << node << std::endl;
    // std::cout << "\t\tparent(leaf) = " << parent(leaf) << std::endl;
    // std::cout << "\t\theap_[parent(leaf)] = " << heap_[parent(leaf)] << std::endl;
    while ((leaf > 0) && (node->cost < heap_[parent(leaf)]->cost)) {
        // sift hole up
        heap_[leaf] = heap_[parent(leaf)];
        leaf = parent(leaf);
    }
    heap_[leaf] = node;
}

// expands heap by "realloc"
void BinaryHeap::expand_heap_if_full() {
    // std::cout << "\texpand_heap_if_full" << std::endl;
    if (heap_tail_ >= heap_size_) { /* Heap is full */
        heap_ = (heap_item_t*)realloc(heap_, prune_limit_ * sizeof(heap_item_t));
    }
}

// adds an element to the back of heap and expand if necessary, but does not maintain heap property
void BinaryHeap::push_back(t_heap* const hptr) {
    // std::cout << "\tpush_back" << std::endl;
    expand_heap_if_full();
    heap_[heap_tail_] = hptr;
    ++heap_tail_;

    check_prune_limit();
}

bool BinaryHeap::is_valid() const {
    // std::cout << "\tis_valid" << std::endl;
    if (is_empty_heap()) {
        return false;
    }

    for (size_t i = 0; i <= parent(heap_tail_); ++i) {
        size_t leftmost_child = first_child(i);
        size_t num_potential_children = (i == 0) ? 3 : 4; // second layer in heap has 3 nodes

        for (size_t j = 0; j < num_potential_children; ++j) {
            if (leftmost_child + j >= heap_tail_)
                break;
            else if (heap_[leftmost_child + j]->cost < heap_[i]->cost)
                return false;
        }
    }

    return true;
}

void BinaryHeap::free_all_memory() {
    // std::cout << "\tfree_all_memory" << std::endl;
    if (!is_empty_heap()) {
        empty_heap();

        // coverity[offset_free : Intentional]
        free_heap(heap_);
    }

    //  heap_ = nullptr; /* Defensive coding:  crash hard if I use these. */

    storage_.free_all_memory();
}

bool BinaryHeap::check_prune_limit() {
    if (heap_tail_ > prune_limit_) {
        prune_heap();
        return true;
    }

    return false;
}

void BinaryHeap::prune_heap() {
    VTR_ASSERT(max_index_ < prune_limit_);

    std::vector<t_heap*> best_heap_item(max_index_, nullptr);

    // Find the cheapest instance of each index and store it.
    for (size_t i = 0; i < heap_tail_; i++) {
        if (heap_[i] == nullptr) {
            continue;
        }

        if (!heap_[i]->index.is_valid()) {
            free(heap_[i]);
            heap_[i] = nullptr;
            continue;
        }

        auto idx = size_t(heap_[i]->index);

        VTR_ASSERT(idx < max_index_);

        if (best_heap_item[idx] == nullptr || best_heap_item[idx]->cost > heap_[i]->cost) {
            best_heap_item[idx] = heap_[i];
        }
    }

    // Free unused nodes.
    for (size_t i = 0; i < heap_tail_; i++) {
        if (heap_[i] == nullptr) {
            continue;
        }

        auto idx = size_t(heap_[i]->index);

        if (best_heap_item[idx] != heap_[i]) {
            free(heap_[i]);
            heap_[i] = nullptr;
        }
    }

    heap_tail_ = 0;

    for (size_t i = 0; i < max_index_; ++i) {
        if (best_heap_item[i] != nullptr) {
            heap_[heap_tail_++] = best_heap_item[i];
        }
    }
}

#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FINAL)

static inline size_t parent(size_t i) { return (i + 2) >> 2; }
static inline size_t first_child(size_t i) { return (i << 2) - 2; } // weird expressions because first index is 1

size_t BinaryHeap::smallest_child(size_t i) {
    // Returns heap_tail_ if i has no children

    size_t child_1 = first_child(i);
    size_t child_2 = child_1 + 1;
    size_t child_3 = child_1 + 2;
    size_t child_4 = child_1 + 3;

    int num_children = (((int)heap_tail_ - (int)child_1) > 4) ? 4 : (int)heap_tail_ - (int)child_1;

    switch (num_children) {
        case 4: {
            size_t minA = (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
            size_t minB = (heap_[child_3].cost < heap_[child_4].cost) ? child_3 : child_4;
            return (heap_[minA].cost < heap_[minB].cost) ? minA : minB;
        }
        case 3: {
            size_t minA = (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
            return (heap_[minA].cost < heap_[child_3].cost) ? minA : child_3;
        }
        case 2:
            return (heap_[child_1].cost < heap_[child_2].cost) ? child_1 : child_2;
        case 1:
            return child_1;
        default:
            return heap_tail_;
    }
}

BinaryHeap::BinaryHeap()
    : heap_()
    , heap_size_(0)
    , heap_tail_(0)
    , max_index_(std::numeric_limits<size_t>::max())
    , prune_limit_(std::numeric_limits<size_t>::max()) {
    VTR_LOG("VPR router use customized 4-ary heap with large heap_elem struct\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    size_t target_heap_size = (grid.width() - 1) * (grid.height() - 1);
    if (heap_.empty() || heap_size_ < target_heap_size) {
        if (!heap_.empty()) {
            // coverity[offset_free : Intentional]
            heap_.clear();
        }
        heap_size_ = (grid.width() - 1) * (grid.height() - 1);
        heap_.resize(heap_size_ + 1); /* heap_size_ + 1 because heap stores from [1..heap_size] */
    }
    heap_tail_ = 1;
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    expand_heap_if_full();
    // start with undefined hole
    ++heap_tail_;
    heap_elem new_elem = {hptr, hptr->cost};
    sift_up(heap_tail_ - 1, new_elem);

    // If we have pruned, rebuild the heap now.
    if (check_prune_limit()) {
        build_heap();
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_tail_ == 1);
}

t_heap* BinaryHeap::get_heap_head() {
    /* Returns a pointer to the smallest element on the heap, or NULL if the     *
     * heap is empty.  Invalid (index == OPEN) entries on the heap are never     *
     * returned -- they are just skipped over.                                   */

    t_heap* cheapest;
    size_t hole, child;

    do {
        if (heap_tail_ == 1) { /* Empty heap. */
            VTR_LOG_WARN("Empty heap occurred in get_heap_head.\n");
            return nullptr;
        }

        cheapest = heap_[1].elem_ptr;

        hole = 1;
        child = smallest_child(hole);
        --heap_tail_;

        while (child < heap_tail_) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        }

        sift_up(hole, heap_[heap_tail_]);
    } while (!cheapest->index.is_valid()); /* Get another one if invalid entry. */

    return (cheapest);
}

void BinaryHeap::empty_heap() {
    for (size_t i = 1; i < heap_tail_; i++)
        free(heap_[i].elem_ptr);

    heap_tail_ = 1;
}

size_t BinaryHeap::size() const { return heap_tail_; }

// make a heap rooted at index hole by **sifting down** in O(lgn) time
void BinaryHeap::sift_down(size_t hole) {
    heap_elem head{heap_[hole]};
    size_t child{smallest_child(hole)};
    while (child < heap_tail_) {
        if (heap_[child].cost < head.cost) {
            heap_[hole] = heap_[child];
            hole = child;
            child = smallest_child(child);
        } else
            break;
    }
    heap_[hole] = head;
}

// runs in O(n) time by sifting down; the least work is done on the most elements: 1 swap for bottom layer, 2 swap for 2nd, ... lgn swap for top
// 1*(n/4) + 2*(n/16) + 3*(n/64) + ... + lgn*1 = (4/9) * n (sum of i/4^i)
void BinaryHeap::build_heap() {
    // start at highest index branch node
    for (size_t i = (int)parent(heap_tail_); i != 0; --i)
        sift_down(i);
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    if (prune_limit != std::numeric_limits<size_t>::max()) {
        VTR_ASSERT(max_index < prune_limit);
    }
    max_index_ = max_index;
    prune_limit_ = prune_limit;
}

// O(lgn) sifting up to maintain heap property after insertion (should sift down when building heap)
void BinaryHeap::sift_up(size_t leaf, heap_elem const& node) {
    while ((leaf > 1) && (node.cost < heap_[parent(leaf)].cost)) {
        // sift hole up
        heap_[leaf] = heap_[parent(leaf)];
        leaf = parent(leaf);
    }
    heap_[leaf] = node;
}

// expands heap by "realloc"
void BinaryHeap::expand_heap_if_full() {
    if (heap_tail_ >= heap_size_) { /* Heap is full */
        heap_size_ *= 2;
        heap_.resize(heap_size_ + 1);
    }
}

// adds an element to the back of heap and expand if necessary, but does not maintain heap property
void BinaryHeap::push_back(t_heap* const hptr) {
    expand_heap_if_full();
    heap_elem new_elem = {hptr, hptr->cost};
    heap_[heap_tail_] = new_elem;
    ++heap_tail_;

    check_prune_limit();
}

bool BinaryHeap::is_valid() const {
    if (heap_.empty()) {
        return false;
    }

    for (size_t i = 1; i <= parent(heap_tail_); ++i) {
        size_t leftmost_child = first_child(i);

        for (size_t j = 0; j < 4; ++j) {
            if (leftmost_child + j >= heap_tail_)
                break;
            else if (heap_[leftmost_child + j].cost < heap_[i].cost)
                return false;
        }
    }

    return true;
}

void BinaryHeap::free_all_memory() {
    if (!heap_.empty()) {
        empty_heap();

        // coverity[offset_free : Intentional]
        heap_.clear();
    }

    //  heap_ = nullptr; /* Defensive coding:  crash hard if I use these. */

    storage_.free_all_memory();
}

bool BinaryHeap::check_prune_limit() {
    if (heap_tail_ > prune_limit_) {
        prune_heap();
        return true;
    }

    return false;
}

void BinaryHeap::prune_heap() {
    VTR_ASSERT(max_index_ < prune_limit_);

    std::vector<heap_elem> best_heap_item(max_index_, {nullptr, 0.0});

    // Find the cheapest instance of each index and store it.
    for (size_t i = 1; i < heap_tail_; i++) {
        if (heap_[i].elem_ptr == nullptr) {
            continue;
        }

        if (!heap_[i].elem_ptr->index.is_valid()) {
            free(heap_[i].elem_ptr);
            heap_[i].elem_ptr = nullptr;
            continue;
        }

        auto idx = size_t(heap_[i].elem_ptr->index);

        VTR_ASSERT(idx < max_index_);

        if (best_heap_item[idx].elem_ptr == nullptr || best_heap_item[idx].cost > heap_[i].cost) {
            best_heap_item[idx] = heap_[i];
        }
    }

    // Free unused nodes.
    for (size_t i = 1; i < heap_tail_; i++) {
        if (heap_[i].elem_ptr == nullptr) {
            continue;
        }

        auto idx = size_t(heap_[i].elem_ptr->index);

        if (best_heap_item[idx].elem_ptr != heap_[i].elem_ptr) {
            free(heap_[i].elem_ptr);
            heap_[i].elem_ptr = nullptr;
        }
    }

    heap_tail_ = 1;

    for (size_t i = 0; i < max_index_; ++i) {
        if (best_heap_item[i].elem_ptr != nullptr) {
            heap_[heap_tail_++] = best_heap_item[i];
        }
    }
}

#elif defined(VPR_ROUTER_USE_STL_BINARY_HEAP_WITH_PRIORITY_QUEUE_WRAPPER)

BinaryHeap::BinaryHeap() {
    VTR_LOG("VPR router use STL binary heap with priority queue wrapper\n");
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}

void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    heap_.clear();
    // `heap_` is used as underlying container to build STL priority queue (binary heap)
    heap_.reserve((grid.width() - 1) * (grid.height() - 1));
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    pq_.push(hptr);
}

// Add an element to the back of heap and expand if necessary, but do not maintain heap property
void BinaryHeap::push_back(t_heap* const hptr) {
    heap_.push_back(hptr);
}

bool BinaryHeap::is_empty_heap() const {
    return pq_.empty();
}

bool BinaryHeap::is_valid() const {
    return true;
}

void BinaryHeap::empty_heap() {
    while (!pq_.empty()) {
        free(pq_.top());
        pq_.pop();
    }
}

t_heap* BinaryHeap::get_heap_head() {
    heap_item_t top_item = pq_.top();
    pq_.pop();
    return top_item;
}

void BinaryHeap::build_heap() {
    pq_ = std::move(std::priority_queue<heap_item_t, std::vector<heap_item_t>, heap_cmp>(
        heap_cmp(), std::move(heap_)));
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    (void)max_index;
    (void)prune_limit;
}

void BinaryHeap::free_all_memory() {
    storage_.free_all_memory();
}

#else // D-Ary Heaps, STL Binary Heap

#    include "d_ary_heap.hpp"

#    if defined(VPR_ROUTER_USE_GITHUB_TWO_ARY_HEAP)
constexpr int D_ARITY = 2;
#    elif defined(VPR_ROUTER_USE_GITHUB_FOUR_ARY_HEAP)
constexpr int D_ARITY = 4;
#    elif defined(VPR_ROUTER_USE_GITHUB_EIGHT_ARY_HEAP)
constexpr int D_ARITY = 8;
#    elif defined(VPR_ROUTER_USE_GITHUB_SIXTEEN_ARY_HEAP)
constexpr int D_ARITY = 16;
#    elif defined(VPR_ROUTER_USE_STL_BINARY_HEAP)
constexpr int D_ARITY = -1;
#    endif

BinaryHeap::BinaryHeap() {
    if constexpr (D_ARITY == -1) {
        VTR_LOG("VPR router use STL binary heap\n");
    } else {
        VTR_LOG("VPR router use GitHub %d-ary heap\n", D_ARITY);
    }
}

BinaryHeap::~BinaryHeap() {
    free_all_memory();
}

t_heap* BinaryHeap::alloc() {
    return storage_.alloc();
}
void BinaryHeap::free(t_heap* hptr) {
    storage_.free(hptr);
}

void BinaryHeap::init_heap(const DeviceGrid& grid) {
    heap_.clear();
    heap_.reserve((grid.width() - 1) * (grid.height() - 1));
}

void BinaryHeap::add_to_heap(t_heap* hptr) {
    heap_.push_back(hptr);
    if constexpr (D_ARITY == -1) {
        std::push_heap(heap_.begin(), heap_.end(), heap_cmp());
    } else {
        push_dary_heap<D_ARITY>(heap_.begin(), heap_.end(), heap_cmp());
    }
}

bool BinaryHeap::is_empty_heap() const {
    return (bool)(heap_.empty());
}

t_heap* BinaryHeap::get_heap_head() {
    t_heap* cheapest = heap_.front();
    if constexpr (D_ARITY == -1) {
        std::pop_heap(heap_.begin(), heap_.end(), heap_cmp());
    } else {
        pop_dary_heap<D_ARITY>(heap_.begin(), heap_.end(), heap_cmp());
    }
    heap_.pop_back();
    return (cheapest);
}

void BinaryHeap::empty_heap() {
    size_t heap_size_ = heap_.size();
    for (size_t i = 0; i < heap_size_; i++) {
        free(heap_[i]);
    }
    heap_.clear();
}

void BinaryHeap::build_heap() {
    if constexpr (D_ARITY == -1) {
        std::make_heap(heap_.begin(), heap_.end(), heap_cmp());
    } else {
        make_dary_heap<D_ARITY>(heap_.begin(), heap_.end(), heap_cmp());
    }
}

void BinaryHeap::set_prune_limit(size_t max_index, size_t prune_limit) {
    (void)max_index;
    (void)prune_limit;
}

void BinaryHeap::push_back(t_heap* const hptr) {
    heap_.push_back(hptr);
}

bool BinaryHeap::is_valid() const {
    return true;
}

void BinaryHeap::free_all_memory() {
    if (!heap_.empty()) {
        empty_heap();
    }
    storage_.free_all_memory();
}

#endif
