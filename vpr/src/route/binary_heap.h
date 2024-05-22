#ifndef _BINARY_HEAP_H
#define _BINARY_HEAP_H

// Note: We use the same name header/class name `Binary Heap` for different
// implementations of the router heaps for convenience. Some implementations
// may not be binary heaps at all.

// VTR original binary heap
// #define VPR_ROUTER_USE_CUSTOMIZED_BINARY_HEAP

// CUSTOMIZED_BINARY_HEAP, with indexing scheme changed so first index is 0 instead of 1
// #define VPR_ROUTER_USE_CUSTOMIZED_BINARY_HEAP_FIRST_INDEX_ZERO

// Modified VTR original binary heap to make it a 4-ary heap
// #define VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP

// CUSTOMIZED_FOUR_ARY_HEAP, but with indexing scheme changed so first index is 0 instead of 1; crucially, the second
// layer of the heap (starting at index 1) only has 3 nodes, while the rest of the heap has 4 like normal
// #define VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FIRST_INDEX_ZERO

// DO NOT USE. CUSTOMIZED_FOUR_ARY_HEAP_FIRST_INDEX_ZERO, but uses heap_elem, a struct which contains a 32-bit ID and
// the node cost; this cannot be used in a real implementation currently due the use of the arbitrary ID
// #define VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM

// CUSTOMIZED_FOUR_ARY_HEAP_FIRST_INDEX_ZERO, but uses heap_elem, a struct which contains the node pointer and its cost
// #define VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM_LARGE

// CUSTOMIZED_FOUR_ARY_HEAP, but using a heap_elem as described in CUSTOMIZED_FOUR_ARY_HEAP_ELEM_LARGE, with this struct
// using 16-byte alignment
#define VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FINAL

// #define VPR_ROUTER_USE_GITHUB_TWO_ARY_HEAP
// #define VPR_ROUTER_USE_GITHUB_FOUR_ARY_HEAP
// #define VPR_ROUTER_USE_GITHUB_EIGHT_ARY_HEAP
// #define VPR_ROUTER_USE_GITHUB_SIXTEEN_ARY_HEAP

// #define VPR_ROUTER_USE_STL_BINARY_HEAP
// #define VPR_ROUTER_USE_STL_BINARY_HEAP_WITH_PRIORITY_QUEUE_WRAPPER

#include "heap_type.h"
#include <vector>
#include <queue>
#include <boost/align/aligned_allocator.hpp>

class BinaryHeap : public HeapInterface {
  public:
    using heap_item_t = t_heap*;

#if defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM)
    struct heap_elem {
        uint32_t ID;
        float cost;
    } __attribute__((aligned(8)));
#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FIRST_INDEX_ZERO)
    struct heap_elem {
        heap_item_t elem_ptr;
    } __attribute__((aligned(8)));
#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM_LARGE) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FINAL)
    struct heap_elem {
        heap_item_t elem_ptr;
        float cost;
    } __attribute__((aligned(16)));
#endif

    struct heap_cmp {
        bool operator()(const heap_item_t& u, const heap_item_t& v) {
            return u->cost >= v->cost;
        }
    };

    BinaryHeap();
    ~BinaryHeap();

    t_heap* alloc() final;
    void free(t_heap* hptr) final;

    void init_heap(const DeviceGrid& grid) final;

    bool is_empty_heap() const final;
    bool is_valid() const final;
    void empty_heap() final;
    void build_heap() final;
    void free_all_memory() final;

#if defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM)
    void add_to_heap_non_virtual(heap_elem elem);
    void push_back_non_virtual(heap_elem const& elem);
    uint32_t get_heap_head_non_virtual();
#endif

    void set_prune_limit(size_t max_index, size_t prune_limit) final;
    void add_to_heap(t_heap* hptr) final;
    void push_back(t_heap* const hptr) final;
    t_heap* get_heap_head() final;

  private:
    HeapStorage storage_;

#if defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FIRST_INDEX_ZERO)
    std::vector<heap_elem, boost::alignment::aligned_allocator<heap_elem, 8>> heap_;
#elif defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM_LARGE) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FINAL)
    std::vector<heap_elem, boost::alignment::aligned_allocator<heap_elem, 16>> heap_;
#else
    std::vector<heap_item_t> heap_;
#endif

#if defined(VPR_ROUTER_USE_CUSTOMIZED_BINARY_HEAP) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP) || defined(VPR_ROUTER_USE_CUSTOMIZED_BINARY_HEAP_FIRST_INDEX_ZERO) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FIRST_INDEX_ZERO) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM_LARGE) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FINAL)

    size_t size() const;
    void sift_down(size_t hole);
    void expand_heap_if_full();

#    if defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM_LARGE) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FINAL)
    void sift_up(size_t leaf, heap_elem const& node);
#    else
    void sift_up(size_t leaf, t_heap* const node);
#    endif

#    if !defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM)
    bool check_prune_limit();
    void prune_heap();
#    endif

#    if defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FIRST_INDEX_ZERO) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_ELEM_LARGE) || defined(VPR_ROUTER_USE_CUSTOMIZED_FOUR_ARY_HEAP_FINAL)
    size_t smallest_child(size_t i);
#    endif

    size_t heap_size_; /* Number of slots in the heap array */
    size_t heap_tail_; /* Index of first unused slot in the heap array */

    size_t max_index_;
    size_t prune_limit_;

#elif defined(VPR_ROUTER_USE_STL_BINARY_HEAP_WITH_PRIORITY_QUEUE_WRAPPER)

    std::priority_queue<heap_item_t, std::vector<heap_item_t>, heap_cmp> pq_;

#endif
};

#endif /* _BINARY_HEAP_H */
