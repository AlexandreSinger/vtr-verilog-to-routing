#ifndef _BINARY_HEAP_H
#define _BINARY_HEAP_H

// Note: We use the same name header/class name `Binary Heap` for different
// implementations of the router heaps for convenience. Some implementations
// may not be binary heaps at all.

// #define VPR_ROUTER_USE_CUSTOMIZED_BINARY_HEAP

// #define VPR_ROUTER_USE_GITHUB_TWO_ARY_HEAP
// #define VPR_ROUTER_USE_GITHUB_FOUR_ARY_HEAP

// #define VPR_ROUTER_USE_STL_BINARY_HEAP
#define VPR_ROUTER_USE_STL_BINARY_HEAP_WITH_PRIORITY_QUEUE_WRAPPER


#include "heap_type.h"
#include <vector>
#include <queue>

class BinaryHeap : public HeapInterface {
  public:
    using heap_item_t = t_heap*;
    struct heap_cmp {
        bool operator()(const heap_item_t &u, const heap_item_t &v) {
            return u->cost >= v->cost;
        }
    };

    BinaryHeap();
    ~BinaryHeap();

    t_heap* alloc() final;
    void free(t_heap* hptr) final;

    void init_heap(const DeviceGrid& grid) final;
    void add_to_heap(t_heap* hptr) final;
    void push_back(t_heap* const hptr) final;
    bool is_empty_heap() const final;
    bool is_valid() const final;
    void empty_heap() final;
    t_heap* get_heap_head() final;
    void build_heap() final;
    void set_prune_limit(size_t max_index, size_t prune_limit) final;

    void free_all_memory() final;

  private:
    HeapStorage storage_;
    std::vector<heap_item_t> heap_;

#if defined(VPR_ROUTER_USE_CUSTOMIZED_BINARY_HEAP)

    size_t size() const;
    void sift_up(size_t leaf, t_heap* const node);
    void sift_down(size_t hole);
    void expand_heap_if_full();
    bool check_prune_limit();
    void prune_heap();

    size_t heap_size_;          /* Number of slots in the heap array */
    size_t heap_tail_;          /* Index of first unused slot in the heap array */

    size_t max_index_;
    size_t prune_limit_;

#elif defined(VPR_ROUTER_USE_STL_BINARY_HEAP_WITH_PRIORITY_QUEUE_WRAPPER)

    std::priority_queue<heap_item_t, std::vector<heap_item_t>, heap_cmp> pq_;

#endif

};

#endif /* _BINARY_HEAP_H */
