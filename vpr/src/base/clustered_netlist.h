#pragma once
/**
 * @file
 * @brief This file defines the ClusteredNetlist class in the ClusteredContext
 *        created during pre-placement stages of the VTR flow (packing & clustering),
 *        and used downstream.
 *
 * Overview
 * ========
 * The ClusteredNetlist is derived from the Netlist class, and contains some
 * separate information on Blocks, Pins, and Nets. It does not make use of Ports.
 *
 * Blocks
 * ------
 * The pieces of unique block information are:
 *       block_pbs_:         Physical block describing the clustering and internal hierarchy
 *                           structure of each CLB.
 *       block_types_:       The type of physical block the block is mapped to, e.g. logic
 *                           block, RAM, DSP (Can be user-defined types).
 *       block_nets_:        Based on the block's pins (indexed from [0...num_pins - 1]),
 *                           lists which pins are used/unused with the net using it.
 *       block_pin_nets_:    Returns the index of a pin relative to the net, when given a block and a pin's
 *                           index on that block (from the type descriptor).
 *                           Differs from block_nets_.
 *
 * Differences between block_nets_ & block_pin_nets_
 * --------------------------------------------------
 * ```
 *           +-----------+
 *       0-->|           |-->3
 *       1-->|   Block   |-->4
 *       2-->|           |-->5
 *           +-----------+
 * ```
 * <!-- NOTE: Markdown ``` are used to display the code properly in the documentation -->
 *
 * block_nets_ tracks all pins on a block, and returns the ClusterNetId to which a pin is connected to.
 * If the pin is unused/open, ClusterNetId::INVALID() is stored.
 *
 * block_pin_nets_ tracks whether the nets connected to the block are drivers/receivers of that net.
 * Driver/receiver nets are determined by the pin_class of the block's pin.
 * A net connected to a driver pin in the block has a 0 is stored. A net connected to a receiver
 * has a counter (from [1...num_sinks - 1]).
 *
 * The net is connected to multiple blocks. Each block_pin_nets_ has a unique number in that net.
 *
 * E.g.
 * ```
 *           +-----------+                   +-----------+
 *       0-->|           |-->3   Net A   0-->|           |-->3
 *       1-->|  Block 1  |---4---------->1-->|  Block 2  |-->4
 *       2-->|           |-->5           2-->|           |-->5
 *           +-----------+               |   +-----------+
 *                                       |
 *                                       |   +-----------+
 *                                       |   |           |-->1
 *                                       0-->|  Block 3  |
 *                                           |           |-->2
 *                                           +-----------+
 * ```
 *
 * In the example, Net A is driven by Block 1, and received by Blocks 2 & 3.
 * For Block 1, block_pin_nets_ of pin 4 returns 0, as it is the driver.
 * For Block 2, block_pin_nets_ of pin 1 returns 1 (or 2), non-zero as it is a receiver.
 * For Block 3, block_pin_nets_ of pin 0 returns 2 (or 1), non-zero as it is also a receiver.
 *
 * The block_pin_nets_ data structure exists for quick indexing, rather than using a linear search
 * with the available functions from the base Netlist, into the net_delay_ structure in the
 * PostClusterDelayCalculator of inter_cluster_delay(). net_delay_ is a 2D array, where the indexing
 * scheme is [net_id] followed by [pin_index on net].
 *
 * Pins
 * ----
 * The only piece of unique pin information is:
 *       logical_pin_index_
 *
 * Example of logical_pin_index_
 * ---------------------
 * Given a ClusterPinId, logical_pin_index_ will return the index of the pin within its block
 * relative to the t_logical_block_type (logical description of the block).
 *
 * ```
 *           +-----------+
 *       0-->|O         X|-->3
 *       1-->|O  Block  O|-->4
 *       2-->|X         O|-->5 (e.g. ClusterPinId = 92)
 *           +-----------+
 * ```
 *
 * The index skips over unused pins, e.g. CLB has 6 pins (3 in, 3 out, numbered [0...5]), where
 * the first two ins, and last two outs are used. Indices [0,1] represent the ins, and [4,5]
 * represent the outs. Indices [2,3] are unused. Therefore, logical_pin_index_[92] = 5.
 *
 *
 * Implementation
 * ==============
 * For all create_* functions, the ClusteredNetlist will wrap and call the Netlist's version as it contains
 * additional information that the base Netlist does not know about.
 *
 * All functions with suffix *_impl() follow the Non-Virtual Interface (NVI) idiom.
 * They are called from the base Netlist class to simplify pre/post condition checks and
 * prevent Fragile Base Class (FBC) problems.
 *
 * Refer to netlist.h for more information.
 *
 */
#include "vpr_types.h"

#include "netlist.h"
#include "clustered_netlist_fwd.h"

class ClusteredNetlist : public Netlist<ClusterBlockId, ClusterPortId, ClusterPinId, ClusterNetId> {
  public:
    /**
     * @brief Constructs a netlist
     *
     *   @param name   the name of the netlist (e.g. top-level module)
     *   @param id     a unique identifier for the netlist (e.g. a secure digest of the input file)
     */
    ClusteredNetlist(std::string name = "", std::string id = "");

  public: //Public Accessors
    /*
     * Blocks
     */

    ///@brief Returns the physical block
    t_pb* block_pb(const ClusterBlockId id) const;

    ///@brief Returns the type of CLB (Logic block, RAM, DSP, etc.)
    t_logical_block_type_ptr block_type(const ClusterBlockId id) const;

    ///@brief Returns the net of the block attached to the specific pin index
    ClusterNetId block_net(const ClusterBlockId blk_id, const int pin_index) const;

    ///@brief Returns the count on the net of the block attached
    int block_pin_net_index(const ClusterBlockId blk_id, const int pin_index) const;

    ///@brief Returns the logical pin Id associated with the specified block and logical pin index
    ClusterPinId block_pin(const ClusterBlockId blk, const int logical_pin_index) const;

    ////@brief Returns true if the specified block contains a primary input (e.g. BLIF .input primitive)
    bool block_contains_primary_input(const ClusterBlockId blk, const LogicalModels& models) const;

    ///@brief Returns true if the specified block contains a primary output (e.g. BLIF .output primitive)
    bool block_contains_primary_output(const ClusterBlockId blk, const LogicalModels& models) const;

    /*
     * Pins
     */

    /**
     * @brief Returns the logical pin index (i.e. pin index on the
     *        t_logical_block_type) of the cluster pin
     */
    int pin_logical_index(const ClusterPinId pin_id) const;

    /**
     * @brief Finds the net_index'th net pin (e.g. the 6th pin of the net) and
     *        returns the logical pin index (i.e. pin index on the
     *        t_logical_block_type) of the block to which the pin belongs
     *
     *   @param net_id          The net
     *   @param net_pin_index   The index of the pin in the net
     */
    int net_pin_logical_index(const ClusterNetId net_id, int net_pin_index) const;

    /*
     * Nets
     */

  public: //Public Mutators
    /**
     * @brief Create a new block in the netlist.
     *
     *   @param name   The unique name of the block
     *   @param pb     The physical representation of the block
     *   @param type   The type of the CLB
     */
    ClusterBlockId create_block(const char* name, t_pb* pb, t_logical_block_type_ptr type);

    /**
     * @brief Create a new port in the netlist.
     *
     *   @param blk_id   The block the port is associated with
     *   @param name     The name of the port (must match the name of a port in the block's model)
     *   @param width    The width (number of bits) of the port
     *   @param type     The type of the port (INPUT, OUTPUT, or CLOCK)
     */
    ClusterPortId create_port(const ClusterBlockId blk_id, const std::string& name, BitIndex width, PortType type);
    /**
     * @brief Create a new pin in the netlist.
     * @note If a pin with the specified ID already exists, the function will crash.
     *
     *   @param port_id    The port this pin is associated with
     *   @param port_bit   The bit index of the pin in the port
     *   @param net_id     The net the pin drives/sinks
     *   @param pin_type   The type of the pin (driver/sink)
     *   @param pin_index  The index of the pin relative to its block, excluding OPEN pins)
     *   @param is_const   Indicates whether the pin holds a constant value (e. g. vcc/gnd)
     */
    ClusterPinId create_pin(const ClusterPortId port_id, BitIndex port_bit, const ClusterNetId net_id, const PinType pin_type, int pin_index, bool is_const = false);

    /**
     * @brief Create a net in the netlist
     *
     *   @param name  The unique name of the net
     */
    ClusterNetId create_net(const std::string& name);

    /**
     * @brief Given a name of a block and vector of possible cluster blocks
     *        that are candidates to match the block name, go through 
     *        the vector of cluster blocks and return the id of the block
     *        where the block name matches the provided name.
     * 
     *        Given a string pattern representing a block name and a vector of
     *        poissble cluster blocks that are candidates to match to the block
     *        name pattern, go through the vector of cluster blocks and return 
     *        the id of the block where the block name matches to the provided
     *        input pattern.
     * 
     */

    /**
     * @brief The intented use is to find the block id of a 
     *        hard block without knowing its name in the netlist. Instead
     *        a pattern can be created that we know the block's name will
     *        match to. Generally, we expect the pattern to be constructed
     *        using the block's module name in the HDL design, since we can
     *        expect the netlist name of the block to include the module
     *        name within it.
     * 
     *        For example, suppose a RAM block was named in the netlist as
     *        "top|alu|test_ram|out". The user instantiated the ram module
     *        in the HDL design as "test_ram". So instead of going through 
     *        the netlist and finding the ram block's full name, this
     *        function can be used by just providing the string pattern as
     *        ".*test_ram.*". We know that the module name should be somewhere
     *        within the string, so the pattern we provide says that the netlist
     *        name of the block contains arbritary characters then the module
     *        name and then some other arbritary characters after.
     *        This pattern will then be used to match to the block in the
     *        netlist. The matched cluster block id is returned, and if no
     *        block matched to the input string then an invalid block id
     *        is returned.
     * 
     *        The function 
     *        here additionally requires a vector of possible cluster blocks
     *        that can match to the input pattern. This vector can be the
     *        entire netlist or a subset. For example, if the intended use 
     *        is to find hard blocks, it is quite inefficient
     *        to go through all the blocks to find a matching one. Instead, a
     *        a filtered vector can be passed that is a subset of the entire 
     *        netlist and can only contain blocks of a certain logical block 
     *        type (blocks that can be placed on a specific hard block).
     *        The idea here is that the filtered vector should be
     *        considereably smaller that a vector of every block in the netlist
     *        and this should help improve run time.
     * 
     *        This function runs in linear time (O(N)) as it goes over the
     *        vector of 'cluster_block_candidates'. 'cluster_block_candidates'
     *        could be the whole netlist or a subset as explained above.
     *        Additionally, if there are multiple
     *        blocks that match to the provided input pattern, then the
     *        first block found is returned.
     * 
     * 
     * @param name_pattern A regex string pattern that can be used to match to  
     *             a clustered block name within the netlist.
     * @param cluster_block_candidates A vector of clustered block ids that
     *                                 represent a subset of the clustered
     *                                 netlist. The blocks in this vector
     *                                 will be used to compare and match
     *                                 to the input string pattern. 
     * @return A cluster block id representing a unique cluster block that 
     *         matched to the input string pattern.
     *         
     */
    ClusterBlockId find_block_by_name_fragment(const std::string& name_pattern, const std::vector<ClusterBlockId>& cluster_block_candidates) const;

  private: //Private Members
    /*
     * Netlist compression/optimization
     */

    ///@brief Removes invalid components and reorders them
    void clean_blocks_impl(const vtr::vector_map<ClusterBlockId, ClusterBlockId>& block_id_map) override;
    void clean_ports_impl(const vtr::vector_map<ClusterPortId, ClusterPortId>& port_id_map) override;
    void clean_pins_impl(const vtr::vector_map<ClusterPinId, ClusterPinId>& pin_id_map) override;
    void clean_nets_impl(const vtr::vector_map<ClusterNetId, ClusterNetId>& net_id_map) override;

    ///@brief Shrinks internal data structures to required size to reduce memory consumption
    void shrink_to_fit_impl() override;

    /*
     * Component removal
     */

    /**
     * @brief Removes a block from the netlist.
     *
     * This will also remove the associated ports and pins.
     *   @param blk_id   The block to be removed
     */
    void remove_block_impl(const ClusterBlockId blk_id) override;
    void remove_port_impl(const ClusterPortId port_id) override;
    void remove_pin_impl(const ClusterPinId pin_id) override;
    void remove_net_impl(const ClusterNetId net_id) override;

    void rebuild_block_refs_impl(const vtr::vector_map<ClusterPinId, ClusterPinId>& pin_id_map, const vtr::vector_map<ClusterPortId, ClusterPortId>& port_id_map) override;
    void rebuild_port_refs_impl(const vtr::vector_map<ClusterBlockId, ClusterBlockId>& block_id_map, const vtr::vector_map<ClusterPinId, ClusterPinId>& pin_id_map) override;
    void rebuild_pin_refs_impl(const vtr::vector_map<ClusterPortId, ClusterPortId>& port_id_map, const vtr::vector_map<ClusterNetId, ClusterNetId>& net_id_map) override;
    void rebuild_net_refs_impl(const vtr::vector_map<ClusterPinId, ClusterPinId>& pin_id_map) override;

    /*
     * Sanity Checks
     */

    //Verify the internal data structure sizes match
    bool validate_block_sizes_impl(size_t num_blocks) const override;
    bool validate_port_sizes_impl(size_t num_ports) const override;
    bool validate_pin_sizes_impl(size_t num_pins) const override;
    bool validate_net_sizes_impl(size_t num_nets) const override;

  private: //Private Data
    //Blocks
    vtr::vector_map<ClusterBlockId, t_pb*> block_pbs_;                              ///<Physical block representing the clustering & internal hierarchy of each CLB
    vtr::vector_map<ClusterBlockId, t_logical_block_type_ptr> block_types_;         ///<The type of logical block this user circuit block is mapped to
    vtr::vector_map<ClusterBlockId, std::vector<ClusterPinId>> block_logical_pins_; ///<The logical pin associated with each physical tile pin
    std::unordered_map<int, std::vector<ClusterBlockId>> blocks_per_type_;          ///<Block IDs associated with each physical block type, Used in placement to move specific block type

    //Pins
    /**
     * @brief The logical pin index of this block (i.e. pin index in t_logical_block_type)
     *        corresponding to the clustered pin
     */
    vtr::vector_map<ClusterPinId, int> pin_logical_index_;

    //Nets
};
