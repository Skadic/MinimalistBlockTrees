#ifndef BLOCKTREE_PBLOCKTREE_H
#define BLOCKTREE_PBLOCKTREE_H

#include <array>
#include <blocktree.utils/HashString.h>
#include <blocktree.utils/RabinKarp.h>
#include <pointer_based/blocks/BackBlock.h>
#include <pointer_based/blocks/Block.h>

#include <string>
#include <unordered_map>

class BlockTree {

  public:
    using Level        = std::vector<Block *>;
    using BlockMap     = std::unordered_map<HashString, std::vector<Block *>>;
    using BlockPairMap = std::unordered_map<HashString, std::vector<std::pair<Block *, Block *>>>;

  private:
    ///
    /// @brief Scans through this level's blocks, for each block hashes its represented string and inserts it into the
    /// hashtable, mapping from the hashed string to all pairs of blocks that have the same hash.
    ///
    /// The hashes are Rabin-Karp fingerprints of the section of the input that the block spans.
    ///
    /// @param level A list of block pointers that represents a level in the tree
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-hashed string to all blocks that match its hash.
    ///
    void hash_blocks(Level &level, int N, BlockMap &hashtable);

    ///
    /// @brief Scans through this level with a window of two blocks, hashes the represented string for each block pair
    /// and inserts it into a hashtable, mapping from the hashed string to all pairs of blocks that have the same hash.
    ///
    /// The hashes are Rabin-Karp fingerprints of the section of the input that the block pair spans.
    ///
    /// @param level A list of block pointers that represents a level in the tree
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-hashed string to all block pairs that match its hash.
    ///
    void hash_block_pairs(Level &level, int N, BlockPairMap &hashtable);

  public:
    /// This tree's arity (except for the root node)
    int arity_;
    /// The root node's arity
    int root_arity_;
    /// The maximum length of a leaf node.
    int leaf_length_;
    /// Input sequence of the Tree
    std::string input_;
    /// The root block
    Block *root_block_;
    /// Whether this tree is built with rank and select support
    bool rank_select_support_;

    BlockTree(std::string &source,
              int          arity,
              int          root_block_arity,
              int          leaf_length,
              bool         process_block_tree  = false,
              bool         rank_select_support = false);
    ~BlockTree();

    int  access(int);
    void add_rank_select_support(int);
    int  rank(int, int);
    int  select(int, int);

    void process_back_pointers_heuristic();

    ///
    /// @brief Actually build the block tree and insert back-links for repeated blocks.
    ///
    void process_block_tree();
    void clean_unnecessary_expansions();

  private:
    void process_level_heuristic(Level &level);

    ///
    /// @brief Calculate the necessary data, populate the fields in the blocks and create backlinks in this level.
    ///
    /// @param level The level to process.
    ///
    void process_level(Level &level);

    ///
    /// @brief Scan for occurrences of strings (of size window_size) that also exist in the hashtable and set the first
    /// occurrences and backlinks of the blocks that match these strings.
    ///
    /// The hashtable should contain a mapping for each block that exists in the tree. Each mapping should map from the
    /// string represented by the block, to a vector of all blocks that share the same hash value.
    /// This can be calculated with `hash_blocks`.
    ///
    /// This scans through all blocks on this level with a window of the given size and upon recognizing a string that
    /// exists in the hashtable, gets all blocks that match that string and for each block creates a backlink to its
    /// source block(s).
    ///
    /// @param level A list of block pointers that represents a level in the tree.
    /// @param window_size The window size used for scanning.
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-Hashed string to all blocks that match its hash.
    ///
    void window_block_scan(Level &level, int window_size, int N, BlockMap &hashtable);

    ///
    /// @brief Scan for occurrences of strings (of size window_size) that also exist in the hashtable
    /// and for all pairs of blocks that match these strings, they are marked.
    ///
    /// The hashtable should contain mappings for each contiguous block pair that exists in the tree.
    /// Each mapping should map from the string represented by the block pair, to a vector of all block pairs that share
    /// the same hash value.
    /// This can be calculated with `hash_block_pairs`.
    ///
    /// This scans through all blocks on this level with a window of the given size and upon recognizing a string that
    /// exists in the hashtable, gets all block pairs that match that string and for each of these pairs marks in both
    /// blocks that the string represented by them both appeared somewhere in the text.
    ///
    /// @param level A list of block pointers that represents a level in the tree.
    /// @param window_size The window size used for scanning.
    /// @param N A large prime number
    /// @param hashtable A hashtable mapping a Rabin-Karp-Hashed string to all pairs of blocks that match its hash.
    ///
    void window_block_pair_scan(Level &level, int pair_window_size, int N, BlockPairMap &pair_hashtable);

    ///
    /// @brief Create back blocks from the information saved in the blocks on the given level.
    ///
    /// This requires the data in the blocks to be populated. See `window_block_scan` and `window_block_pair_scan` for
    /// more.
    ///
    /// @param level A list of block pointers that represents a level in the tree.
    ///
    void create_back_blocks(Level &level);

  public:
    ///
    /// @brief This takes a vector of blocks and returns a vector of their children in order.
    ///
    /// As this method's name says, this is supposed to be used with input vectors
    /// which consist of blocks of a single level in the tree from left to right.
    ///
    /// @param level A vector of pointers to blocks whose children to return.
    ///
    Level next_level(Level &level);

    /// @brief Returns a vector of levels of nodes of the tree where each level is represented by a vector of its nodes
    /// (left-to-right).
    ///
    /// A simple levelwise (left-to-right) traversal of the tree would be:
    ///     for (std::vector<Block*> level : bt->levelwise_iterator()) {
    ///         for (Block* b : level) {
    ///             ...
    std::vector<Level> levelwise_iterator();

    void print();
};

#endif // BLOCKTREE_PBLOCKTREE_H
