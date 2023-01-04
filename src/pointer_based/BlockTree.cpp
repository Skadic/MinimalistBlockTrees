#include "pointer_based/BlockTree.h"

#include <algorithm>
#include <cstdint>
#include <iostream>
#include <pointer_based/blocks/InternalBlock.h>
#include <pointer_based/blocks/LeafBlock.h>

#include <stack>
#include <unordered_set>
#include <vector>

using Level        = BlockTree::Level;
using BlockMap     = BlockTree::BlockMap;
using BlockPairMap = BlockTree::BlockPairMap;

BlockTree::BlockTree(std::string &input,
                     int          arity,
                     int          root_block_arity,
                     int          leaf_length,
                     bool         process_block_tree,
                     bool         rs_support) :
    arity_(arity),
    root_arity_(root_block_arity),
    input_(input),
    leaf_length_(leaf_length),
    rank_select_support_(rs_support) {

    // If a leaf is supposed to be greater than the text, or if the arity is greater than the text is long, just make
    // the entire tree one leaf node
    if (input_.size() <= leaf_length_ || input_.size() < arity) {
        root_block_ = new LeafBlock(nullptr, 0, input_.size() - 1, input_);
    } else {
        const int number_of_leaves =
            (input_.size() % leaf_length_ == 0) ? input_.size() / leaf_length_ : input_.size() / leaf_length_ + 1;
        // The current height of the tree
        int height = 0;

        // The current max. index of nodes at the current height
        // So if there are n nodes at the current level, nl is n-1
        int     nl           = number_of_leaves - 1;
        int64_t block_length = leaf_length_;

        // In total this calculates the block length of the first level of blocks
        while (nl > 0) {
            height++;
            block_length *= arity_;
            nl /= arity_;
        }

        root_block_ = new InternalBlock(nullptr, 0, block_length - 1, input_);
    }

    // Actually build the block tree by actually building backlinks etc.
    if (process_block_tree) {
        this->process_block_tree();
        this->clean_unnecessary_expansions();
    }

    // Create rank-select support for all characters if desired
    if (rank_select_support_) {
        std::unordered_set<int> characters;
        for (char c : input) {
            characters.insert(c);
        }
        for (char c : characters) {
            this->add_rank_select_support(c);
        }
    }
}

BlockTree::~BlockTree() { delete root_block_; }

void BlockTree::add_rank_select_support(int c) { root_block_->add_rank_select_support(c); }

int BlockTree::rank(int c, int i) { return root_block_->rank(c, i); }

int BlockTree::select(int c, int j) { return root_block_->select(c, j); }

std::vector<Level> BlockTree::levelwise_iterator() {
    std::vector<Level> result = {{root_block_}};
    while (!dynamic_cast<LeafBlock *>(result.back()[0])) {
        Level next_level = {};
        for (Block *b : result.back()) {
            const auto block_arity = b == root_block_ ? root_arity_ : arity_;
            for (Block *child : b->children(leaf_length_, block_arity)) {
                next_level.push_back(child);
            }
        }
        result.push_back(next_level);
    }

    return result;
}

void BlockTree::clean_unnecessary_expansions() {
    root_block_->clean_unnecessary_expansions();
    for (Level level : levelwise_iterator()) {
        for (int i = 0; i < level.size(); ++i) {
            level[i]->level_index_                  = i;
            level[i]->first_occurrence_level_index_ = level[i]->first_block_->level_index_;
        }
    }
}

int BlockTree::access(int i) { return root_block_->access(i); }

Level BlockTree::next_level(Level &level) {
    // The arity of the current level. If we're at the root, then the arity is of course the root arity.
    const auto block_arity = (level.size() == 1 && level[0] == root_block_) ? root_arity_ : arity_;

    Level next_level;
    // Take all the children of this level's blocks and add them to the new level
    for (int i = 0; i < level.size(); ++i) {
        Block *b = level[i];
        for (Block *child : b->children(leaf_length_, block_arity)) { // Do it in order
            child->level_index_                  = next_level.size();
            child->first_occurrence_level_index_ = next_level.size();
            next_level.push_back(child);
        }
    }
    return next_level;
}

///
/// @brief Returns true if the blocks are adjacent in the text.
///
/// @param left The left block. This must be the block appearing earlier in the text
/// @param left The right block. This must be the block appearing later in the text
///
inline bool blocks_adjacent(const Block *left, const Block *right) {
    return left->end_index_ == right->start_index_ - 1;
}

void BlockTree::window_block_pair_scan(Level        &level,
                                               int           pair_window_size,
                                               int           N,
                                               BlockPairMap &pair_hashtable) {
    // Iterate through all blocks on this level
    // Each iteration of this for loop will start at a segment on this level which consists of blocks that appear
    // contiguously in the input
    for (auto it = level.begin(); it != level.end();) {
        Block *segment_start_block  = *it;
        segment_start_block->right_ = true;
        int offset                  = 0;
        // offset is always 0 here
        RabinKarp rk(input_, segment_start_block->start_index_ + offset, pair_window_size, N);
        // Iterate through the blocks of this contiguous segment
        for (; it != level.end() && (*it == segment_start_block || blocks_adjacent(*(it - 1), *it)); it++) {
            Block *current = *it;
            // This is true, if the current block has no next block or the next block is not adjacent to it
            bool is_segment_end_block = (it + 1) == level.end() || !blocks_adjacent(current, *(it + 1));
            for (offset = 0; offset < current->length(); ++offset) {
                // If this block is the last block of the segment, and the window does not fit in this block anymore,
                // break, as we need to include the next block for this to fit
                if (is_segment_end_block && current->length() - offset < pair_window_size) {
                    break;
                }
                // Hash the window
                HashString hS(rk.hash(),
                              input_,
                              current->start_index_ + offset,
                              current->start_index_ + offset + pair_window_size - 1);
                // Try to find blocks that match the hashed string
                auto result = pair_hashtable.find(hS);
                if (result != pair_hashtable.end()) {
                    // v Idk I didn't write this v
                    // Here, It could be that the scanning should have finished with the penultimate, but it never
                    // should enter this ''if'' when We're on the penultimate block and the window exceeds the last
                    // block because if that is a first occurrence should have been occured before in a pair of blocks
                    // maybe use a condition more like rk's coandition below could work fine too
                    // Same logic: for when passing a window of size 2l + 2 over 2 block of length l
                    // ^ Idk I didn't write this ^
                    for (auto [left_block, right_block] : result->second) {
                        // left_ and right_ should be true if the string described by left_block and right_block was
                        // found somewhere else in the text.
                        // Thus we do not want to set left_ and right_ to true, if the string we hashed before is the
                        // same (not just "equal" but "the very same") string that is described by left_block and
                        // right_block.
                        if (current->start_index_ + offset < left_block->start_index_) {
                            left_block->left_   = true;
                            right_block->right_ = true;
                        }
                    }
                    // Delete the pair from the hashtable
                    pair_hashtable.erase(hS);
                }
                // If we can move the window to the right without exceeding the input, get the next hash
                if (current->start_index_ + offset + pair_window_size < input_.size()) {
                    rk.next();
                }
            }
        }
        (*(it - 1))->left_ = true;
    }
}

void BlockTree::window_block_scan(Level &level, int window_size, int N, BlockMap &hashtable) {
    // The current block index inside this level
    int i = 0;
    // Iterate through all blocks on this level
    // Each iteration of this for loop will start at a segment on this level which consists of blocks that appear
    // contiguously in the input
    for (auto it = level.begin(); it != level.end();) {
        const Block *segment_start_block = *it;
        int          offset              = 0;
        RabinKarp    rk(input_, segment_start_block->start_index_ + offset, window_size, N);
        for (; it != level.end() && (*it == segment_start_block || blocks_adjacent(*(it - 1), *it)); it++, i++) {
            const Block *current = *it;
            // This is true, if the current block has no next block or the next block is not adjacent to it
            const bool is_segment_end_block = (it + 1) == level.end() || !blocks_adjacent(current, *(it + 1));
            // For every possible offset that does not exceed this block's length
            for (offset = 0; offset < current->length(); ++offset) {
                // If this block is the last block of the segment,
                // and the window does not fit in this block anymore, break
                if (is_segment_end_block && current->length() - offset < window_size) {
                    break;
                }
                // Hash the current window
                HashString hS(rk.hash(),
                              input_,
                              current->start_index_ + offset,
                              current->start_index_ + offset + window_size - 1);
                // Try to find blocks that match this hashed window
                auto result = hashtable.find(hS);
                if (result != hashtable.end()) {
                    std::vector<Block *> &blocks = result->second;
                    // For each block that matches the hashed window, set their first occurrence,
                    // the first block in which the source is contained (since it can span 2 blocks, there is a first
                    // and a second block) and the offset in the source blocks at which the string appears
                    for (Block *b : blocks) {
                        b->first_occurrence_level_index_ = i;
                        b->first_block_                  = const_cast<Block *>(current);
                        b->offset_                       = offset;
                        if (offset + window_size > b->first_block_->length()) {
                            // If the source string exceeds the first block, we need to include the succeeding block as
                            // the second block Note, that we previously excluded the possibility that the next block is
                            // not adjacent *and* the source does not fit into the second block
                            b->second_block_ = *(it + 1);
                        } else {
                            // If the source fits into the first block, there is no need for a second block.
                            b->second_block_ = nullptr;
                        }
                    }
                    // Delete the processed blocks from the hashtable so that we don't process them again
                    hashtable.erase(hS);
                }
                // If the next window does not exceed the input size, we calculate the hash for the next string.
                if (current->start_index_ + offset + window_size < input_.size()) {
                    rk.next();
                }
            }
        }
    }
}

void BlockTree::hash_blocks(Level &level, int N, BlockMap &hashtable) {
    // For every block in this level, calculate its fingerprint and insert it into the hashtable
    for (Block *b : level) {
        RabinKarp  rk(input_, b->start_index_, b->length(), N);
        HashString hS(rk.hash(), input_, b->start_index_, b->end_index_);

        auto result = hashtable.find(hS);

        if (result == hashtable.end()) {
            hashtable[hS] = {b};
        } else {
            hashtable[hS].push_back(b);
        }
    }
}

void BlockTree::hash_block_pairs(Level &level, int N, BlockPairMap &pair_hashtable) {
    for (auto it = level.begin(); it != level.end();) {
        // For every contiguous segment of blocks hash two neighboring blocks and add the block pair to the hashtable
        // with the hashed string as their key
        for (++it; it != level.end() && blocks_adjacent(*(it - 1), *it); ++it) {
            Block     *current          = *(it - 1);
            Block     *next             = *it;
            size_t     pair_length      = current->length() + next->length();
            size_t     pair_start_index = current->start_index_;
            size_t     pair_end_index   = current->start_index_ + pair_length - 1;
            RabinKarp  rk(input_, pair_start_index, pair_length, N);
            HashString hS(rk.hash(), input_, pair_start_index, pair_end_index);

            auto result = pair_hashtable.find(hS);

            if (result == pair_hashtable.end()) {
                pair_hashtable[hS] = {{current, next}};
            } else {
                pair_hashtable[hS].push_back({current, next});
            }
        }
    }
}

void BlockTree::create_back_blocks(Level &level) {
    // BackBlock creation
    for (int i = 0; i < level.size(); ++i) {
        Block *b = level[i];
        // If this block's content (or at least part of it) appears somewhere else in the text both as a right block and
        // as a left block, and its first occurrence is before this block, then we create a backlink to its original
        // source
        if (b->left_ && b->right_ && b->first_occurrence_level_index_ < b->level_index_) {
            // This doesn't have the bug of the dangling reference fixed with first_occurrence_level_index, because it
            // shouldn't happen that A block points back to a BackBlock
            BackBlock *bb =
                new BackBlock(b->parent_,
                              b->start_index_,
                              b->end_index_,
                              input_,
                              level[b->first_occurrence_level_index_],
                              (b->second_block_ == nullptr) ? nullptr : level[b->first_occurrence_level_index_ + 1],
                              b->offset_);
            bb->level_index_                  = b->level_index_;
            bb->first_occurrence_level_index_ = b->first_occurrence_level_index_;
            bb->left_                         = true;
            bb->right_                        = true;
            b->parent_->replace_child(b, bb);
            delete b;
            level[i] = bb;
        }
    }
}

void BlockTree::process_level(Level &level) {
    // Large prime
    const int N = 6700417;
    // The length of blocks on this level
    const int block_length = level.front()->length();

    // Handle single blocks
    {
        // Hash the blocks such that a string maps to all blocks whose blocks have the same hash
        BlockMap hashtable;
        hash_blocks(level, N, hashtable);

        // Window block scan
        // Establishes first occurrences of blocks
        window_block_scan(level, block_length, N, hashtable);
    }

    // Handle block pairs if there is more than one block in this level
    if (level.size() > 1) {
        // We take a pair of blocks, hash the string they represent
        // and insert the block pair into the hashtable at the appropriate place
        BlockPairMap pair_hashtable;
        hash_block_pairs(level, N, pair_hashtable);

        // Window Pair of blocks scans
        window_block_pair_scan(level, block_length * 2, N, pair_hashtable);
    }

    // With the data in the blocks populated, create the back blocks
    create_back_blocks(level);
}

void BlockTree::process_block_tree() {
    Level               current_level = {root_block_};
    std::stack<Block *> none_blocks;
    // While there is a new level in the tree to process
    while ((current_level = next_level(current_level)).size() != 0) {
        const auto current_level_block_length = current_level[0]->length();
        // If the nodes at this levels are at leaf-size,
        // we stop, as there are no back-pointers that can exist on this level
        if (current_level_block_length < arity_ || current_level_block_length <= leaf_length_)
            break;
        // If we have at least 1 block to process on this level and the last node
        // exceeds the length of the text (thus the block is not "full") we take such blocks out before processing
        // as they cannot exist anywhere further to the left on this level (all nodes to the left would be full after
        // all)
        while (current_level.size() != 0 && current_level.back()->end_index_ >= input_.size()) {
            none_blocks.push(current_level.back());
            current_level.pop_back();
        }
        // Process the backlinks on this level
        process_level(current_level);
        // Put the blocks back that we took out before
        while (!none_blocks.empty()) {
            current_level.push_back(none_blocks.top());
            none_blocks.pop();
        }
    }
}

void BlockTree::process_level_heuristic(std::vector<Block *> &level) {

    int N            = 6700417; // Large prime
    int level_length = level.front()->length();

    // Block scan
    std::unordered_map<HashString, std::vector<Block *>> hashtable;
    hash_blocks(level, N, hashtable);

    // Window block scan
    // This is almost the same as forward_window_block_scan, as well as the BackBlock creation
    window_block_scan(level, level_length, N, hashtable);

    // BackBlock creation
    for (int i = 0; i < level.size(); ++i) {
        Block *b = level[i];
        if (b->first_occurrence_level_index_ < b->level_index_) {

            BackBlock *bb =
                new BackBlock(b->parent_,
                              b->start_index_,
                              b->end_index_,
                              input_,
                              level[b->first_occurrence_level_index_],
                              (b->second_block_ == nullptr) ? nullptr : level[b->first_occurrence_level_index_ + 1],
                              b->offset_);
            bb->level_index_                  = b->level_index_;
            bb->first_occurrence_level_index_ = b->first_occurrence_level_index_;
            b->parent_->replace_child(b, bb);
            delete b;
            level[i] = bb;
        }
    }
}

void BlockTree::process_back_pointers_heuristic() {
    std::vector<Block *> current_level = {root_block_};
    std::stack<Block *>  none_blocks;
    while ((current_level = next_level(current_level)).size() != 0) {
        if (current_level[0]->length() < arity_ || current_level[0]->length() <= leaf_length_)
            break;
        while (current_level.size() != 0 && current_level.back()->end_index_ >= input_.size()) {
            none_blocks.push(current_level.back());
            current_level.pop_back();
        }
        process_level_heuristic(current_level);
        while (!none_blocks.empty()) {
            current_level.push_back(none_blocks.top());
            none_blocks.pop();
        }
    }
}

void BlockTree::print() {
    auto it = levelwise_iterator();
    for (auto &level : it) {
        for (Block *b : level) {
            std::cout << "[" << b->start_index_ << ", " << b->end_index_ << "], ";
        }
        std::cout << std::endl;
    }
}
