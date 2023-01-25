#include "gtest/gtest.h"
#include <compressed/CBlockTree.h>
#include <pointer_based/blocks/LeafBlock.h>
#include <unordered_set>

#include "pointer_based/BlockTree.h"

using ::testing::Combine;
using ::testing::Values;

using Arity            = int;
using RootArity        = int;
using MaxLeafLength    = int;
using PrefixSuffixSize = int;
using InputText        = std::string;

typedef BlockTree *CreateBlockTreeFunc(Arity, RootArity, MaxLeafLength, PrefixSuffixSize, const InputText &);

using TestParameters =
    ::testing::tuple<Arity, RootArity, MaxLeafLength, PrefixSuffixSize, InputText, CreateBlockTreeFunc *>;

BlockTree *
block_tree(int arity, int root_arity, int max_leaf_length, int prefix_suffix_size, const std::string &input) {
    auto *block_tree_ = new BlockTree(input, arity, root_arity, max_leaf_length);
    block_tree_->process_block_tree();
    block_tree_->clean_unnecessary_expansions();
    block_tree_->add_fast_substring_support(prefix_suffix_size);
    return block_tree_;
}

BlockTree *block_tree_without_cleaning(int                arity,
                                       int                root_arity,
                                       int                max_leaf_length,
                                       int                prefix_suffix_size,
                                       const std::string &input) {
    auto *block_tree_ = new BlockTree(input, arity, root_arity, max_leaf_length);
    block_tree_->process_block_tree();
    block_tree_->add_fast_substring_support(prefix_suffix_size);
    return block_tree_;
}

BlockTree *
heuristic_block_tree(int arity, int root_arity, int max_leaf_length, int prefix_suffix_size, const std::string &input) {
    auto *block_tree_ = new BlockTree(input, arity, root_arity, max_leaf_length);
    block_tree_->process_back_pointers_heuristic();
    block_tree_->add_fast_substring_support(prefix_suffix_size);
    return block_tree_;
}

class CBlockTreeFixture : public ::testing::TestWithParam<TestParameters> {

  public:
    BlockTree *block_tree_;
    BlockTree *block_tree_rs_;

    CBlockTree *c_block_tree_;
    CBlockTree *c_block_tree_rs_;

    std::string input_;

    int arity_;
    int root_arity_;
    int max_leaf_length_;
    int prefix_suffix_size_;

    /// Characters in the input and its select results
    std::unordered_map<int, std::vector<int>> characters_;

    CBlockTreeFixture() :
        ::testing::TestWithParam<TestParameters>(),
        block_tree_(nullptr),
        block_tree_rs_(nullptr),
        c_block_tree_(nullptr),
        c_block_tree_rs_(nullptr),
        input_(),
        arity_(0),
        root_arity_(0),
        max_leaf_length_(0),
        prefix_suffix_size_(0) {}

    ~CBlockTreeFixture() override = default;

  protected:
    void TearDown() override {
        delete block_tree_;
        delete block_tree_rs_;
        delete c_block_tree_;
        delete c_block_tree_rs_;
    }

    void SetUp() override {
        arity_              = ::testing::get<0>(GetParam());
        root_arity_         = ::testing::get<1>(GetParam());
        max_leaf_length_    = ::testing::get<2>(GetParam());
        prefix_suffix_size_ = ::testing::get<3>(GetParam());
        std::ifstream        t(::testing::get<4>(GetParam()));
        CreateBlockTreeFunc *create_blocktree = ::testing::get<5>(GetParam());

        std::stringstream buffer;
        buffer << t.rdbuf();
        input_        = buffer.str();
        block_tree_   = (*create_blocktree)(arity_, root_arity_, max_leaf_length_, prefix_suffix_size_, input_);
        c_block_tree_ = new CBlockTree(block_tree_);

        block_tree_rs_ = (*create_blocktree)(arity_, root_arity_, max_leaf_length_, prefix_suffix_size_, input_);

        std::unordered_set<int> characters;
        for (char c : input_) characters.insert(c);
        for (int c : characters) {
            characters_[c] = {};
            block_tree_rs_->add_rank_select_support(c);
        }

        for (int i = 0; i < input_.size(); ++i) characters_[input_[i]].push_back(i);

        c_block_tree_rs_ = new CBlockTree(block_tree_rs_);
    }
};

INSTANTIATE_TEST_SUITE_P(PCBlockTreeTest,
                         CBlockTreeFixture,
                         Combine(Values(2),
                                 Values(8),
                                 Values(4),
                                 Values(16),
                                 Values("../../../tests/data/as",
                                        "../../../tests/data/dna",
                                        "../../../tests/data/dna.par",
                                        "../../../tests/data/einstein"),
                                 Values(&block_tree, &block_tree_without_cleaning, &heuristic_block_tree)));

// This test checks if the fields and
// number_of_levels_ are correct
TEST_P(CBlockTreeFixture, general_fields_check) {
    EXPECT_EQ(c_block_tree_->arity_, arity_) << "Incorrect tree arity";
    EXPECT_EQ(c_block_tree_->root_arity_, root_arity_) << "Incorrect root arity";
    auto                 iterator = block_tree_->levelwise_iterator();
    std::vector<Block *> level;
    bool                 contains_back_block = false;
    int                  i;
    for (i = 0; i < iterator.size(); ++i) {
        level = iterator[i];
        for (Block *b : level) {
            if (dynamic_cast<BackBlock *>(b))
                contains_back_block = true;
        }
        if (contains_back_block)
            break;
    }

    EXPECT_EQ(iterator.size() - i, c_block_tree_->number_of_levels_) << "Incorrect number of levels";
}

// This test checks if the CBlockTree has the same structure that its correspondent BlockTree in particular the
// is_internal_ field is checked
TEST_P(CBlockTreeFixture, is_internal_check) {
    auto iterator   = block_tree_->levelwise_iterator();
    auto [level, i] = block_tree_->get_lowest_complete_level();

    for (int j = 0; j < c_block_tree_->number_of_levels_ - 1; ++j) {
        level            = iterator[i + j];
        auto is_internal = *(c_block_tree_->is_internal_[j]);
        EXPECT_EQ(level.size(), is_internal.size());
        for (int k = 0; k < level.size(); ++k) {
            Block *b = level[k];
            if (dynamic_cast<BackBlock *>(b))
                EXPECT_FALSE(is_internal[k]);
            else
                EXPECT_TRUE(is_internal[k]);
        }
    }
}

// This test checks if the CBlockTree has the same structure that its correspondent BlockTree in particular the offsets_
// field is checked
TEST_P(CBlockTreeFixture, offsets_check) {
    auto iterator   = block_tree_->levelwise_iterator();
    auto [level, i] = block_tree_->get_lowest_complete_level();

    for (int j = 0; j < c_block_tree_->number_of_levels_ - 1; ++j) {
        level        = iterator[i + j];
        auto offsets = *(c_block_tree_->offsets_[j]);

        int max_size_level = level.front()->length();
        int l              = 0;
        for (Block *b : level) {
            if (dynamic_cast<BackBlock *>(b)) {
                EXPECT_EQ(offsets[l], max_size_level * b->first_block_->level_index_ + b->offset_);
                ++l;
            }
        }
        EXPECT_EQ(l, offsets.size());
    }
}

// This test checks if the CBlockTree has the same
// structure that its correspondent BlockTree
// in particular the leaf_string_ field is checked
TEST_P(CBlockTreeFixture, leaf_string_check) {
    auto        iterator = block_tree_->levelwise_iterator();
    std::string leaf_string;
    for (Block *b : iterator.back()) {
        leaf_string += b->represented_string();
    }

    std::string leaf_c_string;
    for (int i : (*c_block_tree_->leaf_string_)) {
        leaf_c_string += (char) (*c_block_tree_->alphabet_)[i];
    }

    EXPECT_EQ(leaf_c_string, leaf_string);
}

TEST_P(CBlockTreeFixture, prefix_suffix_symbols_check) {
    ASSERT_EQ(prefix_suffix_size_, c_block_tree_->prefix_suffix_size_) << "prefix_suffix_size_ field incorrect";
}

// This test checks if the CBlockTree has the same
// structure that its correspondent BlockTree
// in particular the first_block_pop_counts_ field are checked
TEST_P(CBlockTreeFixture, first_block_pop_counts_check) {
    auto iterator   = block_tree_rs_->levelwise_iterator();
    auto [level, i] = block_tree_->get_lowest_complete_level();

    for (int j = 0; j < c_block_tree_->number_of_levels_ - 1; ++j) {
        level = iterator[i + j];
        for (const auto &[c, _] : characters_) {
            auto level_first_block_pop_counts = *(c_block_tree_rs_->first_block_pop_counts_[c][j]);

            int l = 0;
            for (Block *b : level) {
                if (dynamic_cast<BackBlock *>(b)) {
                    EXPECT_EQ(level_first_block_pop_counts[l], b->pop_counts_in_first_block_[c]);
                    ++l;
                }
            }
            EXPECT_EQ(l, level_first_block_pop_counts.size());
        }
    }
}

// This test checks if the CBlockTree has the same
// structure that its correspondent BlockTree
// in particular the pop_counts_ is checked
TEST_P(CBlockTreeFixture, pop_counts_check) {
    auto iterator   = block_tree_rs_->levelwise_iterator();
    auto [level, i] = block_tree_->get_lowest_complete_level();

    for (const auto &[c, _] : characters_) {
        level                 = iterator[i];
        auto level_pop_counts = *(c_block_tree_rs_->pop_counts_[c][0]);
        EXPECT_EQ(level.size(), level_pop_counts.size());

        for (int k = 0; k < level.size(); ++k) {
            Block *b = level[k];
            EXPECT_EQ(b->pop_counts_[c], level_pop_counts[k]);
        }
    }

    for (int j = 1; j < c_block_tree_->number_of_levels_; ++j) {
        for (const auto &[c, _] : characters_) {
            level               = iterator[i + j];
            auto level_bt_ranks = *(c_block_tree_rs_->pop_counts_[c][j]);
            EXPECT_EQ(level.size(), level_bt_ranks.size());

            for (int k = 0; k < level.size(); ++k) {
                Block *b = level[k];
                EXPECT_EQ(b->pop_counts_[c], level_bt_ranks[k]);
            }
        }
    }
}

// This test checks if the CBlockTree has the same structure that its correspondent BlockTree in particular the first
// level for first_level_ranks_, is checked
TEST_P(CBlockTreeFixture, first_level_ranks_check) {
    auto iterator   = block_tree_rs_->levelwise_iterator();
    auto [level, i] = block_tree_->get_lowest_complete_level();

    for (const auto &[c, _] : characters_) {
        level                  = iterator[i];
        auto first_level_ranks = *(c_block_tree_rs_->first_level_ranks_[c]);
        int  r                 = 0;

        EXPECT_EQ(first_level_ranks.size(), level.size());
        for (int k = 0; k < level.size(); ++k) {
            EXPECT_EQ(r, first_level_ranks[k]);
            r += level[k]->pop_counts_[c];
        }
    }
}

// Check if the saved prefix/suffixes are correct
TEST_P(CBlockTreeFixture, prefix_suffix_check) {
    const auto [_, lowest] = block_tree_->get_lowest_complete_level();
    auto  levels           = block_tree_rs_->levelwise_iterator();
    auto &alphabet         = *c_block_tree_->alphabet_;

    for (size_t level_index = lowest; level_index < levels.size(); level_index++) {
        const BlockTree::Level   &level      = levels[level_index];
        const size_t              block_size = level[0]->length();
        const size_t              level_size = level.size();
        const size_t              saved_data = std::min((int) block_size, prefix_suffix_size_ << 1);
        const sdsl::int_vector<> &symbols    = *c_block_tree_->prefix_suffix_symbols_[level_index - lowest];

        for (size_t block_index = 0; block_index < level_size; block_index++) {
            Block *block  = level[block_index];
            auto  &prefix = block->prefix_;
            auto  &suffix = block->suffix_;

            // Read the prefix from the compressed block tree
            char prefix_buf[prefix.size() + 1];
            prefix_buf[prefix.size()] = 0;
            for (int i = 0; i < prefix.size(); ++i) {
                const auto unmapped_symbol = symbols[block_index * saved_data + i];
                const char mapped_symbol   = alphabet[unmapped_symbol];
                prefix_buf[i]              = mapped_symbol;
            }
            ASSERT_EQ(prefix, std::string_view(prefix_buf, prefix.size()))
                << "Incorrect prefix in block (" << block->start_index_ << ", " << block->end_index_ << ")";

            // In the former case prefix and suffix are identical
            // In the latter case, we don't care about the suffix
            if (saved_data == block_size || block_index == level_size - 1) {
                continue;
            }

            // If the suffix started beyond the end of the input string we don't expect anything worthwhile to be in
            // there anyway
            if (block->end_index_ - prefix_suffix_size_ + 1 <= input_.size()) {
                char suffix_buf[suffix.size() + 1];
                prefix_buf[suffix.size()] = 0;
                for (int i = 0; i < suffix.size(); ++i) {
                    const auto unmapped_symbol = symbols[block_index * saved_data + prefix.size() + i];
                    const char mapped_symbol   = alphabet[unmapped_symbol];
                    suffix_buf[i]              = mapped_symbol;
                }

                ASSERT_EQ(suffix, std::string_view(suffix_buf, suffix.size()))
                    << "Incorrect suffix in block (" << block->start_index_ << ", " << block->end_index_ << ")";
            }
        }
    }
}

// This test checks if the mapping and alphabet fields are correct
TEST_P(CBlockTreeFixture, mapping_alphabet_check) {
    EXPECT_EQ(c_block_tree_rs_->mapping_.size(), characters_.size());
    EXPECT_EQ(c_block_tree_rs_->alphabet_->size(), characters_.size());
    for (int i = 0; i < c_block_tree_->alphabet_->size(); ++i) {
        EXPECT_EQ(c_block_tree_rs_->mapping_[(*c_block_tree_->alphabet_)[i]], i);
    }
}

// This test checks if the CBlockTree represents its input string and if the access method is correct
TEST_P(CBlockTreeFixture, access) {
    for (int i = 0; i < input_.size(); ++i) EXPECT_EQ(c_block_tree_->access(i), input_[i]);
}

/// Test substr method on substrings larger than the prefixes and suffixes stored in the blocks
TEST_P(CBlockTreeFixture, substr_larger) {
    constexpr size_t size = 50;
    char             buf[size + 1];

    for (size_t i = 0; i < input_.size() - size; i++) {
        std::ranges::fill(buf, 0);
        c_block_tree_->substr(buf, i, size);
        std::string_view actual   = {buf};
        std::string_view expected = {input_.begin() + i, input_.begin() + i + size};
        ASSERT_EQ(expected, actual) << "incorrect substring at index " << i;
    }
}

/// Test substr method on substrings smaller than the prefixes and suffixes stored in the blocks
TEST_P(CBlockTreeFixture, substr_smaller) {
    constexpr size_t size = 2;
    char             buf[size + 1];

    for (size_t i = 99000; i < input_.size() - size; i++) {
        std::ranges::fill(buf, 0);
        c_block_tree_->substr(buf, i, size);
        std::string_view actual   = {buf};
        std::string_view expected = {input_.begin() + i, input_.begin() + i + size};
        ASSERT_EQ(expected, actual) << "incorrect substring at index " << i;
    }
}

// This test checks the rank method for every character and position in the input
TEST_P(CBlockTreeFixture, rank) {
    for (const auto &[c, _] : characters_) {
        int r = 0;
        for (int i = 0; i < input_.size(); ++i) {
            if (input_[i] == c)
                ++r;
            EXPECT_EQ(c_block_tree_rs_->rank(c, i), r);
        }
    }
}

// This test checks the select method for every character and rank
TEST_P(CBlockTreeFixture, select) {
    for (const auto &[c, char_select_results] : characters_) {
        for (int j = 1; j <= char_select_results.size(); ++j)
            EXPECT_EQ(c_block_tree_rs_->select(c, j), char_select_results[j - 1]);
    }
}
