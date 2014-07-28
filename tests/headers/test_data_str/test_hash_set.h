//
//  test_circular_buffer.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "icp/data_str/hash_set.h"
#include "test_unit/test_unit.h"
#include "icp/data_str/hash_funcs.h"

#define TEST_HS_START_SIZE 101  // A "bigish" prime
#define TEST_HS_NUM_VALUES 2048

using icp::data_str::HashSet;
using icp::data_str::HashUInt;

// TEST 1: Create a hash table, insert items from 1:N
TEST(HashSet, CreationAndInsertion) {
  HashSet<uint32_t> ht(TEST_HS_START_SIZE, &HashUInt);

  for (int i = 0; i < TEST_HS_NUM_VALUES; i += 1) {
    EXPECT_EQ(ht.insert(i), true);  // Insert key = i, value = i*2
  }
  // Now check that they're all there
  for (uint32_t i = 0; i < TEST_HS_NUM_VALUES; i += 1) {
    EXPECT_TRUE(ht.lookup(i));  // find key = i (value = i)
  }
  for (uint32_t i = TEST_HS_NUM_VALUES; i < 2*TEST_HS_NUM_VALUES; i += 1) {
    EXPECT_FALSE(ht.lookup(i));  // find key = i (value = i)
  }
  EXPECT_EQ(ht.count(), TEST_HS_NUM_VALUES);
}

//// TEST 2: Create a hash table, insert items from 1:N and delete some
//TEST(HashSet, CreationAndDeletion) {
//  HashSet<uint32_t> ht(TEST_HS_START_SIZE, &HashUInt);
//  EXPECT_FALSE(ht.remove(1));  // Nothing there: should fail
//  EXPECT_TRUE(ht.insert(0));
//  EXPECT_FALSE(ht.remove(1));  // 0 is there, but 1 is not: should fail
//  EXPECT_TRUE(ht.insert(1));
//  EXPECT_TRUE(ht.remove(1));  // should work this time
//  EXPECT_FALSE(ht.remove(1));  // Make sure the value is actually gone
//  EXPECT_TRUE(ht.lookup(0));  // key=0 should still be there
//}

// TEST 3: Duplicate insertion
TEST(HashSet, DuplicateInsertion) {
  HashSet<uint32_t> ht(TEST_HS_START_SIZE, &HashUInt);
  EXPECT_TRUE(ht.insert(0));
  EXPECT_FALSE(ht.insert(0));  // Try twice with two different values
  EXPECT_FALSE(ht.insert(0));  // (since value shouldn't matter)
}

// TEST 4: Create a hash table, clear it and make sure it cleared
TEST(HashSet, CreationAndClear) {
  HashSet<uint32_t> ht(TEST_HM_START_SIZE, &HashUInt);
  for (uint32_t i = 0; i < TEST_HM_START_SIZE; i += 1) {
    EXPECT_EQ(ht.insert(i), true);  // Insert key = i, value = i*2
  }

  ht.clear();
  EXPECT_EQ(ht.count(), 0);

  // Now check that they're all not there  (TEST 1, ensures that they were)
  for (uint32_t i = 0; i < 2*TEST_HM_START_SIZE; i += 1) {
    EXPECT_FALSE(ht.lookup(i));  // find key = i (value = i)
  }
}

// TEST 5: String hashing
TEST(HashSet, StringCreationAndInsertion) {
  HashSet<std::string> ht(TEST_HM_START_SIZE, &HashString);
  std::stringstream ss;

  // Try a few dynamic vs constant hashes and make sure that they're the same
  EXPECT_EQ(CONSTANT_HASH(4294967295u, "0"), DYNAMIC_HASH(4294967295u, "0"));
  EXPECT_EQ(CONSTANT_HASH(4294967295u, "hello world"), 
    DYNAMIC_HASH(4294967295u, "hello world"));
  EXPECT_NEQ(CONSTANT_HASH(4294967295u, "hello world"), 
    DYNAMIC_HASH(4294967295u, "hello worle"));
    
  // Try a constant hash (check assembly to make sure hashing was constant)
  ht.insertPrehash(CONSTANT_HASH(ht.size(), "0"), "0");

  for (uint32_t i = 1; i < TEST_HM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // Insert key = "i", value = i*2
    EXPECT_EQ(ht.insert(ss.str().c_str()), true);  
  }
  // Now check that they're all there
  for (uint32_t i = 0; i < TEST_HM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_TRUE(ht.lookup(ss.str().c_str()));  
  }
  for (uint32_t i = TEST_HM_NUM_VALUES; i < 2*TEST_HM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_FALSE(ht.lookup(ss.str().c_str()));  
  }
  EXPECT_EQ(ht.count(), TEST_HM_NUM_VALUES);
}
