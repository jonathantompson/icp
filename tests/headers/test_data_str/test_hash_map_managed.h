//
//  test_circular_buffer.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "icp/data_str/hash_map_managed.h"
#include "test_unit/test_unit.h"
#include "icp/data_str/hash_funcs.h"

#define TEST_HMM_START_SIZE 101  // A "bigish" prime
#define TEST_HMM_NUM_VALUES 52  // Enough to force at least one re-hash

using icp::data_str::HashMapManaged;
using icp::data_str::HashUInt;

// TEST 1: Create a hash table, insert items from 1:N
TEST(HashMapManaged, CreationAndInsertion) {
  HashMapManaged<uint32_t, uint32_t> ht(TEST_HMM_START_SIZE, &HashUInt);

  for (uint32_t i = 0; i < TEST_HMM_NUM_VALUES; i += 1) {
    uint32_t val = i*2;
    EXPECT_EQ(ht.insert(i, val), true);  // Insert key = i, value = i*2
  }
  // Now check that they're all there
  uint32_t val;
  for (uint32_t i = 0; i < TEST_HMM_NUM_VALUES; i += 1) {
    EXPECT_TRUE(ht.lookup(i, val));  // find key = i (value = i)
    EXPECT_EQ(i*2, val);
  }
  for (uint32_t i = TEST_HMM_NUM_VALUES; i < 2*TEST_HMM_NUM_VALUES; i += 1) {
    EXPECT_FALSE(ht.lookup(i, val));  // find key = i (value = i)
  }
  EXPECT_EQ(ht.count(), TEST_HMM_NUM_VALUES);
}

//// TEST 2: Create a hash table, insert items from 1:N and delete some
//TEST(HashMapManaged, CreationAndDeletion) {
//  HashMapManaged<uint32_t, uint32_t> ht(TEST_HMM_START_SIZE, &HashUInt);
//  EXPECT_FALSE(ht.remove(1));  // Nothing there: should fail
//  EXPECT_TRUE(ht.insert(0, 1));
//  EXPECT_FALSE(ht.remove(1));  // 0 is there, but 1 is not: should fail
//  EXPECT_TRUE(ht.insert(1, 0));
//  EXPECT_TRUE(ht.remove(1));  // should work this time
//  EXPECT_FALSE(ht.remove(1));  // Make sure the value is actually gone
//  uint32_t val;
//  EXPECT_TRUE(ht.lookup(0, val));  // key=0 should still be there
//  EXPECT_EQ(val, 1);  // key=0's value shouldn't have been corrupted
//}

// TEST 3: Duplicate insertion
TEST(HashMapManaged, DuplicateInsertion) {
  HashMapManaged<uint32_t, uint32_t> ht(TEST_HMM_START_SIZE, &HashUInt);
  EXPECT_TRUE(ht.insert(0, 1));
  EXPECT_FALSE(ht.insert(0, 0));  // Try twice with two different values
  EXPECT_FALSE(ht.insert(0, 1));  // (since value shouldn't matter)
}

// TEST 4: Create a hash table, clear it and make sure it cleared
TEST(HashMapManaged, CreationAndClear) {
  HashMapManaged<uint32_t, uint32_t> ht(TEST_HMM_START_SIZE, &HashUInt);
  for (uint32_t i = 0; i < TEST_HMM_NUM_VALUES; i += 1) {
    uint32_t val = i*2;
    EXPECT_EQ(ht.insert(i, val), true);  // Insert key = i, value = i*2
  }

  ht.clear();
  EXPECT_EQ(ht.count(), 0);

  // Now check that they're all not there  (TEST 1, ensures that they were)
  uint32_t val;
  for (uint32_t i = 0; i < 2*TEST_HMM_NUM_VALUES; i += 1) {
    EXPECT_FALSE(ht.lookup(i, val));  // find key = i (value = i)
  }
}

// TEST 5: String hashing
TEST(HashMapManaged, StringCreationAndInsertion) {
  HashMapManaged<std::string, uint32_t> ht(TEST_HM_START_SIZE, &HashString);
  std::stringstream ss;
    
  // Try a constant hash (check assembly to make sure hashing was constant)
  uint32_t cur_val;
  cur_val = 0;
  ht.insertPrehash(CONSTANT_HASH(ht.size(), "0"), "0", cur_val);

  for (uint32_t i = 1; i < TEST_HMM_NUM_VALUES; i += 1) {
    cur_val = i*2;
    ss.str("");
    ss << i;
    // Insert key = i, value = i*2
    EXPECT_EQ(ht.insert(ss.str().c_str(), cur_val), true);  
  }
  // Now check that they're all there
  uint32_t val;
  for (uint32_t i = 0; i < TEST_HMM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_TRUE(ht.lookup(ss.str().c_str(), val));  
    EXPECT_EQ(i*2, val);
  }
  for (uint32_t i = TEST_HMM_NUM_VALUES; i < 2*TEST_HMM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_FALSE(ht.lookup(ss.str().c_str(), val));  
  }
  EXPECT_EQ(ht.count(), TEST_HMM_NUM_VALUES);
}

// TEST 1Ptr: Create a hash table, insert items from 1:N
TEST(HashMapManagedPtr, CreationAndInsertion) {
  HashMapManaged<uint32_t, uint32_t*> ht(TEST_HMM_START_SIZE, &HashUInt);
  uint32_t* cur_val;

  for (uint32_t i = 0; i < TEST_HMM_NUM_VALUES; i += 1) {
    cur_val = new uint32_t[1];
    *cur_val = i*2;
    // Insert key = i, value = i*2
    EXPECT_EQ(ht.insert(i, cur_val), true);  
  }
  // Now check that they're all there
  uint32_t* val = NULL;
  for (uint32_t i = 0; i < TEST_HMM_NUM_VALUES; i += 1) {
    EXPECT_TRUE(ht.lookup(i, val));  // find key = i (value = i)
    EXPECT_EQ(i*2, *val);
  }
  for (uint32_t i = TEST_HMM_NUM_VALUES; i < 2*TEST_HMM_NUM_VALUES; i += 1) {
    EXPECT_FALSE(ht.lookup(i, val));  // find key = i (value = i)
  }
  EXPECT_EQ(ht.count(), TEST_HMM_NUM_VALUES);
}

//// TEST 2Ptr: Create a hash table, insert items from 1:N and delete some
//TEST(HashMapManagedPtr, CreationAndDeletion) {
//  HashMapManaged<uint32_t, uint32_t*> ht(TEST_HMM_START_SIZE, &HashUInt);
//  uint32_t* cur_val;
//  EXPECT_FALSE(ht.remove(1));  // Nothing there: should fail
//  cur_val = new uint32_t[1]; 
//  *cur_val = 1;
//  EXPECT_TRUE(ht.insert(0, cur_val));
//  EXPECT_FALSE(ht.remove(1));  // key 0 is there, but 1 is not: should fail
//  cur_val = new uint32_t[1]; 
//  *cur_val = 0;
//  EXPECT_TRUE(ht.insert(1, cur_val));
//  EXPECT_TRUE(ht.remove(1));  // should work this time
//  EXPECT_FALSE(ht.remove(1));  // Make sure the value is actually gone
//  EXPECT_TRUE(ht.lookup(0, cur_val));  // key=0 should still be there
//  EXPECT_EQ(*cur_val, 1);  // key=0's value shouldn't have been corrupted
//}

// TEST 3Ptr: Duplicate insertion
TEST(HashMapManagedPtr, DuplicateInsertion) {
  HashMapManaged<uint32_t, uint32_t*> ht(TEST_HMM_START_SIZE, &HashUInt);
  uint32_t* cur_val;
  cur_val = new uint32_t[1]; 
  *cur_val = 1;
  EXPECT_TRUE(ht.insert(0, cur_val));
  cur_val = new uint32_t[1]; 
  *cur_val = 0;
  EXPECT_FALSE(ht.insert(0, cur_val));  // Try twice with two different values
  delete cur_val;
  cur_val = new uint32_t[1]; 
  *cur_val = 1;
  EXPECT_FALSE(ht.insert(0, cur_val));  // (since value shouldn't matter)
  delete cur_val;
}

// TEST 4Ptr: Create a hash table, clear it and make sure it cleared
TEST(HashMapManagedPtr, CreationAndClear) {
  HashMapManaged<uint32_t, uint32_t*> ht(TEST_HMM_START_SIZE, &HashUInt);
  uint32_t* cur_val;
  for (uint32_t i = 0; i < TEST_HMM_NUM_VALUES; i += 1) {
    cur_val = new uint32_t[1];
    *cur_val = i*2;
    // Insert key = i, value = i*2
    EXPECT_EQ(ht.insert(i, cur_val), true);  
  }

  ht.clear();
  EXPECT_EQ(ht.count(), 0);

  // Now check that they're all not there  (TEST 1, ensures that they were)
  uint32_t* val = NULL;
  for (uint32_t i = 0; i < 2*TEST_HMM_NUM_VALUES; i += 1) {
    EXPECT_FALSE(ht.lookup(i, val));  // find key = i (value = i)
  }
}

// TEST 5Ptr: String hashing
TEST(HashMapManagedPtr, StringCreationAndInsertion) {
  HashMapManaged<std::string, uint32_t*> ht(TEST_HM_START_SIZE, &HashString);
  std::stringstream ss;
    
  // Try a constant hash (check assembly to make sure hashing was constant)
  uint32_t* cur_val;
  cur_val = new uint32_t[1];
  *cur_val = 0;
  ht.insertPrehash(CONSTANT_HASH(ht.size(), "0"), "0", cur_val);

  for (uint32_t i = 1; i < TEST_HMM_NUM_VALUES; i += 1) {
    cur_val = new uint32_t[1];
    *cur_val = i*2;
    ss.str("");
    ss << i;
    // Insert key = i, value = i*2
    EXPECT_EQ(ht.insert(ss.str().c_str(), cur_val), true);  
  }
  // Now check that they're all there
  uint32_t* val = NULL;
  for (uint32_t i = 0; i < TEST_HMM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_TRUE(ht.lookup(ss.str().c_str(), val));  
    EXPECT_EQ(i*2, *val);
  }
  for (uint32_t i = TEST_HMM_NUM_VALUES; i < 2*TEST_HMM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_FALSE(ht.lookup(ss.str().c_str(), val));  
  }
  EXPECT_EQ(ht.count(), TEST_HMM_NUM_VALUES);
}
