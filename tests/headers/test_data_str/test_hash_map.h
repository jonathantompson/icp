//
//  test_circular_buffer.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include <sstream>
#include "icp/data_str/hash_map.h"
#include "test_unit/test_unit.h"
#include "icp/data_str/hash_funcs.h"

#define TEST_HM_START_SIZE 101  // A "bigish" prime
#define TEST_HM_NUM_VALUES 52  // Enough to force at least one re-hash

using icp::data_str::HashMap;
using icp::data_str::HashUInt;
using icp::data_str::HashString;

// TEST 1: Create a hash table, insert items from 1:N
TEST(HashMap, CreationAndInsertion) {
  HashMap<uint32_t, uint32_t> ht(TEST_HM_START_SIZE, &HashUInt);

  for (uint32_t i = 0; i < TEST_HM_NUM_VALUES; i += 1) {
    uint32_t val = i*2;
    EXPECT_EQ(ht.insert(i, val), true);  // Insert key = i, value = i*2
  }
  // Now check that they're all there
  uint32_t val = 0;
  for (uint32_t i = 0; i < TEST_HM_NUM_VALUES; i += 1) {
    EXPECT_TRUE(ht.lookup(i, val));  // find key = i (value = i)
    EXPECT_EQ(i*2, val);
  }
  for (uint32_t i = TEST_HM_NUM_VALUES; i < 2*TEST_HM_NUM_VALUES; i += 1) {
    EXPECT_FALSE(ht.lookup(i, val));  // find key = i (value = i)
  }
  EXPECT_EQ(ht.count(), TEST_HM_NUM_VALUES);
}

// TEST 2: Create a hash table, insert items from 1:N and delete some
//TEST(HashMap, CreationAndDeletion) {
//  HashMap<uint32_t, uint32_t> ht(TEST_HM_START_SIZE, &HashUInt);
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
TEST(HashMap, DuplicateInsertion) {
  HashMap<uint32_t, uint32_t> ht(TEST_HM_START_SIZE, &HashUInt);
  EXPECT_TRUE(ht.insert(0, 1));
  EXPECT_FALSE(ht.insert(0, 0));  // Try twice with two different values
  EXPECT_FALSE(ht.insert(0, 1));  // (since value shouldn't matter)
}

// TEST 4: Create a hash table, clear it and make sure it cleared
TEST(HashMap, CreationAndClear) {
  HashMap<uint32_t, uint32_t> ht(TEST_HM_START_SIZE, &HashUInt);
  for (uint32_t i = 0; i < TEST_HM_START_SIZE; i += 1) {
    uint32_t val = i*2;
    EXPECT_EQ(ht.insert(i, val), true);  // Insert key = i, value = i*2
  }

  ht.clear();
  EXPECT_EQ(ht.count(), 0);

  // Now check that they're all not there  (TEST 1, ensures that they were)
  uint32_t val = 0;
  for (uint32_t i = 0; i < 2*TEST_HM_START_SIZE; i += 1) {
    EXPECT_FALSE(ht.lookup(i, val));  // find key = i (value = i)
  }
}

// TEST 5: String hashing
TEST(HashMap, StringCreationAndInsertion) {
  HashMap<std::string, uint32_t> ht(TEST_HM_START_SIZE, &HashString);
  std::stringstream ss;

  // Try a few dynamic vs constant hashes and make sure that they're the same
  EXPECT_EQ(CONSTANT_HASH(4294967295u, "0"), DYNAMIC_HASH(4294967295u, "0"));
  EXPECT_EQ(CONSTANT_HASH(4294967295u, "hello world"), 
    DYNAMIC_HASH(4294967295u, "hello world"));
  EXPECT_NEQ(CONSTANT_HASH(4294967295u, "hello world"), 
    DYNAMIC_HASH(4294967295u, "hello worle"));
    
  // Try a constant hash (check assembly to make sure hashing was constant)
  ht.insertPrehash(CONSTANT_HASH(ht.size(), "0"), "0", 0);

  for (uint32_t i = 1; i < TEST_HM_NUM_VALUES; i += 1) {
    uint32_t val = i*2;
    ss.str("");
    ss << i;
    // Insert key = "i", value = i*2
    EXPECT_EQ(ht.insert(ss.str().c_str(), val), true);  
  }
  // Now check that they're all there
  uint32_t val = 0;
  for (uint32_t i = 0; i < TEST_HM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_TRUE(ht.lookup(ss.str().c_str(), val));  
    EXPECT_EQ(i*2, val);
  }
  for (uint32_t i = TEST_HM_NUM_VALUES; i < 2*TEST_HM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_FALSE(ht.lookup(ss.str().c_str(), val));  
  }
  EXPECT_EQ(ht.count(), TEST_HM_NUM_VALUES);
}

// TEST 1Ptr: Create a hash table, insert items from 1:N
TEST(HashMapPtr, CreationAndInsertion) {
  HashMap<uint32_t, uint32_t*> ht(TEST_HM_START_SIZE, &HashUInt);
  uint32_t* vals = new uint32_t[TEST_HM_NUM_VALUES];
  uint32_t* cur_val;

  for (uint32_t i = 0; i < TEST_HM_NUM_VALUES; i += 1) {
    cur_val = &vals[i];
    *cur_val = i*2;
    // Insert key = i, value = i*2
    EXPECT_EQ(ht.insert(i, cur_val), true);  
  }
  // Now check that they're all there
  uint32_t* val = NULL;
  for (uint32_t i = 0; i < TEST_HM_NUM_VALUES; i += 1) {
    EXPECT_TRUE(ht.lookup(i, val));  // find key = i (value = i)
    EXPECT_EQ(i*2, *val);
  }
  for (uint32_t i = TEST_HM_NUM_VALUES; i < 2*TEST_HM_NUM_VALUES; i += 1) {
    EXPECT_FALSE(ht.lookup(i, val));  // find key = i (value = i)
  }
  EXPECT_EQ(ht.count(), TEST_HM_NUM_VALUES);
    
  // Now we must explicitly delete all the values since they will not be done
  // for us.
  delete[] vals;
}

//// TEST 2Ptr: Create a hash table, insert items from 1:N and delete some
//TEST(HashMapPtr, CreationAndDeletion) {
//  HashMap<uint32_t, uint32_t*> ht(TEST_HM_START_SIZE, &HashUInt);
//  uint32_t* vals = new uint32_t[2];
//  vals[0] = 1;
//  vals[1] = 0;
//  uint32_t* cur_val;
//  
//  EXPECT_FALSE(ht.remove(1));  // Nothing there: should fail
//  cur_val = &vals[0]; 
//  EXPECT_TRUE(ht.insert(0, cur_val));
//  EXPECT_FALSE(ht.remove(1));  // key 0 is there, but 1 is not: should fail
//  cur_val = &vals[1];
//  EXPECT_TRUE(ht.insert(1, cur_val));
//  EXPECT_TRUE(ht.remove(1));  // should work this time
//  EXPECT_FALSE(ht.remove(1));  // Make sure the value is actually gone
//  EXPECT_TRUE(ht.lookup(0, cur_val));  // key=0 should still be there
//  EXPECT_EQ(*cur_val, 1);  // key=0's value shouldn't have been corrupted
//  
//  // Now we must explicitly delete all the values since they will not be done
//  // for us.
//  delete[] vals;
//}

// TEST 3Ptr: Duplicate insertion
TEST(HashMapPtr, DuplicateInsertion) {
  HashMap<uint32_t, uint32_t*> ht(TEST_HM_START_SIZE, &HashUInt);
  uint32_t* cur_val0;
  uint32_t* cur_val1;
  cur_val0 = new uint32_t[1]; 
  *cur_val0 = 1;
  EXPECT_TRUE(ht.insert(0, cur_val0));
  cur_val1 = new uint32_t[1]; 
  *cur_val1 = 0;
  EXPECT_FALSE(ht.insert(0, cur_val1));  // Try twice with two different vals
  delete[] cur_val1;
  cur_val1 = new uint32_t[1]; 
  *cur_val1 = 1;
  EXPECT_FALSE(ht.insert(0, cur_val1));  // (since value shouldn't matter)
    
  // Now cleanup after ourselves
  delete cur_val0;    
  delete cur_val1;
}

// TEST 5Ptr: String hashing
TEST(HashMapPtr, StringCreationAndInsertion) {
  HashMap<std::string, uint32_t*> ht(TEST_HM_START_SIZE, &HashString);
  std::stringstream ss;

  uint32_t* vals = new uint32_t[TEST_HM_NUM_VALUES];
  uint32_t* cur_val;

  for (uint32_t i = 0; i < TEST_HM_NUM_VALUES; i += 1) {
    cur_val = &vals[i];
    *cur_val = i*2;
    ss.str("");
    ss << i;
    // Insert key = i, value = i*2
    EXPECT_EQ(ht.insert(ss.str().c_str(), cur_val), true);  
  }
  // Now check that they're all there
  uint32_t* val;
  for (uint32_t i = 0; i < TEST_HM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_TRUE(ht.lookup(ss.str().c_str(), val));  
    EXPECT_EQ(i*2, *val);
  }
  for (uint32_t i = TEST_HM_NUM_VALUES; i < 2*TEST_HM_NUM_VALUES; i += 1) {
    ss.str("");
    ss << i;
    // find key = i (value = i)
    EXPECT_FALSE(ht.lookup(ss.str().c_str(), val));  
  }
  EXPECT_EQ(ht.count(), TEST_HM_NUM_VALUES);
    
  // Now we must explicitly delete all the values since they will not be done
  // for us.
  delete[] vals;
}
