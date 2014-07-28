//
//  test_vector.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "icp/data_str/min_heap.h"
#include "test_unit/test_unit.h"

#define TEST_HEAP_SIZE 541  // A "bigish" prime
#define TEST_HEAP_SIZE_BIG 8831

using icp::data_str::MinHeap;

TEST(MinHeap, SlowInsertionAndMinLookup) {
  MinHeap<int>* heap = new MinHeap<int>(TEST_HEAP_SIZE);
  // Insert the numbers in decreasing order (worst case) in O(nlog(n)) time
  for (int i = TEST_HEAP_SIZE-1; i >= 0; i--) {
    heap->insert(i);
    EXPECT_TRUE(heap->validate());
    // heap->print();
  }
    
  // Now extract the mins and make sure they are correct
  for (int i = 0; i < TEST_HEAP_SIZE; i++) {
    EXPECT_EQ(heap->size(), (uint32_t)(TEST_HEAP_SIZE-i));
    EXPECT_EQ(heap->removeMin(), i);
    EXPECT_TRUE(heap->validate());
    // heap->print();
  }
    
  delete heap;
}
  
TEST(MinHeap, FastInsertionAndMinLookup) {
  int* data = new int[TEST_HEAP_SIZE_BIG];
  // Insert the numbers in decreasing order
  for (int i = 0; i < TEST_HEAP_SIZE_BIG; i++) {
    data[TEST_HEAP_SIZE_BIG-i-1] = i;
  }
  // Build a heap from the bottom up in O(n) time
  MinHeap<int>* heap = new MinHeap<int>(data, TEST_HEAP_SIZE_BIG);
  EXPECT_TRUE(heap->validate());
  // heap->print();
    
  // Now extract the mins and make sure they are correct (just the first few)
  for (int i = 0; i < 10; i++) {
    EXPECT_EQ(heap->size(), (uint32_t)(TEST_HEAP_SIZE_BIG-i));
    EXPECT_EQ(heap->removeMin(), i);
  }
    
  delete[] data;
  delete heap;
}
  
