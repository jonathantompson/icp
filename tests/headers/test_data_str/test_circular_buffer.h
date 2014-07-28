//
//  test_circular_buffer.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "icp/data_str/circular_buffer.h"
#include "test_unit/test_unit.h"

#define TEST_CB_SIZE 11
#define TEST_CB_NUM_REPEATS 2

using icp::data_str::CircularBuffer;

// TEST 1: Simple Write and Read, no overflow (but fill up buffer)
// Repeat the test in case previous writes and reads corrupt pointer state
// Also, repeating will check that the same number of read and writes 
// correctly empties the buffer.
TEST(CircularBuffer, ReadWrite) {
  CircularBuffer<int>* b = new CircularBuffer<int>(TEST_CB_SIZE);

  for (int j = 0; j < TEST_CB_NUM_REPEATS; j++) {
    for (int i = 0; i < TEST_CB_SIZE; i++)
      b->write(i);

    // Should get back exactly what we put in
    int value = 0;
    for (int i = 0; i < TEST_CB_SIZE; i++) {
      EXPECT_TRUE(b->read(value));
      EXPECT_EQ(value, i);
    }
  }
  delete b;
}

// TEST 2: Write overflow
// Repeat to check that TEST_BUFFER_SIZE reads will empty the buffer even if
// writes overflow
TEST(CircularBuffer, WriteOverflow) {
  CircularBuffer<int> b(TEST_CB_SIZE);

  for (int j = 0; j < TEST_CB_NUM_REPEATS; j++) {
    for (int i = 0; i < 2 * TEST_CB_SIZE; i++)
      b.write(i);

    // Should get back last TEST_BUFFER_SIZE values only (TEST_BUFFER_SIZE
    // ... 2*TEST_BUFFER_SIZE-1)
    int value = 0;
    for (int i = TEST_CB_SIZE; i < 2 * TEST_CB_SIZE; i++) {
      EXPECT_TRUE(b.read(value));
      EXPECT_EQ(value, i);
    }
  }
}

// TEST 3: Clear test
// Check that clearing the buffer works without corrupting buffer pointers.
TEST(CircularBuffer, FillAndClear) {
  CircularBuffer<int> b(TEST_CB_SIZE);

  for (int j = 0; j < TEST_CB_NUM_REPEATS; j++) {
    for (int i = 0; i < TEST_CB_SIZE; i++)
      b.write(i);
    b.clear();
    b.write(-1);
    int value;
    EXPECT_TRUE(b.read(value));
  }
}

// TEST 4: Read underflow
// Check that underflow case is handled correctly.
// Repeat the test to make sure under-flows don't corrupt the buffer state.
TEST(CircularBuffer, ReadUnderflow) {
  CircularBuffer<int> b(TEST_CB_SIZE);

  int value = 0;
  for (int j = 0; j < TEST_CB_NUM_REPEATS; j++) {
    b.write(1);
    EXPECT_TRUE(b.read(value));
    EXPECT_EQ(1, value);  // Should be 1
    for (int j = 0; j < TEST_CB_NUM_REPEATS; j++) {
      EXPECT_FALSE(b.read(value));
    }
  }
}
