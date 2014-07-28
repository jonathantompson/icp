//
//  test_vec2_mat2x2.h
//
//  Created by Jonathan Tompson on 4/26/12.
//

#include "test_unit/test_unit.h"
#include "icp/math/math_types.h"
#include "icp/math/math_base.h"

using icp::math::IsPrime;
using icp::math::NextPrime;

TEST(MathBase, Primes) {
  EXPECT_FALSE(IsPrime(4));
  EXPECT_FALSE(IsPrime(6));
  EXPECT_FALSE(IsPrime(8));
  EXPECT_FALSE(IsPrime(9));
  EXPECT_FALSE(IsPrime(10));
  EXPECT_FALSE(IsPrime(12));
  EXPECT_FALSE(IsPrime(14));
  EXPECT_FALSE(IsPrime(15));
  EXPECT_FALSE(IsPrime(16));
  EXPECT_FALSE(IsPrime(18));
  EXPECT_FALSE(IsPrime(20));
  EXPECT_FALSE(IsPrime(111546435));  // 1*3*5*7*11*13*17*19*23

  EXPECT_TRUE(IsPrime(1));
  EXPECT_TRUE(IsPrime(2));
  EXPECT_TRUE(IsPrime(3));
  EXPECT_TRUE(IsPrime(5));
  EXPECT_TRUE(IsPrime(7));
  EXPECT_TRUE(IsPrime(11));
  EXPECT_TRUE(IsPrime(13));
  EXPECT_TRUE(IsPrime(17));
  EXPECT_TRUE(IsPrime(19));
  EXPECT_TRUE(IsPrime(15485863));

  EXPECT_EQ(NextPrime(8), 11);
  EXPECT_EQ(NextPrime(7178), 7187);
  EXPECT_EQ(NextPrime(7186), 7187);
}
