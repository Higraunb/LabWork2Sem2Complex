#include "HComplex.h"

#include <gtest.h>

TEST(TComplex, can_create_Complex)
{
  ASSERT_NO_THROW(TComplex<int> a(1, 1));
}

TEST(TComplex, can_create_Complex_copy)
{
  TComplex<int> a(1, 2), d(a);
  EXPECT_NEAR(1, a == d, 0.00001);
}

TEST(TComplex, can_create_Complex_initializer)
{
  TComplex<int> a(1, 2);
  EXPECT_EQ(1, a.GetRe());
  EXPECT_EQ(2, a.GetIm());
}

TEST(TComplex, can_Complex_Set_Get)
{
  TComplex<int> a;
  a.SetIm(2);
  a.SetRe(4);
  EXPECT_EQ(4, a.GetRe());
  EXPECT_EQ(2, a.GetIm());
}
TEST(TComplex, checking_the_operation_functionality_plus_equal)
{
  TComplex<double> a(1.0, 1.0), b(2.0,2.0),
    c(a.GetRe() + b.GetRe(), a.GetIm() + b.GetIm());
  a += b;
  EXPECT_NEAR(1, c == a, 0.00001);
}

TEST(TComplex, checking_the_operation_functionality_plus)
{
  TComplex<double> a(1.0, 1.0), b(2.0, 2.0),
    c(a.GetRe() + b.GetRe(), a.GetIm() + b.GetIm()), d;
  d = b + a;
  EXPECT_NEAR(1, c == d, 0.00001);
}

TEST(TComplex, checking_the_operation_functionality_equal)
{
  TComplex<double> a(1.0, 1.0), b(1.0, 1.0);
  bool e = (a == b);
  EXPECT_TRUE(e);
}

TEST(TComplex, checking_the_operation_functionality_not_equal)
{
  TComplex<double> a(1.0, 1.0), b(1.0, 1.0);
  bool e = (a != b);
  EXPECT_FALSE(e);
}

TEST(TComplex, checking_the_operation_functionality_minus_equal)
{
  TComplex<double> a(1.0, 1.0), b(2.0, 2.0),
    c(a.GetRe() - b.GetRe(), a.GetIm() - b.GetIm());
  a -= b;
  EXPECT_NEAR(1, c == a, 0.00001);
}

TEST(TComplex, checking_the_operation_functionality_minus)
{
  TComplex<double> a(1.0, 1.0), b(2.0, 2.0),
    c(a.GetRe() - b.GetRe(), a.GetIm() - b.GetIm()), d;
  d = a - b;
  EXPECT_NEAR(1, c == d, 0.00001);
}

//TEST(TComplex, checking_the_operation_functionality_equating)
//{
//  TComplex<double> b(2.0, 2.0), d(3.0, 4.0);
//  d = b;
//  EXPECT_EQ(d, b);
//}

TEST(TComplex, checking_the_operation_functionality_multiplication)
{
  TComplex<double> a(1.0, 1.0), b(2.0, 2.0),
    c(a.GetRe() * b.GetRe() - a.GetIm() * b.GetIm(),
      a.GetRe() * b.GetIm() + a.GetIm() * b.GetRe()), d;
  d = a * b;
  EXPECT_NEAR(1, c == d, 0.00001);
}

TEST(TComplex, checking_the_operation_functionality_multiplication_equal)
{
  TComplex<double> a(1.0, 1.0), b(2.0, 2.0),
    c(a.GetRe() * b.GetRe() - a.GetIm() * b.GetIm(),
      a.GetRe() * b.GetIm() + a.GetIm() * b.GetRe()), d;
  a *= b;
  EXPECT_NEAR(1, c == a, 0.00001);
}

TEST(TComplex, checking_the_operation_functionality_division_equal)
{
  TComplex<double> a(1.0, 1.0), b(2.0, 2.0),
    c((a.GetRe() * b.GetRe() + a.GetIm() * b.GetIm()) / (b.GetRe() * b.GetRe() + b.GetIm() * b.GetIm()),
      (a.GetIm() * b.GetRe() - a.GetRe() * b.GetIm()) / (b.GetIm() * b.GetIm() + b.GetIm() * b.GetIm()));
  a /= b;
  EXPECT_NEAR(1, c == a, 0.00001);
}

TEST(TComplex, checking_the_operation_functionality_division)
{
  TComplex<double> a(1.0, 1.0), b(2.0, 2.0),
    c((a.GetRe() * b.GetRe() + a.GetIm() * b.GetIm()) / (b.GetRe() * b.GetRe() + b.GetIm() * b.GetIm()),
      (a.GetIm() * b.GetRe() - a.GetRe() * b.GetIm()) / (b.GetIm() * b.GetIm() + b.GetIm() * b.GetIm())), d;
  d = a / b;
  EXPECT_NEAR(1, c == d, 0.00001);
}


