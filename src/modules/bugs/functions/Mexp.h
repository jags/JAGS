#ifndef MEXP_H_
#define MEXP_H_

#include <function/Function.h>

class Mexp : public Function
{
 public:
  Mexp();
  void evaluate (SArray &value,
		 std::vector <SArray const *> const &args) const;
  Index dim (std::vector <Index> const &dims) const;
  bool checkParameterDim (std::vector <Index> const &dims) const;
};

#endif /* MEXP_H_ */
