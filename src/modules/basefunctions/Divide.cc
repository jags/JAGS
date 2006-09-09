#include <config.h>
#include <sarray/SArray.h>
#include "Divide.h"

using std::vector;

Divide::Divide ():ScalarFunc ("/", 2)
{
}

double Divide::eval (vector <SArray const *> const &args) const
{
  return *args[0]->value() / *args[1]->value();
}

bool Divide::checkParameterValue (vector <SArray const *> const &args) const
{
  if (*(args[1]->value ()) == 0)
    return false;
  else
    return true;
}

bool Divide::isLinear(vector<bool> const &mask, vector<bool> const &fix) const
{
    if (mask[1])
	return false; //No reciprocal terms

    if (fix.empty()) {
	return true;
    }
    else {
	return !mask[0] || fix[1]; 
    }
}

/*
bool Divide::isScale(unsigned int index, vector<bool> const &fix) const
{
    if (index == 1)
	return false;

    if (fix.empty())
	return true;
    else
	return fix[1];
}
*/
