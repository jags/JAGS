#include <testlib.h>
#include <cfloat>
#include <cmath>

using std::sqrt;

JAGSFixture::JAGSFixture()
    : tol(sqrt(DBL_EPSILON)), 
      F({false}), T({true}),
      FF({false, false}), FT({false, true}),
      TF({true, false}),  TT({true, true}),
      FFF({false, false, false}),
      FFT({false, false, true}),
      FTF({false, true, false}),
      FTT({false, true, true}),
      TFF({true, false, false}),
      TFT({true, false, true}),
      TTF({true, true, false}),
      TTT({true, true, true})      
{
}

BoolIterator::BoolIterator(unsigned long n) :
    vector<bool>(n, false), atEnd(false) {}

void BoolIterator::next() {
    bool bump = true;
    for (unsigned long i = 0; i < size(); ++i) {
	if (bump) {
	    bool x = operator[](i); //current value
	    bump = x;
	    operator[](i) = !x;
	}
	else return;
    }
    if (bump) atEnd=true;
}
