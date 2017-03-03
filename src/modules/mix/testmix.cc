#include "testmix.h"
#include "distributions/testmixdist.h"
#include <cppunit/extensions/HelperMacros.h>

void init_mix_test() {
    CPPUNIT_TEST_SUITE_REGISTRATION( MixDistTest );
}
