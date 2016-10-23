#include "testglm.h"
#include "samplers/testglmsamp.h"
#include "distributions/testglmdist.h"
#include <cppunit/extensions/HelperMacros.h>

void init_glm_test() {
    CPPUNIT_TEST_SUITE_REGISTRATION( GLMSampTest );
    CPPUNIT_TEST_SUITE_REGISTRATION( GLMDistTest );
}
