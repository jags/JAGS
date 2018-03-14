#include "testbugs.h"
#include "functions/testbugsfun.h"
#include "distributions/testbugsdist.h"
#include <cppunit/extensions/HelperMacros.h>

void init_bugs_test() {
    CPPUNIT_TEST_SUITE_REGISTRATION( BugsFunTest );
    //CPPUNIT_TEST_SUITE_REGISTRATION( BugsDistTest );
}
