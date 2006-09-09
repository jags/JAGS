#include <Module.h>

#include "Add.h"
#include "Divide.h"
#include "Equals.h"
#include "Multiply.h"
#include "Neg.h"
#include "Pow.h"
#include "Subtract.h"

using std::vector;

class BaseFuncModule : public Module {
public:
    BaseFuncModule();
    ~BaseFuncModule();
};

BaseFuncModule::BaseFuncModule() {

    insert(new Add);
    insert(new Divide);
    insert(new Equals);
    insert(new Multiply);
    insert(new Neg);
    insert(new Pow);
    insert(new Subtract);
}

BaseFuncModule::~BaseFuncModule() {

    vector<Function*> const &fvec = functions();
    for (unsigned int i = 0; i < fvec.size(); ++i) {
	delete fvec[i];
    }
}

BaseFuncModule _base_func_module;

