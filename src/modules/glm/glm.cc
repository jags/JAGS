#include <Module.h>

#include "samplers/LinearFactory.h"
#include "samplers/LinearGibbsFactory.h"
#include "samplers/IWLSFactory.h"
#include "samplers/HolmesHeldFactory.h"
#include "samplers/HolmesHeldBlockFactory.h"
#include "samplers/AlbertChibFactory.h"
#include "samplers/AlbertChibGibbsFactory.h"
#include "samplers/AMFactory.h"
#include "samplers/ConjugateFFactory.h"

#include <cholmod.h>

using std::vector;

cholmod_common *glm_wk = 0; /* Workspace for CHOLMOD */

namespace glm {
    
    class GLMModule : public Module {
    public:
	GLMModule();
	~GLMModule();
    };
    
    GLMModule::GLMModule() 
	: Module("glm")
    {
	glm_wk = new cholmod_common;
	cholmod_start(glm_wk);
/*
	glm_wk->final_ll = true; //Use LL' factorisation instead of LDL

	//For debuggin purposes we may choose not to reorder matrices
	//Use only on small problems

	glm_wk->nmethods = 1 ;
	glm_wk->method [0].ordering = CHOLMOD_NATURAL ;
	glm_wk->postorder = 0 ;
*/

	insert(new IWLSFactory);
	insert(new LinearGibbsFactory);
	insert(new LinearFactory);
 	insert(new AMFactory);
	insert(new AlbertChibGibbsFactory);
	insert(new AlbertChibFactory);
	insert(new HolmesHeldFactory);
	insert(new HolmesHeldBlockFactory);
	insert(new ConjugateFFactory);
    }
    
    GLMModule::~GLMModule() {
	
	vector<SamplerFactory*> const &svec = samplerFactories();
	for (unsigned int i = 0; i < svec.size(); ++i) {
	    delete svec[i];
	}

	cholmod_finish(glm_wk);
	delete glm_wk;
    }
}

glm::GLMModule _glm_module;
