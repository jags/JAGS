#include <module/Module.h>

//#include "samplers/IWLSFactory.h"
#include "samplers/HolmesHeldFactory.h"
//#include "samplers/ConjugateFFactory.h"
#include "samplers/GLMGenericFactory.h"
#include "samplers/ScaledGammaFactory.h"
#include "samplers/ScaledWishartFactory.h"
#include "samplers/REScaledGammaFactory.h"
#include "samplers/REScaledWishartFactory.h"
#include "samplers/REGammaFactory.h"

#include "distributions/DScaledGamma.h"
#include "distributions/DScaledWishart.h"
#include "distributions/DOrderedLogit.h"
#include "distributions/DOrderedProbit.h"

#include <cholmod.h>

using std::vector;

cholmod_common *glm_wk = 0; /* Workspace for CHOLMOD */

namespace jags {
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

        //Force use of simplicial factorization. Supernodal factorizations
	//have a completely different data structure, although held in
	//the same object.
        glm_wk->supernodal = CHOLMOD_SIMPLICIAL;

/*	
	//Force use of LL' factorisation instead of LDL
	//glm_wk->final_ll = true; 

	//For debuggin purposes we may choose not to reorder matrices
	//Use only on small problems

	glm_wk->nmethods = 1 ;
	glm_wk->method [0].ordering = CHOLMOD_NATURAL ;
	glm_wk->postorder = 0 ;
*/

	//insert(new IWLSFactory);
	insert(new GLMGenericFactory);
	insert(new HolmesHeldFactory);
	//insert(new ConjugateFFactory);
	insert(new ScaledGammaFactory);
	insert(new ScaledWishartFactory);

	/*
	insert(new REScaledGammaFactory);
	insert(new REScaledWishartFactory);
	insert(new REGammaFactory);
	*/
	
	insert(new DScaledGamma);
	insert(new DScaledWishart);
	insert(new DOrderedLogit);
	insert(new DOrderedProbit);
    }
    
    GLMModule::~GLMModule() {
	
	vector<SamplerFactory*> const &svec = samplerFactories();
	for (unsigned int i = 0; i < svec.size(); ++i) {
	    delete svec[i];
	}

	cholmod_finish(glm_wk);
	delete glm_wk;
    }

}}

jags::glm::GLMModule _glm_module;
