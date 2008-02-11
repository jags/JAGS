#include <config.h>
#include <sampler/SamplerFactory.h>

/* We must have a real destructor, even though it doesn't do anything,
 * otherwise the destructors for child classes never get called */

SamplerFactory::~SamplerFactory()
{
}
