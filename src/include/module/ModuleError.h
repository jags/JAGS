#ifndef MODULE_ERROR_H_
#define MODULE_ERROR_H_

#include <string>

namespace jags {

class Node;
class Distribution;
class Function;

    [[noreturn]] void throwRuntimeError(std::string const &message);

    [[noreturn]]
    void throwLogicError(std::string const &message);

    [[noreturn]]
    void throwNodeError(Node const *node, std::string const &message);

    [[noreturn]]
    void throwDistError(Distribution const *dist, std::string const &message);

    [[noreturn]]
    void throwFuncError(Function const *func, std::string const &message);

} /* namespace jags */

#endif /* MODULE_ERROR_H_ */
